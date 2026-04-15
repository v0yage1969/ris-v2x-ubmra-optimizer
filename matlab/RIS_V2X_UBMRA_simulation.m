%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RIS_V2X_UBMRA_simulation.m
%
% RIS辅助车联网中继系统：基于上界引导最小RIS激活策略（UB-MRA）
%
% 系统模型：
%   - 源节点 S (BS) + DF中继 R + RIS面板 + 4辆车编队
%   - 路径损耗：3GPP TR 38.901 UMi Street Canyon LOS
%   - 小尺度衰落：Rician (K=5)
%   - RIS相位：理想对齐
%
% 对比方案：
%   1. Fixed-256  : 始终激活256个RIS元件
%   2. Fixed-128  : 始终激活128个RIS元件
%   3. Proposed-UBMRA : 自适应选择最小激活规模
%
% 输出：3张图 + 3个CSV + 控制台汇总
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% =========================================================
%% 1. 系统参数
%% =========================================================

% 载波与带宽
fc_GHz = 3.5;                   % 载波频率 (GHz)
BW     = 20e6;                  % 系统带宽 (Hz)
NF_dB  = 7;                     % 接收机噪声系数 (dB)

% 三维节点坐标 [x, y, z]，单位 m
S_pos   = [0,   0,  15];        % 源节点 (基站)
R_pos   = [150, 0,   8];        % DF中继
RIS_pos = [75,  10,  8];        % RIS面板中心

% 编队参数
N_cars = 4;                     % 车辆数量
d_car  = 20;                    % 纵向车间距 (m)
z_car  = 1.5;                   % 车辆天线高度 (m)

% 头车位置扫描
x_head_vec = 60:10:260;         % 头车x坐标扫描范围 (m)
n_pos = length(x_head_vec);

% 候选RIS激活规模集合
N_cand = [64, 128, 192, 256];
n_cand = length(N_cand);

% Rician衰落参数
K_rice = 5;

% 蒙特卡洛仿真次数
N_MC = 500;

% UB-MRA阈值
R_UB_thresh     = 7.5;          % 队尾速率上界预选门限 (bits/s/Hz)
J_thresh        = 0.97;         % Jain公平性指数门限
R_outage_thresh = 6;            % 队尾中断判据 (bits/s/Hz)

% 功率参数（用于归一化能效计算，为工程假设值）
P_s    = 1;                     % 源发射功率 (W)
P_r    = 1;                     % 中继发射功率 (W)
P_c    = 1.5;                   % 电路功耗 (W)
P_elem = 0.02;                  % 每个RIS元件功耗 (W)
eta    = 1;                     % RIS反射幅度（假设理想）

% 噪声功率
k_B    = 1.38e-23;              % 玻尔兹曼常数 (J/K)
T0     = 290;                   % 参考温度 (K)
NF_lin = 10^(NF_dB / 10);
sigma2 = k_B * T0 * BW * NF_lin;  % 热噪声功率 ≈ 4e-13 W

% 发射SNR
rho_s = P_s / sigma2;           % 源端SNR
rho_r = P_r / sigma2;           % 中继SNR

%% =========================================================
%% 2. 路径损耗函数（3GPP TR 38.901 UMi Street Canyon LOS）
%%    PL [dB] = 32.4 + 21*log10(d3D) + 20*log10(fc_GHz)
%% =========================================================
PL_dB  = @(d3D) 32.4 + 21 * log10(max(d3D, 1)) + 20 * log10(fc_GHz);
PL_lin = @(d3D) 10 .^ (-PL_dB(d3D) / 10);

%% =========================================================
%% 3. 辅助函数
%% =========================================================

% Rician复信道系数向量（单位大尺度增益，通过乘以sqrt(PL_lin)引入大尺度）
% 输出 1×N 复向量，满足 E[|h|^2] = 1
gen_rician_ch = @(N) ...
    sqrt(K_rice / (K_rice + 1)) * exp(1j * 2 * pi * rand(1, N)) + ...
    sqrt(1 / (K_rice + 1)) * (randn(1, N) + 1j * randn(1, N)) / sqrt(2);

% Jain公平性指数
jain_idx = @(r) sum(r)^2 / (numel(r) * sum(r .^ 2) + eps);

%% =========================================================
%% 4. 预计算固定链路大尺度参数
%% =========================================================
d_SR  = norm(S_pos - R_pos);
pl_SR = PL_lin(d_SR);

d_SI  = norm(S_pos - RIS_pos);
pl_SI = PL_lin(d_SI);

fprintf('===== UB-MRA 仿真参数 =====\n');
fprintf('载波频率: %.1f GHz | 带宽: %.0f MHz | 噪声功率: %.2e W\n', ...
    fc_GHz, BW/1e6, sigma2);
fprintf('rho_s = %.2e (%.1f dB)\n', rho_s, 10*log10(rho_s));
fprintf('d_SR = %.1f m | d_SI = %.1f m\n', d_SR, d_SI);
fprintf('===========================\n\n');

%% =========================================================
%% 5. 主仿真循环
%% =========================================================

% 方案设置：1=Fixed-256, 2=Fixed-128, 3=Proposed-UBMRA
n_schemes    = 3;
scheme_names = {'Fixed-256', 'Fixed-128', 'Proposed-UBMRA'};

% 位置级结果矩阵
N_act_ubmra = zeros(1, n_pos);     % UB-MRA在每个位置选取的N_act
avg_rate    = zeros(n_schemes, n_pos);
tail_p5     = zeros(n_schemes, n_pos);
fairness    = zeros(n_schemes, n_pos);
outage_pct  = zeros(n_schemes, n_pos);

fprintf('仿真进度（共%d个头车位置）：\n', n_pos);

for pi = 1:n_pos
    x_head = x_head_vec(pi);

    % ---------- 编队车辆位置（头=1，尾=N_cars）----------
    x_cars  = x_head - (0 : N_cars - 1) * d_car;
    car_pos = [x_cars', zeros(N_cars, 1), ones(N_cars, 1) * z_car];

    % ---------- 各车大尺度路径增益 ----------
    pl_ID = zeros(1, N_cars);
    pl_RD = zeros(1, N_cars);
    for k = 1:N_cars
        d_ID_k   = norm(RIS_pos - car_pos(k, :));
        d_RD_k   = norm(R_pos   - car_pos(k, :));
        pl_ID(k) = PL_lin(d_ID_k);
        pl_RD(k) = PL_lin(d_RD_k);
    end

    % ======================================================
    % Step 2-3: UB-MRA预选最小激活规模
    %
    % 使用Jensen不等式上界：
    %   R_UB = log2(1 + E[gamma])
    %   E[gamma_RIS_UB] <= rho_s * eta^2 * N^2 * pl_SI * pl_ID(k)
    %
    % 选取满足 R_tail_UB >= R_UB_thresh 且 J_UB >= J_thresh 的最小N
    % ======================================================
    N_act_sel  = 256;           % 默认满配（兜底）
    gamma_SR_ls = rho_s * pl_SR;

    for ni = 1:n_cand
        N_try = N_cand(ni);

        % 各车RIS速率上界
        gamma_RIS_UB_vec = rho_s * eta^2 * (N_try^2) * pl_SI * pl_ID;
        R_RIS_UB_vec     = log2(1 + gamma_RIS_UB_vec);

        % 各车DF速率（大尺度确定值，无小尺度波动）
        R_DF_UB_vec = zeros(1, N_cars);
        for k = 1:N_cars
            gamma_RD_ls     = rho_r * pl_RD(k);
            R_DF_UB_vec(k)  = 0.5 * log2(1 + min(gamma_SR_ls, gamma_RD_ls));
        end

        R_sum_UB_vec   = R_RIS_UB_vec + R_DF_UB_vec;
        R_tail_UB_val  = R_sum_UB_vec(end);        % 尾车速率上界
        J_UB_val       = jain_idx(R_sum_UB_vec);   % 编队公平性上界

        if R_tail_UB_val >= R_UB_thresh && J_UB_val >= J_thresh
            N_act_sel = N_try;
            break;
        end
    end
    N_act_ubmra(pi) = N_act_sel;

    % ======================================================
    % Step 4: 三种方案的蒙特卡洛仿真
    % ======================================================
    N_act_per_scheme = [256, 128, N_act_sel];

    for si = 1:n_schemes
        N_act  = N_act_per_scheme(si);
        rates_MC = zeros(N_MC, N_cars);

        for mc = 1:N_MC
            % S -> R 信道（Rician）
            h_SR_mc    = sqrt(pl_SR) * gen_rician_ch(1);
            gamma_SR_mc = rho_s * abs(h_SR_mc)^2;

            % S -> RIS 信道（N_act个元件）
            h_SI_arr = sqrt(pl_SI) * gen_rician_ch(N_act);  % 1 × N_act

            for k = 1:N_cars
                % RIS -> 车k 信道
                h_ID_arr = sqrt(pl_ID(k)) * gen_rician_ch(N_act);

                % R -> 车k 信道
                h_RD_mc    = sqrt(pl_RD(k)) * gen_rician_ch(1);
                gamma_RD_mc = rho_r * abs(h_RD_mc)^2;

                % RIS速率（理想相位对齐：各反射元件相位使所有路径同相叠加）
                ris_eff_gain  = sum(abs(h_SI_arr) .* abs(h_ID_arr));
                gamma_RIS_mc  = rho_s * eta^2 * ris_eff_gain^2;
                R_RIS_mc      = log2(1 + gamma_RIS_mc);

                % DF中继速率（半双工，两时隙）
                R_DF_mc = 0.5 * log2(1 + min(gamma_SR_mc, gamma_RD_mc));

                rates_MC(mc, k) = R_DF_mc + R_RIS_mc;
            end
        end

        % ---------- 统计指标 ----------
        % 平均编队速率（对所有MC轮次和所有车辆取均值）
        avg_rate(si, pi)   = mean(mean(rates_MC, 2));

        % 队尾5%分位速率（尾车 = 第N_cars辆）
        tail_p5(si, pi)    = prctile(rates_MC(:, end), 5);

        % 队尾中断率
        outage_pct(si, pi) = mean(rates_MC(:, end) < R_outage_thresh) * 100;

        % Jain公平性指数（每次MC后计算，取均值）
        J_mc = arrayfun(@(m) jain_idx(rates_MC(m, :)), 1:N_MC);
        fairness(si, pi)   = mean(J_mc);
    end

    fprintf('  [%2d/%d] x_head=%3dm | N_act_UB=%3d | AvgRate: F256=%.2f F128=%.2f UBMRA=%.2f\n', ...
        pi, n_pos, x_head, N_act_sel, ...
        avg_rate(1, pi), avg_rate(2, pi), avg_rate(3, pi));
end

%% =========================================================
%% 6. 汇总统计与归一化能效
%% =========================================================
avg_N_act = [256, 128, mean(N_act_ubmra)];

% 归一化能效 EE = AvgRate / P_total，以Fixed-256为基准归一化
P_total_vec = P_s + P_r + P_c + avg_N_act * P_elem;
EE_raw      = mean(avg_rate, 2)' ./ P_total_vec;
EE_norm     = EE_raw / EE_raw(1);

fprintf('\n========= 汇总结果 =========\n');
fprintf('%-22s %9s %9s %9s %9s %9s\n', ...
    '方案', 'Avg N', 'AvgRate', 'TailP5', 'Outage%', 'Norm.EE');
for si = 1:n_schemes
    fprintf('%-22s %9.2f %9.4f %9.4f %9.4f %9.4f\n', ...
        scheme_names{si}, avg_N_act(si), ...
        mean(avg_rate(si, :)), mean(tail_p5(si, :)), ...
        mean(outage_pct(si, :)), EE_norm(si));
end

% 与Fixed-256对比
delta_N   = (avg_N_act(3) - avg_N_act(1)) / avg_N_act(1) * 100;
delta_rate = (mean(avg_rate(3,:)) - mean(avg_rate(1,:))) / mean(avg_rate(1,:)) * 100;
delta_tail = (mean(tail_p5(3,:)) - mean(tail_p5(1,:))) / mean(tail_p5(1,:)) * 100;
delta_ee   = (EE_norm(3) - EE_norm(1)) / EE_norm(1) * 100;
fprintf('\n[UBMRA vs Fixed-256] N: %+.1f%% | AvgRate: %+.1f%% | TailP5: %+.1f%% | Norm.EE: %+.1f%%\n', ...
    delta_N, delta_rate, delta_tail, delta_ee);

%% =========================================================
%% 7. 导出CSV文件
%% =========================================================

% --- 位置级结果 ---
T_pos = table(...
    x_head_vec', N_act_ubmra', ...
    avg_rate(1,:)', avg_rate(2,:)', avg_rate(3,:)', ...
    tail_p5(1,:)', tail_p5(2,:)', tail_p5(3,:)', ...
    fairness(1,:)', fairness(2,:)', fairness(3,:)', ...
    outage_pct(1,:)', outage_pct(2,:)', outage_pct(3,:)', ...
    'VariableNames', {...
    'x_head', 'N_act_UBMRA', ...
    'AvgRate_F256', 'AvgRate_F128', 'AvgRate_UBMRA', ...
    'TailP5_F256',  'TailP5_F128',  'TailP5_UBMRA', ...
    'Fairness_F256','Fairness_F128', 'Fairness_UBMRA', ...
    'Outage_F256',  'Outage_F128',  'Outage_UBMRA'});
writetable(T_pos, 'RIS_V2X_UBMRA_position_results.csv');

% --- 汇总结果 ---
T_sum = table(scheme_names', avg_N_act', ...
    mean(avg_rate, 2), mean(tail_p5, 2), ...
    mean(fairness, 2), mean(outage_pct, 2), EE_norm', ...
    'VariableNames', {'Scheme','Avg_N_act','Avg_Rate','Tail_P5', ...
    'Fairness','Outage_pct','Norm_EE'});
writetable(T_sum, 'RIS_V2X_UBMRA_summary_results.csv');

% --- 配对t检验 ---
[~, p1, ~, s1] = ttest(avg_rate(3,:)', avg_rate(1,:)');  % UBMRA vs F256, avgrate
[~, p2, ~, s2] = ttest(avg_rate(3,:)', avg_rate(2,:)');  % UBMRA vs F128, avgrate
[~, p3, ~, s3] = ttest(outage_pct(3,:)', outage_pct(2,:)'); % UBMRA vs F128, outage
[~, p4, ~, s4] = ttest(tail_p5(3,:)', tail_p5(2,:)');   % UBMRA vs F128, tailP5

T_stats = table(...
    {'Avg rate: Proposed vs Fixed-256'; ...
     'Avg rate: Proposed vs Fixed-128'; ...
     'Outage: Proposed vs Fixed-128';   ...
     'Tail P5: Proposed vs Fixed-128'}, ...
    [s1.tstat; s2.tstat; s3.tstat; s4.tstat], ...
    [p1; p2; p3; p4], ...
    'VariableNames', {'Test', 't_stat', 'p_value'});
writetable(T_stats, 'RIS_V2X_UBMRA_evaluation_stats.csv');

fprintf('\nCSV文件已保存：\n');
fprintf('  RIS_V2X_UBMRA_position_results.csv\n');
fprintf('  RIS_V2X_UBMRA_summary_results.csv\n');
fprintf('  RIS_V2X_UBMRA_evaluation_stats.csv\n');

%% =========================================================
%% 8. 绘图
%% =========================================================
c = lines(3);           % 三个方案颜色
mk = {'o', 's', 'd'};  % 标记形状
lw = 2;  ms = 7;

% ---- 图1: 平均编队速率 ----
fig1 = figure('Name', 'Avg Platoon Rate', 'NumberTitle', 'off', ...
    'Position', [80, 420, 680, 430]);
hold on;
for si = 1:n_schemes
    plot(x_head_vec, avg_rate(si, :), ...
        'Color', c(si,:), 'Marker', mk{si}, ...
        'LineWidth', lw, 'MarkerSize', ms, ...
        'DisplayName', scheme_names{si});
end
xlabel('Lead vehicle position x_{head} (m)', 'FontSize', 11);
ylabel('Average platoon rate (bits/s/Hz)',    'FontSize', 11);
title('Average platoon rate versus lead-vehicle position', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 10);
grid on; box on;
saveas(fig1, 'RIS_V2X_UBMRA_avg_rate.png');

% ---- 图2: 队尾5%分位速率 ----
fig2 = figure('Name', 'Tail P5 Rate', 'NumberTitle', 'off', ...
    'Position', [780, 420, 680, 430]);
hold on;
for si = 1:n_schemes
    plot(x_head_vec, tail_p5(si, :), ...
        'Color', c(si,:), 'Marker', mk{si}, ...
        'LineWidth', lw, 'MarkerSize', ms, ...
        'DisplayName', scheme_names{si});
end
yline(R_outage_thresh, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Outage threshold (6 b/s/Hz)');
xlabel('Lead vehicle position x_{head} (m)', 'FontSize', 11);
ylabel('Tail-user 5th percentile rate (bits/s/Hz)', 'FontSize', 11);
title('Tail-user P5 rate versus lead-vehicle position', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 10);
grid on; box on;
saveas(fig2, 'RIS_V2X_UBMRA_tail_p5.png');

% ---- 图3: UB-MRA选择的激活规模 ----
fig3 = figure('Name', 'UB-MRA Schedule', 'NumberTitle', 'off', ...
    'Position', [80, -50, 680, 380]);
stairs(x_head_vec, N_act_ubmra, 'b-', 'LineWidth', lw);
hold on;
plot(x_head_vec, N_act_ubmra, 'bo', 'MarkerSize', ms, 'MarkerFaceColor', 'b');
xlabel('Lead vehicle position x_{head} (m)', 'FontSize', 11);
ylabel('Selected active RIS elements',       'FontSize', 11);
title('UB-MRA selected active RIS elements', 'FontSize', 12);
yticks(N_cand); ylim([0, 290]);
grid on; box on;
saveas(fig3, 'RIS_V2X_UBMRA_schedule.png');

% ---- 图4: 队尾中断率 ----
fig4 = figure('Name', 'Outage Probability', 'NumberTitle', 'off', ...
    'Position', [780, -50, 680, 380]);
hold on;
for si = 1:n_schemes
    plot(x_head_vec, outage_pct(si, :), ...
        'Color', c(si,:), 'Marker', mk{si}, ...
        'LineWidth', lw, 'MarkerSize', ms, ...
        'DisplayName', scheme_names{si});
end
xlabel('Lead vehicle position x_{head} (m)',      'FontSize', 11);
ylabel('Tail-user outage probability (%)',          'FontSize', 11);
title('Outage probability (R_{tail} < 6 bits/s/Hz)', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 10);
grid on; box on;
saveas(fig4, 'RIS_V2X_UBMRA_outage.png');

fprintf('图形已保存：avg_rate / tail_p5 / schedule / outage (PNG)\n');

%% =========================================================
%% 9. 输出统计检验摘要
%% =========================================================
fprintf('\n========= 配对t检验结果 =========\n');
fprintf('%-40s %10s %10s\n', '检验项目', 't统计量', 'p值');
for i = 1:height(T_stats)
    fprintf('%-40s %10.4f %10.4f\n', ...
        T_stats.Test{i}, T_stats.t_stat(i), T_stats.p_value(i));
end

fprintf('\n===== 仿真完成 =====\n');
fprintf('答辩关键数字：\n');
fprintf('  平均激活单元数   : %.1f（较Fixed-256降低 %.1f%%）\n', ...
    avg_N_act(3), -delta_N);
fprintf('  全位置平均速率   : %.4f bits/s/Hz\n', mean(avg_rate(3,:)));
fprintf('  队尾5%%分位速率   : %.4f bits/s/Hz\n', mean(tail_p5(3,:)));
fprintf('  队尾中断率       : %.2f%%\n', mean(outage_pct(3,:)));
fprintf('  归一化能效       : %.4f（较Fixed-256提升 %.1f%%）\n', ...
    EE_norm(3), delta_ee);
