%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB 脚本：参数化分析 (A方案：系统级速率仿真)
%
% 目标：复现专利中的 图3、图4、图5
%
% *** 变更：***
% *** 1. 仿真 1 (图2) 横坐标 x_U 范围扩大到 [-100, 800] ***
% *** 2. 仿真 2 (图3) 横坐标 SNR 范围扩大 ***
% *** 3. 仿真 3 (图4) 已重写，复现专利图5 (平衡点 vs 中继-RIS间距) ***
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; 
clc;
close all;

%% 1. 定义系统参数 (基于专利)
eta = 1;          % RIS 反射幅度 (假设为1)
S_pos_base = [250, 833];  % 源 (S) 基础位置
I_pos_base = [0, 400];    % RIS (I) 基础位置
R_pos_base = [500, 400];  % 中继 (R) 基础位置

%% 2. 定义编队 (Platoon) 参数
K = 4;            % 编队中的车辆数 (例如 4 辆)
d_platoon = 20;   % 车辆间距 (例如 20 米)

%% 3. 假设的信道参数
rho_lin_base = 2.7e11;   % 基础：线性信噪比 (P/sigma^2)
alpha = 3.0;      % 假设：路径损耗指数 alpha

%% 4. 仿真一：分析 RIS 单元数 N 的影响 (复现专利图3)
% -------------------------------------------------------------------------
fprintf('正在运行仿真 1：分析 N 的影响 (图3) ...\n');

N_values = [250, 500, 750, 1000];
% x_U_lead = 0:10:500; % (旧范围)
x_U_lead = -100:10:800; % (新范围：扩大横坐标)
num_N_values = length(N_values);
num_positions = length(x_U_lead);

R_Total_avg_vs_N = zeros(num_positions, num_N_values);

% 使用基础坐标
d_SI = norm(S_pos_base - I_pos_base);
d_SR = norm(S_pos_base - R_pos_base);
E_beta_SI = d_SI ^ (-alpha);
E_beta_SR = d_SR ^ (-alpha);

for n_idx = 1:num_N_values
    N = N_values(n_idx);
    ris_gain_term = (N * (16 + (N-1)*pi^2)) / 16;
    R_Total_temp = zeros(num_positions, K);

    for i = 1:num_positions  % 遍历位置
        for k = 1:K  % 遍历车辆
            D_k_x_pos = x_U_lead(i) - (k-1) * d_platoon;
            D_k_pos = [D_k_x_pos, 0];
            d_ID_k = norm(I_pos_base - D_k_pos);
            d_RD_k = norm(R_pos_base - D_k_pos);
            E_beta_ID_k = d_ID_k ^ (-alpha);
            E_beta_RD_k = d_RD_k ^ (-alpha);

            min_Relay_Path_Gain_k = min(E_beta_SR, E_beta_RD_k);
            R_Relay = 0.5 * log2(1 + rho_lin_base * min_Relay_Path_Gain_k);
            E_gamma_RIS_k = rho_lin_base * (eta^2) * E_beta_SI * E_beta_ID_k * ris_gain_term;
            R_RIS = log2(1 + E_gamma_RIS_k);
            R_Total_temp(i, k) = R_Relay + R_RIS;
        end
    end
    R_Total_avg_vs_N(:, n_idx) = mean(R_Total_temp, 2);
end 

% 绘制图 2 (复现专利图3)
figure(2);
plot(x_U_lead, R_Total_avg_vs_N, 'LineWidth', 2);
title('图2：编队平均速率 vs. 头车位置 (按 N 区分)');
xlabel('头车位置 x_U [m]');
ylabel('编队平均速率 [bits/s/Hz]');
legend(arrayfun(@(n) sprintf('N = %d', n), N_values, 'UniformOutput', false), 'Location', 'best');
grid on;
xlim([x_U_lead(1) x_U_lead(end)]); % (新范围)
ylim([0, max(R_Total_avg_vs_N, [], 'all')*1.1]);
fprintf('仿真 1 (图3复现) 完成。\n');

%% 5. 仿真二：分析 SNR (rho) 的影响 (复现专利图4)
% -------------------------------------------------------------------------
fprintf('正在运行仿真 2：分析 SNR 的影响 (图4) ...\n');

% rho_lin_values = logspace(10.5, 12.5, 20); % (旧范围)
rho_lin_values = logspace(10, 13, 20); % (新范围：扩大横坐标)
N_snr_values = [250, 500, 750, 1000];
x_U_fixed = 50;
num_rho_values = length(rho_lin_values);
num_N_snr_values = length(N_snr_values);

R_Total_avg_vs_P = zeros(num_rho_values, num_N_snr_values);

% 使用基础坐标
d_SI_snr = norm(S_pos_base - I_pos_base);
d_SR_snr = norm(S_pos_base - R_pos_base);
E_beta_SI_snr = d_SI_snr ^ (-alpha);
E_beta_SR_snr = d_SR_snr ^ (-alpha);

for n_idx = 1:num_N_snr_values
    N = N_snr_values(n_idx);
    ris_gain_term = (N * (16 + (N-1)*pi^2)) / 16;
    
    R_Total_temp = zeros(num_rho_values, K);
    
    E_beta_ID_k_fixed = zeros(1, K);
    E_beta_RD_k_fixed = zeros(1, K);
    for k = 1:K
        D_k_x_pos = x_U_fixed - (k-1) * d_platoon;
        D_k_pos = [D_k_x_pos, 0];
        d_ID_k = norm(I_pos_base - D_k_pos);
        d_RD_k = norm(R_pos_base - D_k_pos);
        E_beta_ID_k_fixed(k) = d_ID_k ^ (-alpha);
        E_beta_RD_k_fixed(k) = d_RD_k ^ (-alpha);
    end
    
    for r_idx = 1:num_rho_values
        rho_lin = rho_lin_values(r_idx);
        for k = 1:K
            min_Relay_Path_Gain_k = min(E_beta_SR_snr, E_beta_RD_k_fixed(k));
            R_Relay = 0.5 * log2(1 + rho_lin * min_Relay_Path_Gain_k);
            E_gamma_RIS_k = rho_lin * (eta^2) * E_beta_SI_snr * E_beta_ID_k_fixed(k) * ris_gain_term;
            R_RIS = log2(1 + E_gamma_RIS_k);
            R_Total_temp(r_idx, k) = R_Relay + R_RIS;
        end
    end
    R_Total_avg_vs_P(:, n_idx) = mean(R_Total_temp, 2);
end 

% 绘制图 3 (复现专利图4)
figure(3);
rho_db_axis = 10*log10(rho_lin_values / 1e10);
plot(rho_db_axis, R_Total_avg_vs_P, '-s', 'LineWidth', 2);
title(sprintf('图3：编队平均速率 vs. 信噪比 (复现专利图4, x_U = %d m)', x_U_fixed));
xlabel('相对信噪比 (dB)');
ylabel('编队平均速率 [bits/s/Hz]');
legend(arrayfun(@(n) sprintf('N = %d', n), N_snr_values, 'UniformOutput', false), 'Location', 'best');
grid on;
ylim([0, max(R_Total_avg_vs_P, [], 'all')*1.1]);
fprintf('仿真 2 (图4复现) 完成。\n');

%% 6. 仿真三：复现专利图5 (平衡点 vs 中继-RIS间距)
% -------------------------------------------------------------------------
% 目标：找到 R_RIS = R_Relay 时的 x_U 位置
% X轴: d_IR (中继-RIS间距)
% Y轴: x_U (平衡点位置)
% 曲线: 不同的 N
% -------------------------------------------------------------------------
fprintf('正在运行仿真 3：复现专利图5 (平衡点) ...\n');

% 定义要遍历的参数 (基于专利 [0128] 描述)
d_IR_values = 0:20:500; % X轴: 中继-RIS间距, 0到500m
N_values_fig5 = [250, 500, 750, 1000]; % 4条曲线
num_d_IR = length(d_IR_values);
num_N_fig5 = length(N_values_fig5);

% 定义仿真坐标系 (基于专利 [0128])
d_S_line = 100; % 源 S 到 I-R 线的垂直距离
d_D_line = 400; % 用户 D 到 I-R 线的垂直距离 (基于主模型)

x_U_search = -200:2:1000; % 搜索平衡点的 x_U 范围 - 扩大范围
num_x_search = length(x_U_search);

% 初始化存储结果
balance_points_fig5 = zeros(num_d_IR, num_N_fig5);

% 循环遍历 N (每条曲线)
for n_idx = 1:num_N_fig5
    N = N_values_fig5(n_idx);
    ris_gain_term = (N * (16 + (N-1)*pi^2)) / 16;

    % 循环遍历 d_IR (X轴)
    for d_idx = 1:num_d_IR
        d_IR = d_IR_values(d_idx);
        
        % 1. 定义当前循环的坐标
        I_pos = [0, 0];
        R_pos = [d_IR, 0];
        S_pos = [d_IR/2, d_S_line]; % S 位于 I 和 R 中间，垂直距离 d_S_line
        
        % 计算 S->I 和 S->R 增益 (固定)
        d_SI = norm(S_pos - I_pos);
        d_SR = norm(S_pos - R_pos);
        E_beta_SI = d_SI ^ (-alpha);
        E_beta_SR = d_SR ^ (-alpha);
        
        % 2. 搜索 x_U 来找到平衡点
        R_RIS_temp   = zeros(num_x_search, 1);
        R_Relay_temp = zeros(num_x_search, 1);

        for i = 1:num_x_search
            x_U = x_U_search(i);
            D_pos = [x_U, -d_D_line]; % D 位于平行线上
            
            d_ID = norm(I_pos - D_pos);
            d_RD = norm(R_pos - D_pos);
            E_beta_ID = d_ID ^ (-alpha);
            E_beta_RD = d_RD ^ (-alpha);
            
            min_Relay_Path_Gain = min(E_beta_SR, E_beta_RD);
            R_Relay_temp(i) = 0.5 * log2(1 + rho_lin_base * min_Relay_Path_Gain);

            E_gamma_RIS = rho_lin_base * (eta^2) * E_beta_SI * E_beta_ID * ris_gain_term;
            R_RIS_temp(i) = log2(1 + E_gamma_RIS);
        end
        
        % 3. 找到 R_RIS 和 R_Relay 之间的交叉点 (使用符号变化检测)
        rate_diff = R_RIS_temp - R_Relay_temp;
        
        % 识别符号变化
        sign_change_mask = rate_diff(1:end-1) .* rate_diff(2:end) <= 0;
        cross_idx = find(sign_change_mask, 1, 'first');
        
        if ~isempty(cross_idx)
            % 线性插值交点位置
            x1 = x_U_search(cross_idx);
            x2 = x_U_search(cross_idx + 1);
            y1 = rate_diff(cross_idx);
            y2 = rate_diff(cross_idx + 1);
            
            if abs(y2 - y1) < 1e-9
                X_E_val = (x1 + x2) / 2;
            else
                X_E_val = x1 - y1 * (x2 - x1) / (y2 - y1);
            end
            
            % 边界检查
            balance_points_fig5(d_idx, n_idx) = max(x_U_search(1), min(x_U_search(end), X_E_val));
        else
            % 无交叉，选择最接近的点
            [min_diff, min_idx] = min(abs(rate_diff));
            if min_diff < 1e-3
                balance_points_fig5(d_idx, n_idx) = x_U_search(min_idx);
            else
                balance_points_fig5(d_idx, n_idx) = NaN; % 无平衡点
            end
        end
    end
end % 结束 N 循环

% 绘制图 4 (复现专利图5)
figure(4);
plot(d_IR_values, balance_points_fig5, '-s', 'LineWidth', 2);
title('图4：平衡点 vs. 中继-RIS间距 (复现专利图5)');
xlabel('中继-RIS 间距 [m]');
ylabel('平衡点位置 x_U [m]');
legend(arrayfun(@(n) sprintf('N = %d', n), N_values_fig5, 'UniformOutput', false), 'Location', 'best');
grid on;
ylim([0, max(balance_points_fig5, [], 'all', 'omitnan')*1.1]);
fprintf('仿真 3 (图5复现) 完成。\n');