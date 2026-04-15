%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB 脚本：复现专利 图 5 (平衡点 vs. 距离, 对比 N)
%
% 对应专利：图 5
% 仿真内容：平衡位置 (X_E) vs. 中继-RIS 间距 (d_IR)
% 变化参数：RIS 单元数 (N)
%
% 注意：
% 1. "平衡位置" (X_E 或 L) 定义为 R_RIS = R_Relay 的用户位置。
% 2. "中继-RIS 间距" (d_IR) 是图 5 的 x 轴，我们将通过改变 R 的 x 坐标
%    来实现，同时保持 I 的 x 坐标为 0。
% 3. 此仿真计算量较大。
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; 
clc;
close all;
%% 1. 定义系统参数 (基于专利)
% -------------------------------------------------------------------------
% *** (修改) 调整 N 值的范围 - 选择能产生交叉点的范围 ***
N_values = [2500, 3500, 4500, 5500]; % RIS 单元数量 (数组) - 精细调整
eta = 1;          % RIS 反射幅度 (假设为1)
% 节点坐标 (基于 "数值结果与分析")
% *** (修改) 调整S位置使其更靠近中继区域 ***
S_pos = [200, 500];  % 源 (S) 位置 - 降低高度，向右移
I_pos = [0, 400];    % RIS (I) 位置
% R_pos 将在循环中变化
%% 1.5. 定义 5G 和通信物理参数 (新)
% -------------------------------------------------------------------------
fc = 3.5e9;          % 载波频率 (Hz) - 3.5 GHz (Sub-6 GHz)
bw = 100e6;         % 系统带宽 (Hz) - 100 MHz
c = 299792458;       % 光速 (m/s)
% --- 源 (S) 和 用户 (D) 参数 ---
txPower_S_dBm = 30;   % 源 (S) 发射功率 (dBm)
txGain_S_dBi = 20;    % 源 (S) 发射天线增益 (dBi)
rxGain_D_dBi = 10;    % 用户 (D) 接收天线增益 (dBi)
% --- (新) 中继 (R) 节点参数 ---
% *** (修改) 大幅提高中继功率和增益以产生交叉点 ***
txPower_R_dBm = 40;   % 中继 (R) 发射功率 (dBm) - 提升至40dBm
txGain_R_dBi = 20;    % 中继 (R) 发射增益 (dBi) - 提升至20
rxGain_R_dBi = 20;    % 中继 (R) 接收增益 (dBi) - 提升至20
rxNoiseFigure_dB = 7; % 接收机噪声系数 (dB) (假设 D 和 R 相同)
% 计算热噪声功率 (dBm)
thermalNoise_dBm_per_Hz = -174; 
noisePower_dBm = thermalNoise_dBm_per_Hz + 10*log10(bw) + rxNoiseFigure_dB;
fprintf('5G 物理参数:\n');
fprintf('  载波频率: %.1f GHz\n', fc/1e9);
fprintf('  源 (S) 发射功率: %.0f dBm\n', txPower_S_dBm);
fprintf('  中继 (R) 发射功率: %.0f dBm\n', txPower_R_dBm);
fprintf('  接收机底噪: %.2f dBm\n', noisePower_dBm);
%% 2. 定义仿真范围
% -------------------------------------------------------------------------
K = 1; % 仅考虑单用户 (D)
% 搜索平衡点的用户位置范围
x_U_search_range = -200:2:1000; % [m] - 大幅扩大范围以覆盖所有可能的平衡点
num_positions = length(x_U_search_range);
% 变化的中继-RIS 间距 (专利图 5 的 x 轴)
Relay_to_RIS_dist_range = 0:10:500; % [m]
num_distances = length(Relay_to_RIS_dist_range);
%% 4. 计算固定距离和增益
% -------------------------------------------------------------------------
% S -> I 链路是固定的
d_SI = norm(S_pos - I_pos); % S -> I
pl_SI_dB = calculatePathLoss_dB(d_SI, fc, c);
fprintf('S->I 距离: %.1f m, 路径损耗: %.2f dB\n', d_SI, pl_SI_dB);
% I -> D_k 链路取决于 x_U，在内层循环计算
R_RIS_vec = zeros(num_positions, 1);
R_Relay_vec = zeros(num_positions, 1);
% 初始化存储 (使用 NaN 便于后续检测)
X_E_results = nan(length(N_values), num_distances);
%% 5. 循环计算平衡点 (X_E)
% -------------------------------------------------------------------------
for n_idx = 1:length(N_values)
    N = N_values(n_idx);
    
    % RIS 增益项
    ris_gain_term = (N * (16 + (N-1)*pi^2)) / 16;
    ris_gain_dB = 10*log10(ris_gain_term);
    fprintf('正在计算 N = %d (增益 %.2f dB)...\n', N, ris_gain_dB);
    
    % --- 预计算 R_RIS 速率 (不随 R_pos 变化) ---
    for i = 1:num_positions
        D_k_pos = [x_U_search_range(i), 0];
        d_ID_k = norm(I_pos - D_k_pos);
        pl_IDk_dB = calculatePathLoss_dB(d_ID_k, fc, c);
        
        % 链路预算: S -> I -> D
        Prx_RISk_dBm = txPower_S_dBm + txGain_S_dBi + rxGain_D_dBi - pl_SI_dB - pl_IDk_dB + ris_gain_dB;
        snr_RIS_dB = Prx_RISk_dBm - noisePower_dBm;
        R_RIS_vec(i) = log2(1 + 10^(snr_RIS_dB/10));
    end
    
    % --- 遍历 R_pos 来查找平衡点 ---
    for dist_idx = 1:num_distances
        d_IR = Relay_to_RIS_dist_range(dist_idx);
        % 定义中继位置 R_pos
        R_pos = [I_pos(1) + d_IR, I_pos(2)]; % [d_IR, 400]
        
        % S -> R 链路取决于 R_pos
        d_SR = norm(S_pos - R_pos);
        pl_SR_dB = calculatePathLoss_dB(d_SR, fc, c);
        
        % --- 计算 R_Relay 速率 vs. x_U ---
        for i = 1:num_positions
            D_k_pos = [x_U_search_range(i), 0];
            d_RD_k = norm(R_pos - D_k_pos);
            pl_RDk_dB = calculatePathLoss_dB(d_RD_k, fc, c);
            
            % 链路预算: S -> R
            Prx_SR_dBm = txPower_S_dBm + txGain_S_dBi + rxGain_R_dBi - pl_SR_dB;
            % 链路预算: R -> D
            Prx_RDk_dBm = txPower_R_dBm + txGain_R_dBi + rxGain_D_dBi - pl_RDk_dB;
            
            snr_SR_dB = Prx_SR_dBm - noisePower_dBm;
            snr_RDk_dB = Prx_RDk_dBm - noisePower_dBm;
            snr_Relay_lin = min( 10^(snr_SR_dB/10), 10^(snr_RDk_dB/10) );
            R_Relay_vec(i) = 0.5 * log2(1 + snr_Relay_lin);
        end
        
        % --- 查找 R_RIS_vec 和 R_Relay_vec 的交叉点 ---
        rate_diff = R_RIS_vec - R_Relay_vec;

        % 识别任意符号变化 (包含等于零的情形)
        sign_change_mask = rate_diff(1:end-1) .* rate_diff(2:end) <= 0;
        cross_idx = find(sign_change_mask, 1, 'first');

        if ~isempty(cross_idx)
            % 线性插值交点位置，提高精度
            x1 = x_U_search_range(cross_idx);
            x2 = x_U_search_range(cross_idx + 1);
            y1 = rate_diff(cross_idx);
            y2 = rate_diff(cross_idx + 1);
            % 避免除零：若 y1 == y2，直接取区间中点
            if abs(y2 - y1) < 1e-9
                X_E_val = (x1 + x2) / 2;
            else
                X_E_val = x1 - y1 * (x2 - x1) / (y2 - y1);
            end
            
            % 检查插值结果是否在合理范围内
            if X_E_val < x_U_search_range(1) || X_E_val > x_U_search_range(end)
                % 超出搜索范围，使用最接近的边界点
                X_E_results(n_idx, dist_idx) = max(x_U_search_range(1), min(x_U_search_range(end), X_E_val));
                if dist_idx <= 3
                    fprintf('  d_IR=%.0f m -> X_E=%.1f m (插值超界，已限制)\n', d_IR, X_E_results(n_idx, dist_idx));
                end
            else
                X_E_results(n_idx, dist_idx) = X_E_val;
            end
            
            % 调试输出：显示部分找到的平衡点
            if dist_idx <= 3 || mod(dist_idx, 10) == 0
                fprintf('  d_IR=%.0f m -> X_E=%.1f m (R_RIS=%.2f, R_Relay=%.2f)\n', ...
                        d_IR, X_E_results(n_idx, dist_idx), ...
                        R_RIS_vec(cross_idx), R_Relay_vec(cross_idx));
            end
        else
            % 无符号变化，选择最接近的点作为近似平衡
            [min_diff, min_idx] = min(abs(rate_diff));
            if min_diff < 1e-3 % 速率差足够小，认为近似相等
                X_E_results(n_idx, dist_idx) = x_U_search_range(min_idx);
            else
                % 无交叉且差值大 -> 留空 (NaN)，稍后绘图时跳过
                X_E_results(n_idx, dist_idx) = nan;
            end
        end
    end
    fprintf('N = %d 计算完成。\n', N);
    
    % 调试输出：显示该 N 值下找到的平衡点统计
    valid_points = X_E_results(n_idx, :);
    num_valid = sum(~isnan(valid_points));
    if num_valid > 0
        fprintf('  找到 %d 个有效平衡点，范围: [%.1f, %.1f] m\n', ...
                num_valid, min(valid_points, [], 'omitnan'), max(valid_points, [], 'omitnan'));
    else
        fprintf('  警告：未找到任何有效平衡点！\n');
    end
end

% 总体统计
fprintf('\n=== 总体结果统计 ===\n');
for n_idx = 1:length(N_values)
    valid_count = sum(~isnan(X_E_results(n_idx, :)));
    fprintf('N = %d: %d/%d 个距离点有平衡解\n', N_values(n_idx), valid_count, num_distances);
end
fprintf('\n');

%% 6. 绘图 (复现 图 5)
% -------------------------------------------------------------------------
figure('Name', '图 5 (复现): 平衡点 vs. 中继-RIS 间距 (对比 N)');
hold on;
markers = {'-o', '-s', '-d', '-^'};
for n_idx = 1:length(N_values)
    plot(Relay_to_RIS_dist_range, X_E_results(n_idx, :), ...
         markers{mod(n_idx-1, length(markers)) + 1}, ...
         'LineWidth', 1.5, 'MarkerSize', 5, ...
         'DisplayName', sprintf('N = %d', N_values(n_idx)));
end
title('平衡位置 (X_E) vs. 中继-RIS 水平间距 (对比 N 值)');
xlabel('中继-RIS 水平间距 (d_{IR}) [m]');
ylabel('平衡位置 X_E [m]');
grid on;
box on;
legend('Location', 'best');
xlim([Relay_to_RIS_dist_range(1) Relay_to_RIS_dist_range(end)]);
ylim_max = max(X_E_results, [], 'all', 'omitnan') * 1.1;
% *** (新) 修复：如果最大值仍为 0 或 NaN，则设置一个默认的 Y 轴高度 ***
if isnan(ylim_max) || ylim_max <= 0
    ylim_max = 100; % 默认 Y 轴最大值 (例如 100m)
end
ylim([0, ylim_max]);
hold off;
%% 辅助函数：计算 FSPL
function pl_dB = calculatePathLoss_dB(d_m, fc_Hz, c)
    if d_m < 1.0
        d_m = 1.0; 
    end
    lambda = c / fc_Hz;
    pl_dB = 20 * log10(4 * pi * d_m / lambda);
end