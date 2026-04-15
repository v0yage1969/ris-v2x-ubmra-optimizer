%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB 脚本：多车辆编队行驶模型 (5G 物理链路预算版)
%
% 扩展自 专利 CN 114584587 B 模型
%
% ========================= 变更说明 ====================================
%
% 此脚本已从原始的 "d^-alpha" 抽象模型修改为使用 5G 和通信工具箱
% (Communication Toolbox) 概念的物理链路预算模型。
%
% 1. [替换] 抽象信噪比 (rho_lin) 被替换为
%    - 5G 载波频率 (fc), 带宽 (bw)
%    - 发射功率 (txPower_dBm), 天线增益 (txGain_dBi, rxGain_dBi)
%    - 接收机噪声系数 (rxNoiseFigure_dB)
% 2. [替换] 路径损耗模型 (d^-alpha) 被替换为
%    - 自由空间路径损耗 (FSPL)，使用标准物理公式计算。
% 3. [保留] RIS 增益项 (ris_gain_term) 的公式被保留，但被转换为 dB
%    单位 (ris_gain_dB) 并集成到链路预算中。
% 4. [更新] 速率计算现在基于香non公式 C = B * log2(1 + SNR)，
%    其中 SNR 是通过完整的链路预算 (TxPwr + Gains - Losses - Noise) 计算得出的。
%    为了与原图保持一致 (bits/s/Hz)，我们计算 C/B。
%
% =========================================================================
clear; 
clc;
close all;
%% 1. 定义系统参数 (基于专利)
% -------------------------------------------------------------------------
% *** 修改：增加 N 以提高 RIS 增益 ***
N = 10000;         % RIS 反射单元数量 (N=10000)
eta = 1;          % RIS 反射幅度 (假设为1)
% 节点坐标 (基于 "数值结果与分析")
S_pos = [250, 833];  % 源 (S) 位置
I_pos = [0, 400];    % RIS (I) 位置
R_pos = [500, 400];  % 中继 (R) 位置
%% 1.5. 定义 5G 和通信物理参数 (新)
% -------------------------------------------------------------------------
% *** 修改：降低 fc 以减少路径损耗 ***
fc = 3.5e9;          % 载波频率 (Hz) - 3.5 GHz (Sub-6 GHz)
bw = 100e6;         % 系统带宽 (Hz) - 100 MHz
c = 299792458;       % 光速 (m/s)
txPower_dBm = 30;   % 发射功率 (dBm) - (例如 1W)
% *** 修改：增加 G_tx 假设为基站 ***
txGain_dBi = 20;    % 发射天线增益 (dBi)
rxGain_dBi = 10;    % 接收天线增益 (dBi)
rxNoiseFigure_dB = 7; % 接收机噪声系数 (dB)
% 计算热噪声功率 (dBm)
% k = 1.38e-23 (玻尔兹曼常数)
% T = 290 (开尔文)
% thermalNoise_dBm = 10*log10(k*T*1000) = -174 dBm/Hz
thermalNoise_dBm_per_Hz = -174; 
noisePower_dBm = thermalNoise_dBm_per_Hz + 10*log10(bw) + rxNoiseFigure_dB;
fprintf('5G 物理参数 (已修改):\n');
fprintf('  载波频率: %.1f GHz\n', fc/1e9);
fprintf('  RIS 单元数: %d\n', N);
fprintf('  发射增益: %.0f dBi\n', txGain_dBi);
fprintf('  带宽: %.0f MHz\n', bw/1e6);
fprintf('  发射功率: %.0f dBm\n', txPower_dBm);
fprintf('  接收机底噪: %.2f dBm\n', noisePower_dBm);
%% 2. 定义编队 (Platoon) 参数
% -------------------------------------------------------------------------
K = 4;            % 编队中的车辆数 (例如 4 辆)
d_platoon = 20;   % 车辆间距 (例如 20 米)
% 头车 (D_1) 沿x轴移动
x_U_lead = 0:10:500; % 头车位置 x_U [m]
num_positions = length(x_U_lead);
%% 3. 假设的信道参数 (已弃用)
% -------------------------------------------------------------------------
% rho_lin = 2.7e11;   % (已弃用) - 被物理链路预算取代
% alpha = 3.0;      % (已弃用) - 被 FSPL 模型取代
%% 4. 计算距离和信道增益 (已更新)
% -------------------------------------------------------------------------
% 计算固定距离
d_SI = norm(S_pos - I_pos); % S -> I
d_SR = norm(S_pos - R_pos); % S -> R
% 计算固定路径损耗 (dB)
pl_SI_dB = calculatePathLoss_dB(d_SI, fc, c);
pl_SR_dB = calculatePathLoss_dB(d_SR, fc, c);
% RIS 增益项 (来自专利)
ris_gain_term = (N * (16 + (N-1)*pi^2)) / 16;
ris_gain_dB = 10*log10(ris_gain_term); % 将线性增益转换为dB
fprintf('RIS 理论增益: %.2f dB\n', ris_gain_dB);
% 初始化存储速率的矩阵 (行: 位置, 列: 车辆)
R_Relay_platoon = zeros(num_positions, K);
R_RIS_platoon   = zeros(num_positions, K);
R_Total_platoon = zeros(num_positions, K);
%% 5. 循环计算每辆车的速率 (使用物理链路预算更新)
% -------------------------------------------------------------------------
for i = 1:num_positions  % 遍历头车的每个位置
    
    for k = 1:K  % 遍历编队中的每辆车
        
        % 5.1. 计算第k辆车的当前x坐标
        D_k_x_pos = x_U_lead(i) - (k-1) * d_platoon;
        D_k_y_pos = 0; % 假设都在x轴上行驶
        D_k_pos = [D_k_x_pos, D_k_y_pos];
        
        % 5.2. 计算第k辆车的距离
        d_ID_k = norm(I_pos - D_k_pos); % I -> D_k
        d_RD_k = norm(R_pos - D_k_pos); % R -> D_k
        
        % 5.3. 计算第k辆车的路径损耗 (dB)
        pl_IDk_dB = calculatePathLoss_dB(d_ID_k, fc, c);
        pl_RDk_dB = calculatePathLoss_dB(d_RD_k, fc, c);
        
        % 5.4. 计算第k辆车的速率 (bits/s/Hz)
        
        % --- 场景 1: "仅中继" 速率 (半双工) ---
        % 计算 S->R 和 R->D_k 两跳的接收功率和SNR
        Prx_SR_dBm = txPower_dBm + txGain_dBi + rxGain_dBi - pl_SR_dB;
        Prx_RDk_dBm = txPower_dBm + txGain_dBi + rxGain_dBi - pl_RDk_dB;
        
        snr_SR_dB = Prx_SR_dBm - noisePower_dBm;
        snr_RDk_dB = Prx_RDk_dBm - noisePower_dBm;
        
        % 中继链路的 SNR 取决于两跳中的 "瓶颈"
        snr_Relay_lin = min( 10^(snr_SR_dB/10), 10^(snr_RDk_dB/10) );
        
        % 速率 (bits/s/Hz)，0.5 因子来自半双工中继
        R_Relay_platoon(i, k) = 0.5 * log2(1 + snr_Relay_lin);
        
        % --- 场景 2: "仅RIS" 速率 ---
        % 链路预算: P_tx + G_tx - PL_SI - PL_IDk + G_RIS + G_rx
        % 注意: RIS 增益被加到 S-I-D 级联路径损耗上
        Prx_RISk_dBm = txPower_dBm + txGain_dBi + rxGain_dBi - pl_SI_dB - pl_IDk_dB + ris_gain_dB;
        
        snr_RIS_dB = Prx_RISk_dBm - noisePower_dBm;
        snr_RIS_lin = 10^(snr_RIS_dB/10);
        
        % 速率 (bits/s/Hz)，假设为全双工 (无 0.5 因子)
        R_RIS_platoon(i, k) = log2(1 + snr_RIS_lin);
        
        % --- 场景 3: "RIS + 中继" 速率 ---
        % 假设两者正交 (例如不同时隙或频段)，速率相加
        R_Total_platoon(i, k) = R_Relay_platoon(i, k) + R_RIS_platoon(i, k);
    end
end
%% 6. 绘图 (综合对比图)
% -------------------------------------------------------------------------
% *** 修改：恢复绘图样式为线条，因为现在曲线已分离 ***
figure('Name', '5G 物理链路预算仿真 (Sub-6GHz & 高增益)');
hold on;
% 设置颜色顺序，确保 K 辆车颜色循环
colors = get(gca, 'ColorOrder'); 
% 绘制 K 辆车, 3 种方案
for k = 1:K
    % 获取当前车辆的颜色
    current_color = colors(mod(k-1, size(colors, 1)) + 1, :);
    
    % 方案 1: RIS + Relay (实线)
    plot(x_U_lead, R_Total_platoon(:, k), ...
         'LineStyle', '-', 'LineWidth', 2, 'Color', current_color, ...
         'DisplayName', sprintf('车辆 %d - 总速率', k));
         
    % 方案 2: 仅 RIS (虚线)
    plot(x_U_lead, R_RIS_platoon(:, k), ...
         'LineStyle', '--', 'LineWidth', 1.5, 'Color', current_color, ...
         'DisplayName', sprintf('车辆 %d - 仅 RIS', k));
         
    % 方案 3: 仅 Relay (点线)
    plot(x_U_lead, R_Relay_platoon(:, k), ...
         'LineStyle', ':', 'LineWidth', 1.5, 'Color', current_color, ...
         'DisplayName', sprintf('车辆 %d - 仅 Relay', k));
end
% --- 添加图例 ---
% 创建虚拟线条来表示线型
h_style = zeros(3, 1);
h_style(1) = plot(NaN,NaN,'-k', 'LineWidth', 2); % 黑色实线
h_style(2) = plot(NaN,NaN,'--k', 'LineWidth', 1.5); % 黑色虚线
h_style(3) = plot(NaN,NaN,':k', 'LineWidth', 1.5); % 黑色点线
% 创建虚拟线条来表示车辆
h_vehicle = zeros(K, 1);
vehicle_labels = cell(K, 1);
for k = 1:K
    current_color = colors(mod(k-1, size(colors, 1)) + 1, :);
    h_vehicle(k) = plot(NaN, NaN, 'Color', current_color, 'LineWidth', 2);
    if k == 1
        vehicle_labels{k} = '车辆 1 (头车)';
    else
        vehicle_labels{k} = sprintf('车辆 %d', k);
    end
end
% 合并图例
legend([h_style; h_vehicle], ...
       {'RIS + Relay', '仅 RIS', '仅 Relay', vehicle_labels{:}}, ...
       'Location', 'best', 'NumColumns', 2);
title(sprintf('编队性能综合对比 (K=%d 辆, 间距=%.0fm) - 5G 物理模型 (Sub-6GHz)', K, d_platoon));
xlabel('头车位置 x_U [m]');
ylabel('频谱效率 [bits/s/Hz]');
grid on;
box on;
xlim([x_U_lead(1) x_U_lead(end)]);
ylim_max = max(R_Total_platoon, [], 'all') * 1.1;
ylim([0, ylim_max]);
hold off;
%% 辅助函数：计算 FSPL
function pl_dB = calculatePathLoss_dB(d_m, fc_Hz, c)
    % 检查距离是否为零或太近，避免 log(0)
    if d_m < 1.0
        d_m = 1.0; % 最小距离设为 1 米
    end
    
    % FSPL (dB) = 20*log10(d) + 20*log10(f) + 20*log10(4*pi/c)
    % 这是 MATLAB 'fspl' 函数的底层公式 (在 Antenna Toolbox 中)
    lambda = c / fc_Hz;
    pl_dB = 20 * log10(4 * pi * d_m / lambda);
end