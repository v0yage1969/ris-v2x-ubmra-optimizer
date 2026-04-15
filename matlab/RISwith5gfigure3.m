%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB 脚本：复现专利 图 3 (速率 vs. 位置, 对比 N)
%
% 对应专利：图 3
% 仿真内容：总速率 (R_Total) vs. 用户位置 (x_U)
% 变化参数：RIS 单元数 (N)
%
% 注意：
% 1. 为了复现专利图 3 中 x=0 处的 "拐点", 增加了 RIS 盲区 (x < 0)。
% 2. 仅仿真 K=1 (单用户 D)，即头车。
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; 
clc;
close all;
%% 1. 定义系统参数 (基于专利)
% -------------------------------------------------------------------------
N_values = [2500, 5000, 7500, 10000]; % RIS 单元数量 (数组)
eta = 1;          % RIS 反射幅度 (假设为1)
% 节点坐标 (基于 "数值结果与分析")
S_pos = [250, 833];  % 源 (S) 位置
I_pos = [0, 400];    % RIS (I) 位置
R_pos = [500, 400];  % 中继 (R) 位置
%% 1.5. 定义 5G 和通信物理参数 (新)
% -------------------------------------------------------------------------
fc = 3.5e9;          % 载波频率 (Hz) - 3.5 GHz (Sub-6 GHz)
bw = 100e6;         % 系统带宽 (Hz) - 100 MHz
c = 299792458;       % 光速 (m/s)
txPower_dBm = 30;   % 发射功率 (dBm) - (例如 1W)
txGain_dBi = 20;    % 发射天线增益 (dBi)
rxGain_dBi = 10;    % 接收天线增益 (dBi)
rxNoiseFigure_dB = 7; % 接收机噪声系数 (dB)
% 计算热噪声功率 (dBm)
thermalNoise_dBm_per_Hz = -174; 
noisePower_dBm = thermalNoise_dBm_per_Hz + 10*log10(bw) + rxNoiseFigure_dB;
fprintf('5G 物理参数:\n');
fprintf('  载波频率: %.1f GHz\n', fc/1e9);
fprintf('  发射功率: %.0f dBm\n', txPower_dBm);
fprintf('  接收机底噪: %.2f dBm\n', noisePower_dBm);
%% 2. 定义编队 (Platoon) 参数
% -------------------------------------------------------------------------
K = 1;            % 对应专利图，只考虑单用户 (D)
% 头车 (D_1) 沿x轴移动
% *** 增加 x < 0 的部分来复现专利图 3 的盲区 ***
x_U_range = -100:10:500; % 头车位置 x_U [m]
num_positions = length(x_U_range);
%% 4. 计算固定距离和增益
% -------------------------------------------------------------------------
% 计算固定距离
d_SI = norm(S_pos - I_pos); % S -> I
d_SR = norm(S_pos - R_pos); % S -> R
% 计算固定路径损耗 (dB)
pl_SI_dB = calculatePathLoss_dB(d_SI, fc, c);
pl_SR_dB = calculatePathLoss_dB(d_SR, fc, c);
% 初始化存储
R_Total_results = zeros(num_positions, length(N_values));
%% 5. 循环计算速率
% -------------------------------------------------------------------------
for n_idx = 1:length(N_values)
    
    N = N_values(n_idx);
    
    % RIS 增益项 (来自专利)
    ris_gain_term = (N * (16 + (N-1)*pi^2)) / 16;
    ris_gain_dB = 10*log10(ris_gain_term); % 将线性增益转换为dB
    fprintf('正在计算 N = %d, RIS 增益 = %.2f dB\n', N, ris_gain_dB);

    for i = 1:num_positions  % 遍历头车的每个位置
        
        % 5.1. 计算第k辆车的当前x坐标
        D_k_x_pos = x_U_range(i);
        D_k_y_pos = 0; % 假设都在x轴上行驶
        D_k_pos = [D_k_x_pos, D_k_y_pos];
        
        % 5.2. 计算第k辆车的距离
        d_ID_k = norm(I_pos - D_k_pos); % I -> D_k
        d_RD_k = norm(R_pos - D_k_pos); % R -> D_k
        
        % 5.3. 计算第k辆车的路径损耗 (dB)
        pl_IDk_dB = calculatePathLoss_dB(d_ID_k, fc, c);
        pl_RDk_dB = calculatePathLoss_dB(d_RD_k, fc, c);
        
        % 5.4. 计算速率 (bits/s/Hz)
        
        % --- 场景 1: "仅中继" 速率 (半双工) ---
        Prx_SR_dBm = txPower_dBm + txGain_dBi + rxGain_dBi - pl_SR_dB;
        Prx_RDk_dBm = txPower_dBm + txGain_dBi + rxGain_dBi - pl_RDk_dB;
        snr_SR_dB = Prx_SR_dBm - noisePower_dBm;
        snr_RDk_dB = Prx_RDk_dBm - noisePower_dBm;
        snr_Relay_lin = min( 10^(snr_SR_dB/10), 10^(snr_RDk_dB/10) );
        R_Relay = 0.5 * log2(1 + snr_Relay_lin);
        
        % --- 场景 2: "仅RIS" 速率 ---
        Prx_RISk_dBm = txPower_dBm + txGain_dBi + rxGain_dBi - pl_SI_dB - pl_IDk_dB + ris_gain_dB;
        snr_RIS_dB = Prx_RISk_dBm - noisePower_dBm;
        snr_RIS_lin = 10^(snr_RIS_dB/10);
        
        % *** 专利中的盲区 (Fig 3): 假设 x < 0 时 RIS 速率为 0 ***
        if D_k_x_pos < 0
            R_RIS = 0; 
        else
            R_RIS = log2(1 + snr_RIS_lin);
        end
        
        % --- 场景 3: "RIS + 中继" 速率 ---
        R_Total_results(i, n_idx) = R_Relay + R_RIS;
    end
end
%% 6. 绘图 (复现 图 3)
% -------------------------------------------------------------------------
figure('Name', '图 3 (复现): 速率 vs. 位置 (对比 N)');
hold on;
markers = {'-o', '-s', '-d', '-^'};
for n_idx = 1:length(N_values)
    plot(x_U_range, R_Total_results(:, n_idx), ...
         markers{mod(n_idx-1, length(markers)) + 1}, ...
         'LineWidth', 1.5, 'MarkerSize', 5, ...
         'DisplayName', sprintf('N = %d', N_values(n_idx)));
end
title('总速率 vs. 用户位置 (对比 N 值)');
xlabel('用户位置 x_U [m]');
ylabel('频谱效率 [bits/s/Hz]');
grid on;
box on;
legend('Location', 'best');
xlim([x_U_range(1) x_U_range(end)]);
ylim_max = max(R_Total_results, [], 'all') * 1.1;
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