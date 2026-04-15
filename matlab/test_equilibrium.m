% 测试平衡点计算逻辑
clear; clc;

% 简化参数
fc = 3.5e9;
c = 299792458;
bw = 100e6;
thermalNoise_dBm_per_Hz = -174;
rxNoiseFigure_dB = 7;
noisePower_dBm = thermalNoise_dBm_per_Hz + 10*log10(bw) + rxNoiseFigure_dB;

% 节点位置
S_pos = [100, 600];
I_pos = [0, 400];
R_pos = [100, 400];  % 测试：d_IR = 100m

% 功率参数
txPower_S_dBm = 30;
txGain_S_dBi = 20;
rxGain_D_dBi = 10;
txPower_R_dBm = 20;
txGain_R_dBi = 10;
rxGain_R_dBi = 10;

N = 10000;
eta = 1;

% 计算固定链路
d_SI = norm(S_pos - I_pos);
pl_SI_dB = calculatePathLoss_dB(d_SI, fc, c);
d_SR = norm(S_pos - R_pos);
pl_SR_dB = calculatePathLoss_dB(d_SR, fc, c);

% RIS增益
ris_gain_term = (N * (16 + (N-1)*pi^2)) / 16;
ris_gain_dB = 10*log10(ris_gain_term);

fprintf('=== 固定参数 ===\n');
fprintf('d_SI = %.1f m, PL = %.2f dB\n', d_SI, pl_SI_dB);
fprintf('d_SR = %.1f m, PL = %.2f dB\n', d_SR, pl_SR_dB);
fprintf('RIS增益 = %.2f dB\n', ris_gain_dB);
fprintf('底噪 = %.2f dBm\n\n', noisePower_dBm);

% 测试几个用户位置
x_test = [0, 50, 100, 150, 200, 250, 300];

fprintf('=== 速率计算 ===\n');
fprintf('x_U (m) | R_RIS (bps/Hz) | R_Relay (bps/Hz) | 差值\n');
fprintf('--------|----------------|------------------|--------\n');

for x_U = x_test
    D_pos = [x_U, 0];
    
    % RIS路径
    d_ID = norm(I_pos - D_pos);
    pl_ID_dB = calculatePathLoss_dB(d_ID, fc, c);
    Prx_RIS_dBm = txPower_S_dBm + txGain_S_dBi + rxGain_D_dBi - pl_SI_dB - pl_ID_dB + ris_gain_dB;
    snr_RIS_dB = Prx_RIS_dBm - noisePower_dBm;
    R_RIS = log2(1 + 10^(snr_RIS_dB/10));
    
    % 中继路径
    d_RD = norm(R_pos - D_pos);
    pl_RD_dB = calculatePathLoss_dB(d_RD, fc, c);
    Prx_SR_dBm = txPower_S_dBm + txGain_S_dBi + rxGain_R_dBi - pl_SR_dB;
    Prx_RD_dBm = txPower_R_dBm + txGain_R_dBi + rxGain_D_dBi - pl_RD_dB;
    snr_SR_dB = Prx_SR_dBm - noisePower_dBm;
    snr_RD_dB = Prx_RD_dBm - noisePower_dBm;
    snr_Relay_lin = min(10^(snr_SR_dB/10), 10^(snr_RD_dB/10));
    R_Relay = 0.5 * log2(1 + snr_Relay_lin);
    
    fprintf('%7.0f | %14.2f | %16.2f | %+7.2f\n', x_U, R_RIS, R_Relay, R_RIS - R_Relay);
end

function pl_dB = calculatePathLoss_dB(d_m, fc_Hz, c)
    if d_m < 1.0
        d_m = 1.0; 
    end
    lambda = c / fc_Hz;
    pl_dB = 20 * log10(4 * pi * d_m / lambda);
end
