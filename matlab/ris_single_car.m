%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB 脚本：复现专利 CN 114584587 B (修正版)
%
% 对比三种场景 (遵照专利模型，S->D 路径被阻挡):
% 1. RIS + Relay  (辅助中继)
% 2. RIS          (单RIS)
% 3. Relay        (单中继)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; 
clc;
close all;

%% 1. 定义系统参数 (基于专利)
% -------------------------------------------------------------------------
N = 1000;         % RIS 反射单元数量 (N=1000 用于图2)
eta = 1;          % RIS 反射幅度 (假设为1)

% 节点坐标 (基于 "数值结果与分析")
S_pos = [250, 833];  % 源 (S) 位置
I_pos = [0, 400];    % RIS (I) 位置
R_pos = [500, 400];  % 中继 (R) 位置

% 用户 (D) 沿x轴移动
x_U = 0:10:500;      % 用户位置 x_U [m]
D_pos_x = x_U;
D_pos_y = zeros(size(x_U));

%% 2. 假设的信道参数
% -------------------------------------------------------------------------
% 调整 rho 和 alpha 以匹配专利图2的形状和量级
% rho = P/sigma^2
rho_lin = 2.7e11;   % 假设：线性信噪比 (P/sigma^2)
alpha = 3.0;      % 假设：路径损耗指数 alpha

%% 3. 计算距离和信道增益
% -------------------------------------------------------------------------
% 计算固定距离
d_SI = norm(S_pos - I_pos); % S -> I
d_SR = norm(S_pos - R_pos); % S -> R

% 计算随用户移动而变化的距离
d_ID = zeros(size(x_U)); % I -> D
d_RD = zeros(size(x_U)); % R -> D

for i = 1:length(x_U)
    D_pos = [D_pos_x(i), D_pos_y(i)];
    d_ID(i) = norm(I_pos - D_pos); 
    d_RD(i) = norm(R_pos - D_pos); 
end

% 计算路径损耗 (平均信道增益 d^-alpha)
E_beta_SI = d_SI ^ (-alpha);
E_beta_SR = d_SR ^ (-alpha);
E_beta_ID = d_ID .^ (-alpha);
E_beta_RD = d_RD .^ (-alpha);

%% 4. 计算可达速率 (三种场景)
% -------------------------------------------------------------------------

% --- 场景 1: 中继 (Relay) 速率: R_Relay ---
% 专利公式 (3), 0.5 是因为半双工 
% R_DF = 0.5 * log2(1 + rho * min(E[|h_SR|^2], E[|h_RD|^2]))
min_Relay_Path_Gain = min(E_beta_SR, E_beta_RD);
R_Relay = 0.5 * log2(1 + rho_lin * min_Relay_Path_Gain);

% --- 场景 2: RIS 速率: R_RIS ---
% 专利公式 (6) 和 (10)
% E[gamma_RIS] = (1/16) * rho * d_SI^-a * d_ID^-a * eta^2 * N * (16 + (N-1)*pi^2)
ris_gain_term = (N * (16 + (N-1)*pi^2)) / 16;
E_gamma_RIS = rho_lin * (eta^2) * E_beta_SI .* E_beta_ID * ris_gain_term;
R_RIS = log2(1 + E_gamma_RIS); %

% --- 场景 3: 总速率 (RIS + Relay): R_Total ---
% 专利公式 (7) 和 (12): R_Sum = R_RIS + R_DF
R_Total = R_Relay + R_RIS;

%% 5. 绘图 (复现图2)
% -------------------------------------------------------------------------
figure;
plot(x_U, R_Total, '-o', 'LineWidth', 2, 'DisplayName', 'RIS + Relay (辅助中继)');
hold on;
plot(x_U, R_RIS, '-s', 'LineWidth', 2, 'DisplayName', 'RIS (单RIS)');
plot(x_U, R_Relay, '-d', 'LineWidth', 2, 'DisplayName', 'Relay (单中继)');

title('三种车联网通信方案对比 (基于专利 CN 114584587 B)');
xlabel('用户位置 x_U [m]');
ylabel('Rate [bits/s/Hz]');
legend('show', 'Location', 'west');
grid on;
xlim([0 500]);

% *** YLIM 错误修复 ***
% 获取当前y轴的最大值
y_max = max(R_Total); 
% 设置y轴范围，从0到最大值的1.1倍 (留一点空白)
ylim([0, y_max * 1.1]); 

hold off;