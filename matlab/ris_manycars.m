%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB 脚本：多车辆编队行驶模型 (综合对比图)
%
% 扩展自 专利 CN 114584587 B 模型
% 场景: 一个 K 辆车的编队沿 x 轴行驶
% 对比: 在一张图中显示 K 辆车在 3 种方案下的性能
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; 
clc;
close all;

%% 1. 定义系统参数 (基于专利)
% -------------------------------------------------------------------------
N = 1000;         % RIS 反射单元数量 (N=1000)
eta = 1;          % RIS 反射幅度 (假设为1)

% 节点坐标 (基于 "数值结果与分析")
S_pos = [250, 833];  % 源 (S) 位置
I_pos = [0, 400];    % RIS (I) 位置
R_pos = [500, 400];  % 中继 (R) 位置

%% 2. 定义编队 (Platoon) 参数
% -------------------------------------------------------------------------
K = 4;            % 编队中的车辆数 (例如 4 辆)
d_platoon = 20;   % 车辆间距 (例如 20 米)

% 头车 (D_1) 沿x轴移动
x_U_lead = 0:10:500; % 头车位置 x_U [m]
num_positions = length(x_U_lead);

%% 3. 假设的信道参数
% -------------------------------------------------------------------------
rho_lin = 2.7e11;   % 假设：线性信噪比 (P/sigma^2)
alpha = 3.0;      % 假设：路径损耗指数 alpha

%% 4. 计算距离和信道增益
% -------------------------------------------------------------------------
% 计算固定距离和增益
d_SI = norm(S_pos - I_pos); % S -> I
d_SR = norm(S_pos - R_pos); % S -> R
E_beta_SI = d_SI ^ (-alpha);
E_beta_SR = d_SR ^ (-alpha);

% RIS 增益项 (常量)
ris_gain_term = (N * (16 + (N-1)*pi^2)) / 16;

% 初始化存储速率的矩阵 (行: 位置, 列: 车辆)
R_Relay_platoon = zeros(num_positions, K);
R_RIS_platoon   = zeros(num_positions, K);
R_Total_platoon = zeros(num_positions, K);

%% 5. 循环计算每辆车的速率
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
        
        % 5.3. 计算第k辆车的路径增益
        E_beta_ID_k = d_ID_k ^ (-alpha);
        E_beta_RD_k = d_RD_k ^ (-alpha);
        
        % 5.4. 计算第k辆车的速率
        
        % 场景 1: "仅中继" 速率
        min_Relay_Path_Gain_k = min(E_beta_SR, E_beta_RD_k);
        R_Relay_platoon(i, k) = 0.5 * log2(1 + rho_lin * min_Relay_Path_Gain_k);

        % 场景 2: "仅RIS" 速率
        E_gamma_RIS_k = rho_lin * (eta^2) * E_beta_SI * E_beta_ID_k * ris_gain_term;
        R_RIS_platoon(i, k) = log2(1 + E_gamma_RIS_k);

        % 场景 3: "RIS + 中继" 速率
        R_Total_platoon(i, k) = R_Relay_platoon(i, k) + R_RIS_platoon(i, k);
    end
end

%% 6. 绘图 (综合对比图)
% -------------------------------------------------------------------------
figure;
hold on;

% 设置颜色顺序，确保 K 辆车颜色循环
colors = get(gca, 'ColorOrder'); 

% 绘制 K 辆车, 3 种方案
for k = 1:K
    % 获取当前车辆的颜色
    current_color = colors(mod(k-1, size(colors, 1)) + 1, :);
    
    % 方案 1: RIS + Relay (实线)
    plot(x_U_lead, R_Total_platoon(:, k), ...
         'LineStyle', '-', 'LineWidth', 2, 'Color', current_color);
         
    % 方案 2: 仅 RIS (虚线)
    plot(x_U_lead, R_RIS_platoon(:, k), ...
         'LineStyle', '--', 'LineWidth', 1.5, 'Color', current_color);
         
    % 方案 3: 仅 Relay (点线)
    plot(x_U_lead, R_Relay_platoon(:, k), ...
         'LineStyle', ':', 'LineWidth', 1.5, 'Color', current_color);
end

% --- 添加图例 ---
% 创建虚拟线条来表示线型
h = zeros(3, 1);
h(1) = plot(NaN,NaN,'-k', 'LineWidth', 2); % 黑色实线
h(2) = plot(NaN,NaN,'--k', 'LineWidth', 1.5); % 黑色虚线
h(3) = plot(NaN,NaN,':k', 'LineWidth', 1.5); % 黑色点线
legend(h, 'RIS + Relay', '仅 RIS', '仅 Relay', 'Location', 'best');

% --- 添加车辆图例 (在图的右侧) ---
% 这比较复杂，一个更简单的方法是为车辆添加文本
% 在图的末尾为每种颜色添加标签
y_pos_legend = linspace(max(R_Total_platoon, [], 'all')*0.9, ...
                        max(R_Total_platoon, [], 'all')*0.7, K);
for k = 1:K
    current_color = colors(mod(k-1, size(colors, 1)) + 1, :);
    if k == 1
        label = '车辆 1 (头车)';
    else
        label = sprintf('车辆 %d', k);
    end
    text(x_U_lead(end), y_pos_legend(k), label, ...
         'Color', current_color, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
end


title(sprintf('编队性能综合对比 (K=%d 辆, 间距=%.0fm)', K, d_platoon));
xlabel('头车位置 x_U [m]');
ylabel('各车速率 [bits/s/Hz]');
grid on;
xlim([x_U_lead(1) x_U_lead(end)]);
ylim_max = max(R_Total_platoon, [], 'all') * 1.1;
ylim([0, ylim_max]);
hold off;