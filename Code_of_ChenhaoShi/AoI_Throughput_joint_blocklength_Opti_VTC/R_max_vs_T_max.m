clc
clear
m_1 = 1:1000; % 用户 1 的 m 值
m_2 = (1:1000)'; % 用户 2 的 m 值
M = m_1 + m_2; % 计算 M

% 定义不同的 AoI_th, error_th, SNR 参数组合 (d 固定为 64)
combinations = [
    120, 0.1, 3, 5    % AoI_th=120, error_th=0.1, SNR1=2, SNR2=5
    120, 0.5, 2, 5; % AoI_th=150, error_th=0.2, SNR1=2.5, SNR2=5
    120, 0.1,2, 5;   % AoI_th=120, error_th=0.1, SNR1=2, SNR2=5
    150, 0.1, 2, 2;   % AoI_th=150, error_th=0.2, SNR1=2, SNR2=5
    150, 0.5, 2, 5; % AoI_th=150, error_th=0.2, SNR1=2.5, SNR2=5
];

d_1 = 64; % 固定 d_1
d_2 = 64; % 固定 d_2

M_th_range = 4:4:100; % 横坐标 M_th 范围
num_thresholds = length(M_th_range); % M_th 阈值数量
legends = cell(size(combinations, 1), 1); % 存储图例

% 设置图形属性
figure
hold on
set(gca, 'FontSize', 16); % 坐标轴字体大小

for idx = 1:size(combinations, 1)
    % 提取组合参数
    AoI_th = combinations(idx, 1);
    error_th = combinations(idx, 2);
    SNR_1 = combinations(idx, 3);
    SNR_2 = combinations(idx, 4);
    
    % 计算错误概率
    err_1 = error_prob_fbl(SNR_1, m_1, d_1 ./ m_1);
    err_2 = error_prob_fbl(SNR_2, m_2, d_2 ./ m_2);
    
    [M_rows, M_cols] = size(M);
    
    % 将错误率扩展为矩阵
    err_1_matrix = repmat(err_1, M_rows, 1);
    err_2_matrix = repmat(err_2, 1, M_cols);
    
    % 计算 tau 值
    tau_1 = d_1 ./ m_1 .* (1 - err_1_matrix);
    tau_2 = d_2 ./ m_2 .* (1 - err_2_matrix);
    
    % 计算 sum_tau
    sum_tau = tau_1 + tau_2;
    
    % 计算 AOI
    aoi_1 = 0.5 .* M + M ./ (1 - err_1_matrix);
    
    % 计算最大 sum_tau 和对应的 m_1, m_2
    max_sum_tau_values = zeros(1, num_thresholds); % 存储每个 M_th 对应的最大 sum_tau
    m_1_opt_values = zeros(1, num_thresholds); % 存储每个 M_th 下的 m_1 最优值
    m_2_opt_values = zeros(1, num_thresholds); % 存储每个 M_th 下的 m_2 最优值

    % 遍历 M_th 值
    for i = 1:num_thresholds
        M_th = M_th_range(i);
        
        % 应用约束：error <= error_th 和 M < M_th，且 AoI <= AoI_th
        mask = (err_1_matrix <= error_th) & (M < M_th) & (aoi_1 <= AoI_th);
        
        % 计算符合条件的 sum_tau
        sum_tau_masked = sum_tau.*mask; 
        sum_tau_masked(~mask) = NaN; % 设定不符合条件的值为 NaN
        
        % 找到最大值及其索引，忽略 NaN 值
        [max_value, linear_index] = max(sum_tau_masked(:), [], 'omitnan');
        
        % 存储结果
        max_sum_tau_values(i) = max_value;
        [m_1_opt_values(i), m_2_opt_values(i)] = ind2sub(size(sum_tau_masked), linear_index);
    end
    
    % 绘制曲线
    plot(M_th_range, max_sum_tau_values*1e4, 'LineWidth', 3);
    
    % 添加图例
    legends{idx} = sprintf('$\\Delta_{th}=%d$, $\\varepsilon_{max}=%.1f$, $\\gamma_1=%.1f$, $\\gamma_2=%.1f$', ...
        AoI_th, error_th, SNR_1, SNR_2);

end

xlabel('Maximum Transmission Delay T_{max} (ms)', 'FontSize', 16); % 设置横坐标字体大小
ylabel('Maximum Throughput R_{max} (bits/s)', 'FontSize', 16); % 设置纵坐标字体大小
grid on;
legend(legends, 'Location', 'best', 'Interpreter', 'latex');
hold off;
