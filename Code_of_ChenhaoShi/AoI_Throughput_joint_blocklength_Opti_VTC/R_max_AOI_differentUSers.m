clc
clear
m_1 = 1:300; % 用户 1 的 m 值
m_2 = (1:300)'; % 用户 2 的 m 值
M = m_1 + m_2; % 计算 M

% 定义不同的 d1 和 SNR_1 的组合
combinations = [
    150, 3, 150, 5;  % d1=30bits, SNR_1=2; d_2=60bits, SNR=5
    30, 5, 60, 5;  % d1=30bits, SNR_1=5; d_2=60bits, SNR=5
    100, 2, 60, 5;  % d1=100bits, SNR_1=2; d_2=60bits, SNR=5
    64,  2, 64,4  
];

aoi_th_range = 100:5:160; % AOI 阈值范围
num_thresholds = length(aoi_th_range); % 阈值数量

% colors = ['r', 'g', 'b','g']; % 不同曲线的颜色
legends = cell(size(combinations, 1), 1); % 存储图例

figure
hold on

for idx = 1:size(combinations, 1)
    % 提取组合参数
    d_1 = combinations(idx, 1);
    SNR_1 = combinations(idx, 2);
    d_2 = combinations(idx, 3);
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
    max_sum_tau_values = zeros(1, num_thresholds); % 存储每个阈值对应的最大 sum_tau
    m_1_opt_values = zeros(1, num_thresholds); % 存储每个阈值下的 m_1 最优值
    m_2_opt_values = zeros(1, num_thresholds); % 存储每个阈值下的 m_2 最优值

    % 遍历 aoi_th 值
    for i = 1:num_thresholds
        aoi_th = aoi_th_range(i);
        
        % 应用约束：error <= 1e-3 和 M < 1000
        mask = (err_1_matrix <= 1e-2) & (M < 1000) & (aoi_1 <= aoi_th);
        
        % 计算符合条件的 sum_tau
        sum_tau_masked = sum_tau; 
        sum_tau_masked(~mask) = NaN; % 设定不符合条件的值为 NaN
        
        % 找到最大值及其索引，忽略 NaN 值
        [max_value, linear_index] = max(sum_tau_masked(:), [], 'omitnan');
        
        % 存储结果
        max_sum_tau_values(i) = max_value;
        [m_1_opt_values(i), m_2_opt_values(i)] = ind2sub(size(sum_tau_masked), linear_index);
    end
    
    % 绘制曲线
    plot(aoi_th_range*1e-3, max_sum_tau_values*1e4);
    
    % 添加图例
    legends{idx} = sprintf('d1=%dbits, SNR_1=%d; d2=%dbits, SNR_2=%d', d_1, SNR_1, d_2, SNR_2);
end

xlabel('\Delta_{th}');
ylabel('Maximum Throughput R_{max}');
grid on;
legend(legends, 'Location', 'best');
hold off
