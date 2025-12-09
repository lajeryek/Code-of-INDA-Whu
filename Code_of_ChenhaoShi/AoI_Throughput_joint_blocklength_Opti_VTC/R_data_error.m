clc;
clear;
m_1 = 1:300; % 用户 1 的 m 值
m_2 = (1:300)'; % 用户 2 的 m 值
M = m_1 + m_2; % 计算 M
aoi_th=250;
% 定义不同的 error_th 和 SNR 的组合
combinations = [
    1e-2, 3, 3;
    1e-2, 3, 5;
    1e-3, 3, 5;
    1e-3, 3, 3
];

d_range = 24:5:256; % d1 和 d2 的范围
num_d_values = length(d_range); % d 值的数量

legends = cell(size(combinations, 1), 1); % 存储图例

figure
hold on

for idx = 1:size(combinations, 1)
    % 提取组合参数
    error_th = combinations(idx, 1);
    SNR_1 = combinations(idx, 2);
    SNR_2 = combinations(idx, 3);
    
    max_sum_tau_values = zeros(1, num_d_values); % 存储每个 d 值对应的最大 sum_tau

    for d_idx = 1:num_d_values
        d_1 = d_range(d_idx);
        d_2 = d_range(d_idx);
        
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
        
        % 应用约束：error <= error_th 和 M < 1000
        mask = (err_1_matrix <= error_th) & (M < 1000) & (aoi_1 <= aoi_th);
        
        % 计算符合条件的 sum_tau
        sum_tau_masked = sum_tau; 
        sum_tau_masked(~mask) = NaN; % 设定不符合条件的值为 NaN
        
        % 找到最大值，忽略 NaN 值
        max_value = max(sum_tau_masked(:), [], 'omitnan');
        
        % 存储结果
        max_sum_tau_values(d_idx) = max_value;
    end
    
    % 绘制曲线，设置线宽为 3，标记大小为 8
    plot(d_range, max_sum_tau_values*1e4, 'LineWidth', 3, 'MarkerSize', 8);
    
    % 添加图例
    legends{idx} = sprintf('error_th=%.e, SNR_1=%d, SNR_2=%d', error_th, SNR_1, SNR_2);
end

% 设置字号为 18
xlabel('数据包大小 (比特)', 'FontSize', 18);
ylabel('有效吞吐量 (bits/s)', 'FontSize', 18);
grid on;
legend(legends, 'Location', 'best', 'FontSize', 18);
hold off    