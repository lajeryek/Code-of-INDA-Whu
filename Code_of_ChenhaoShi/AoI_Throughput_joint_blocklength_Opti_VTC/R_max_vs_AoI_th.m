clc;
clear;

% 参数初始化
T_S = 3e-5;        % 符号长度
I = 2;             % 传感器数量
SNR_vector = [2, 5];       % 每个用户的信噪比
d_vector = [100, 110];     % 每个用户的数据量
T_max = 40e-3;     % 最大时延
error_max = 1e-1;  % 最大错误率

% 遍历的AOI阈值
AOI_th_array = 0.1:0.5:5;  % 不同的AOI阈值
R_max_values = zeros(1, length(AOI_th_array));  % 存储每个AOI_th对应的最大吞吐量

% 遍历每个AOI阈值
for idx = 1:length(AOI_th_array)
    AOI_th = AOI_th_array(idx);
    
    % 创建 m1 和 m2 的网格
    [m1, m2] = meshgrid(1:200, 1:200);
    
    % 计算每个用户的错误率
    error_array = [erro(m1, SNR_vector(1), d_vector(1)); erro(m2', SNR_vector(2), d_vector(2))];

    % 计算每个用户的AOI
    AOI_array = 0.5 * m1 * T_S + (m1 * T_S) ./ (1 - error_array(1, :));  % 用户 1 的 AOI
    AOI_array_2 = 0.5 * m2 * T_S + (m2 * T_S) ./ (1 - error_array(2, :));  % 用户 2 的 AOI

    % 检查是否满足约束条件
    valid_indices = (AOI_array <= AOI_th) & ...
                    all(error_array <= error_max, 1) & ...
                    ((m1 + m2) * T_S <= T_max);
    
    % 计算当前的 A_mi 和 B_mi
    A_mi = (1 - error_array(1, :)) .* d_vector(1) + (1 - error_array(2, :)) .* d_vector(2);
    B_mi = (m1 + m2) * T_S;

    % 计算吞吐量
    R_current = A_mi ./ B_mi;

    % 检查 valid_indices 是否有有效值
    if any(valid_indices)
        % 更新最大吞吐量
        R_max = max(R_current(valid_indices), [], 'omitnan'); % 忽略 NaN 值
    else
        R_max = 0; % 或者设定为一个合适的值
    end
    
    % 记录当前 AOI_th 对应的最大吞吐量
    R_max_values(idx) = R_max;
end


% 绘制 R_max 随 AOI_th 变化的曲线
figure;
plot(AOI_th_array, R_max_values, '-o', 'LineWidth', 2);
xlabel('\Delta_{th}');
ylabel('Maximum Throughput R_{max}');
grid on;
