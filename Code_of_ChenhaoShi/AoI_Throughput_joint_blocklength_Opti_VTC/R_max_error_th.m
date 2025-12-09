clc;
clear;

% 用户数量和参数设置
I = 2; % 用户数量
d_array = [64, 52]; 
SNR_array = [3, 5];

% 固定 AoI 阈值和 M_th 值
aoi_th = 120;
M_th = 80;

% 设置 error_th 的范围
error_th_range = 0.001:0.005:0.1; % error_th 从 0.001 到 0.1
R_max_values_fixed = zeros(1, length(error_th_range)); % 固定 m 下的 R_max
optimal_aoi_values_fixed = zeros(1, length(error_th_range)); % 固定 m 下的最优 AoI
maximum_aoi_values_fixed = zeros(1, length(error_th_range)); % 固定 m 下的最大 AoI

R_max_values_variable = zeros(1, length(error_th_range)); % 变量 m 下的 R_max
optimal_aoi_values_variable = zeros(1, length(error_th_range)); % 变量 m 下的最优 AoI
maximum_aoi_values_variable = zeros(1, length(error_th_range)); % 变量 m 下的最大 AoI

% 计算每个用户的固定 m 值
fixed_m_values = zeros(1, I);
for i = 1:I
    C = log2(1 + SNR_array(i)); % 计算 C
    fixed_m_values(i) = d_array(i) / C; % 计算固定 m 值
end

% 计算固定 m 值下的 R_max, optimal_aoi, 和 maximum_aoi
for idx = 1:length(error_th_range)
    error_th = error_th_range(idx);
    
    % 固定 m 值的计算
    err_fixed = zeros(1, I);
    tau_fixed = zeros(1, I);
    aoi_fixed = zeros(1, I);

    for i = 1:I
        err_fixed(i) = error_prob_fbl(SNR_array(i), fixed_m_values(i), d_array(i) / fixed_m_values(i));
        tau_fixed(i) = d_array(i) / fixed_m_values(i) * (1 - err_fixed(i));
        aoi_fixed(i) = 0.5 * fixed_m_values(i) + fixed_m_values(i) / (1 - err_fixed(i));
    end

%     if all(err_fixed <= error_th) && all(aoi_fixed <= aoi_th) && (sum(fixed_m_values) <= M_th)
        R_max_values_fixed(idx) = sum(tau_fixed);
        optimal_aoi_values_fixed(idx) = min(aoi_fixed);
        maximum_aoi_values_fixed(idx) = max(aoi_fixed);

%     end
end

% 计算不固定 m 值下的 R_max, optimal_aoi, 和 maximum_aoi
for idx = 1:length(error_th_range)
    error_th = error_th_range(idx);
    
    % 初始化
    max_sum_tau = -inf; % 最大 sum_tau 初始化为负无穷
    optimal_aoi = inf;  % 最小 AOI 初始化为正无穷
    maximum_aoi = -inf; % 最大 AOI 初始化为负无穷

    % 遍历所有可能的 m 值组合
    for m1 = 1:500
        for m2 = 1:500
            m_current = [m1, m2]; % 当前组合的 m 值

            % 初始化
            err = zeros(1, I);
            tau = zeros(1, I);
            aoi = zeros(1, I);

            for i = 1:I
                err(i) = error_prob_fbl(SNR_array(i), m_current(i), d_array(i) / m_current(i));
                tau(i) = d_array(i) / m_current(i) * (1 - err(i));
                aoi(i) = 0.5 * m_current(i) + m_current(i) / (1 - err(i));
            end

            % 检查是否满足约束条件
            if all(err <= error_th) && all(aoi <= aoi_th) && (sum(m_current) <= M_th)
                current_sum_tau = sum(tau);
                
                % 记录 R_max
                if current_sum_tau > max_sum_tau
                    max_sum_tau = current_sum_tau;
                    R_max_values_variable(idx) = max_sum_tau;
                    optimal_aoi_values_variable(idx) = min(aoi);
                    maximum_aoi_values_variable(idx) = max(aoi);
                end
            end
        end
    end
end

% 绘图
figure;
hold on;
yyaxis left;
plot(error_th_range, R_max_values_fixed, '-o', 'Color', 'b', 'DisplayName', 'Fixed m - Sum Throughput R_{max}');
plot(error_th_range, R_max_values_variable, '-x', 'Color', 'g', 'DisplayName', 'Variable m - Sum Throughput R_{max}');
xlabel('Error Threshold \varspsilon');
ylabel('Maximum Throughput R_{max}');
grid on;

% 绘制右轴曲线：optimal aoi 和 maximum aoi
yyaxis right;
plot(error_th_range, optimal_aoi_values_fixed, '--*', 'Color', 'r', 'DisplayName', 'Fixed m - Optimal AoI');
plot(error_th_range, maximum_aoi_values_fixed, ':', 'Color', 'm', 'DisplayName', 'Fixed m - Maximum AoI');
plot(error_th_range, optimal_aoi_values_variable, '--s', 'Color', 'c', 'DisplayName', 'Variable m - Optimal AoI');
plot(error_th_range, maximum_aoi_values_variable, ':', 'Color', 'k', 'DisplayName', 'Variable m - Maximum AoI');
ylabel('AoI');
legend;
hold off;

% 输出固定 m 值
fprintf('Fixed m values: m1 = %.2f, m2 = %.2f\n', fixed_m_values(1), fixed_m_values(2));
