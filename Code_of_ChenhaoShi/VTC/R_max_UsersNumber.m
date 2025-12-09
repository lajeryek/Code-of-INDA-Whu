clc
clear

% 用户数量和蒙特卡洛模拟次数
I = 3; % 用户数量
num_iterations = 1000; % 蒙特卡洛模拟次数

% 固定 AoI 阈值和 error 阈值
aoi_th = 120;
error_th = 0.1;  % 固定 error 阈值
M_th = 2000; % 固定的 M_th 值

% 结果存储
I_range = 1:I; % 用户数量的范围
R_max_values = zeros(1, length(I_range)); % 存储每个 I 对应的最大 sum_tau
optimal_aoi_values = zeros(1, length(I_range)); % 存储每个 I 对应的最优 AoI
maximum_aoi_values = zeros(1, length(I_range));

% 启动并行计算
figure
hold on
yyaxis left

parfor I_idx = I_range
    R_max_iter = zeros(1, num_iterations); 
    optimal_aoi_iter = zeros(1, num_iterations); 
    maximum_aoi_iter = zeros(1, num_iterations);
    
    % 蒙特卡洛模拟
    for k = 1:num_iterations
        [R_max, optimal_aoi, maximum_aoi] = solveUsingCVX(I_idx, M_th, error_th, aoi_th);
        
        % 存储每次模拟结果
        R_max_iter(k) = R_max;
        optimal_aoi_iter(k) = optimal_aoi;
        maximum_aoi_iter(k) = maximum_aoi;
    end
    
    % 计算平均值
    R_max_values(I_idx) = mean(R_max_iter, 'omitnan');
    optimal_aoi_values(I_idx) = mean(optimal_aoi_iter, 'omitnan');
    maximum_aoi_values(I_idx) = mean(maximum_aoi_iter, 'omitnan');
end

% 绘制结果
plot(I_range, R_max_values , '-o', 'Color', 'b', 'DisplayName', ['Sum Throughput and M_{th} = ' num2str(M_th)]);
xlabel('Number of Users (I)');
ylabel('Maximum Throughput R_{max}');
grid on;

yyaxis right
plot(I_range, optimal_aoi_values, '--*', 'Color', 'r', 'DisplayName', ['min\{\Delta^*\} and M_{th}=' num2str(M_th)]);
plot(I_range, maximum_aoi_values, ':', 'Color', 'k', 'DisplayName', ['max\{\Delta\} and M_{th}=' num2str(M_th)]);
ylabel('Optimal AoI');
legend
hold off

% CVX方法计算最大吞吐量和AoI
function [R_max, optimal_aoi, maximum_aoi] = solveUsingCVX(I, M_th, error_th, aoi_th)
    % 随机生成 d 和 SNR 值
    d = randi([24, 256], 1, I); 
    SNR = randi([2, 8], 1, I);
    
    % 初始化变量
    R_max = NaN;
    optimal_aoi = NaN;
    maximum_aoi = NaN;

    % 使用CVX求解优化问题
    cvx_begin quiet
        variables m(I) % 每个用户的码字长度
        variables err(I) tau(I) aoi(I) % 错误概率、吞吐量、AoI
        
        % 错误概率和吞吐量计算公式
        for i = 1:I
            err(i) = error_prob_fbl(SNR(i), m(i), d(i) / m(i));
            tau(i) = d(i) * (1 - err(i)) / m(i);
            aoi(i) = 0.5 * m(i) + m(i) / (1 - err(i));
        end

        % 优化目标：最大化吞吐量
        maximize sum(tau)
        
        % 约束条件
        subject to
            sum(m) <= M_th; % 总资源约束
            err <= error_th; % 错误概率约束
            aoi <= aoi_th; % AoI 约束
            m >= 1; % m 必须是正值
    cvx_end

    % 返回结果
    if strcmp(cvx_status, 'Solved')
        R_max = sum(tau);
        optimal_aoi = min(aoi);
        maximum_aoi = max(aoi);
    end
end