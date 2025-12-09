clc; clear;

% 参数设置
I_values = 1:1:10; % 用户数量范围
aoi_th = 150; % AOI 阈值
error_th = 0.1; % 误码率阈值
M_th = 1500; % 时延约束
MC_times = 10; % 每个I值蒙特卡洛模拟次数

% 初始化存储平均最大吞吐量、最小 AOI 和最大 AOI 的数组
R_max_averages = zeros(length(I_values), 1);
AOI_min_averages = zeros(length(I_values), 1);
AOI_max_averages = zeros(length(I_values), 1);

% 开启并行池（可根据CPU核数调整）
if isempty(gcp('nocreate'))
    parpool('local');
end

% 遍历不同用户数量 I
for idx = 1:length(I_values)
    I = I_values(idx);

    % 并行蒙特卡洛模拟
    R_all = zeros(MC_times, 1);
    AOI_min_all = zeros(MC_times, 1);
    AOI_max_all = zeros(MC_times, 1);

    parfor mc = 1:MC_times
        % 随机生成数据包大小和SNR
        d = randi([50, 100], I, 1);       % 数据包大小
        SNR = 1 + 4 * rand(I, 1);         % SNR 范围 [1, 5]

        x0 = ones(I, 1) * 200; % 初始时隙分配
        x = x0;
        R = 0;
        eta = 10; % 初始学习率

        for k = 1:5000 % 最大迭代限制
            R_old = R;
            grad = -grad_throughput(x, d, SNR, I);
            x = x - eta * grad;
            eta = eta / (1 + 0.001 * k); % 自适应学习率衰减

            x = max(1, floor(x)); % 向下取整并保证非零

            R = throughput(x, d, SNR, I);
            err_probs = error_prob_fbl(SNR, x, d ./ x);
            X = sum(x);
            AoI = 0.5 * X + X ./ (1 - err_probs);

            if abs(R - R_old) <= 1e-2 && max(AoI) <= aoi_th && max(err_probs) <= error_th && sum(x) <= M_th
                break;
            end
        end

        R_all(mc) = R;
        AOI_min_all(mc) = min(AoI);
        AOI_max_all(mc) = max(AoI);
    end

    % 统计平均值
    R_max_averages(idx) = mean(R_all);
    AOI_min_averages(idx) = mean(AOI_min_all);
    AOI_max_averages(idx) = mean(AOI_max_all);
end

% 绘图
figure;
yyaxis left
plot(I_values, R_max_averages, '-o', 'LineWidth', 1.5);
xlabel('Number of Users (I)');
ylabel('Average Maximum Throughput');
title('Throughput and AoI vs. Number of Users');
grid on;

yyaxis right
plot(I_values, AOI_min_averages, '--s', 'LineWidth', 1.5);
hold on;
plot(I_values, AOI_max_averages, '--d', 'LineWidth', 1.5);
ylabel('Average AoI (min/max)');
legend('Throughput', 'Min AOI', 'Max AOI');

%---------------- 辅助函数 -----------------

function grad = grad_error(m_n, SNR_normal, d)
    V_r = 1 - (1 + SNR_normal).^(-2);  % 信道色散
    C_r = log2(1 + SNR_normal);        % 香农容量

    % qfunc的导数部分
    qfunc_derivative = -1 / sqrt(2 * pi) * exp(-(((m_n / V_r)^(1/2)) * (C_r - d / m_n))^2 / 2);

    % 对内部项求导
    term1 = 1 / (2 * sqrt(V_r * m_n));  % 对 sqrt(m_n / V_r) 的求导
    term2 = d / (m_n^2);  % 对 (C_r - d / m_n) 的求导

    % 导数
    grad = qfunc_derivative * (term1 * (C_r - d / m_n) + sqrt(m_n / V_r) * term2) * log(2);
    
end
function grad = grad_aoi(m, SNR, d)
    I = length(m);
    grad = zeros(I, 1); % Initialize gradient as a column vector
    for i = 1:I
        err = error_prob_fbl(SNR(i), m(i), d(i) ./ m(i));
        grad_err = grad_error(m(i), SNR(i), d(i));  % grad_error should return scalar
        grad(i) = 0.5 + 1 / (1 - err) +(m(i) / (1 - err)^2) * grad_err;  
    end
    grad=grad(:);
end



function tau = throughput(m, d, SNR, I)
    numerator = 0;
    denominator = sum(m(:));
    for i = 1:I
        err = error_prob_fbl(SNR(i), m(i), d(i) / m(i));
        numerator = numerator + (1 - err) * d(i);
    end
    tau = numerator / denominator;
end

function grad = grad_throughput(m, d, SNR, I)
    denominator = sum(m,'all');  % 总和 m
    numerator = 0;         % 初始化分子
    grad = zeros(I, 1);    % 梯度初始化为列向量
    %计算分子部分求和（1-err）d
    for i = 1:I
        % 计算误码概率及其梯度
        err = error_prob_fbl(SNR(i), m(i), d(i) / m(i));
        numerator = numerator + d(i) * (1 - err); % 更新分子
    end
    for i=1:I
        % 计算对每一个用户m_i的梯度
        grad_err = grad_error(m(i), SNR(i), d(i));
        grad(i) = -(d(i) * grad_err * denominator + numerator) / denominator^2;
    end
    grad=grad(:);
end



% ======================= Supporting Functions =======================
function [m_opt] = solve_m_given_error(epsilon, SNR_normal, d)
    V_r = 1 - (1 + SNR_normal).^(-2);
    C_r = log2(1 + SNR_normal);
    Q_inv_epsilon = sqrt(2) * erfcinv(2 * epsilon);

    objective_function = @(m_n) sqrt(m_n/V_r) * (C_r - d./m_n) * log(2) - Q_inv_epsilon;
    options = optimoptions('fsolve', 'Display', 'off', 'MaxIterations', 100);

    % 初始化搜索范围，考虑更广的初始值
    m_initial_guess = max(d^2 / C_r, 1); % 避免初始值太小
    m_opt = fsolve(objective_function, m_initial_guess, options);

    % 检查并修正约束边界
    if m_opt < 1
        m_opt = 1;
    elseif m_opt > 500
        m_opt = 500;
    end
end


function [m_opt] = solve_m_given_AOI(AOI_target, SNR_normal, d)
    V_r = 1 - (1 + SNR_normal).^(-2);
    C_r = log2(1 + SNR_normal);
    erro_Q = @(m_n) qfunc(sqrt(m_n/V_r) * (C_r - d./m_n) * log(2));

    objective_function = @(m_n) 0.5 * m_n + m_n / (1 - erro_Q(m_n)) - AOI_target;
    options = optimoptions('fsolve', 'Display', 'off', 'MaxIterations', 100);

    % 初始化搜索范围
    m_initial_guess = max(d^2 / C_r, 1); % 确保正值
    m_opt = fsolve(objective_function, m_initial_guess, options);

    % 约束边界修正
    if m_opt < 1
        m_opt = 1;
    elseif m_opt > 500
        m_opt = 500;
    end
end


