clc
clear

% 参数初始化
k = 1;         % 迭代参数
error_array = [];
AOI_array = [];
T_S = 3e-5;    % 符号长度
I = 2;         % 传感器数量
y_old = 0;     % 初始吞吐量

% 用户信道状态
SNR_vector = [10, 12];    % 每个用户的信噪比
d_vector = [140, 120];  % 每个用户的数据量
m_vector = randi(300, 1, I);  % 初始化传输块长，随机值

% 约束条件
AOI_th = 5;    % 设定的平均AOI阈值
T_max = 50e-3; % 最大时延
error_max = 1e-3; % 最大错误率

% 收敛条件的参数
tolerance = 1e-2; % 收敛条件：吞吐量的变化小于此阈值
max_iters = 10000; % 最大迭代次数

% 二分法的上下限
m_lower = 1 * ones(1, I);  % 每个用户的传输块长最小值
m_upper = 400 * ones(1, I);   % 每个用户的传输块长最大值

% 存储每次迭代的吞吐量值
y_values = [];

% Dinkelbach's Transform 迭代
while k <= max_iters
    m_vector_old = m_vector;  % 记录上一轮的 m_i
    A_mi = 0;  % 总的 A_mi
    B_mi = 0;  % 总的 B_mi
    
    % 计算当前的错误率、AOI、A_mi 和 B_mi
    for i = 1:I
        m_i = m_vector(i);
        d_i = d_vector(i);
        SNR_i = SNR_vector(i);
        error_array(i) = erro(m_i, SNR_i, d_i);  % 计算错误率
        AOI_array(i) = 0.5 * m_i * T_S + (m_i * T_S) / (1 - error_array(i));  % 计算 AOI
        A_mi = A_mi + (1 - error_array(i)) * d_i;  % 累加 A_mi
        B_mi = B_mi + m_i * T_S;  % 累加 B_mi
    end
    
    % 判断约束条件
    if all(AOI_array <= AOI_th) && B_mi <= T_max && all(error_array <= error_max)
        % 更新 y：根据上一轮的 m_i 更新 y
        y_new = A_mi / B_mi;

        % 记录当前迭代的吞吐量
        y_values(k) = y_new;

        % 对每个用户使用二分法更新 m_i
        for i = 1:I
            % 保存当前用户的 m_i
            m_i = m_vector(i);
            d_i = d_vector(i);
            SNR_i = SNR_vector(i);
            
            % 二分法迭代更新 m_i
            while (m_upper(i) - m_lower(i)) > tolerance
                m_mid = (m_upper(i) + m_lower(i)) / 2;  % 二分法中点
                
                % 计算中点的错误率和平均 AOI
                error_mid = erro(m_mid, SNR_i, d_i);  % 计算中点处的错误率
                AOI_mid = 0.5 * m_mid * T_S + (m_mid * T_S) / (1 - error_mid);  % 计算 AOI
                
                % 计算总的 A_mi 和 B_mi
                A_total = 0;
                B_total = 0;
                for j = 1:I
                    if j ~= i  % 排除当前用户
                        error_j = erro(m_vector(j), SNR_vector(j), d_vector(j));  % 计算其他用户的错误率
                        A_total = A_total + d_vector(j) * (1 - error_j);  % 累加 A
                        B_total = B_total + m_vector(j) * T_S;  % 累加 B
                    end
                end
                
                % 计算中点的目标函数
                A_mid = A_total + d_i * (1 - error_mid);
                B_mid = B_total + m_mid * T_S;  % 当前 m_i 也需要计算在内
                f_mid = A_mid - y_new * B_mid;  % 目标函数

                % 计算当前传输块长的目标函数
                error_i = erro(m_i, SNR_i, d_i);
                A_i = A_total + d_i * (1 - error_i);
                B_i = B_total + m_i * T_S;  % 再加上当前用户
                f_i = A_i - y_new * B_i;

                % 根据目标函数 f(m) 来调整 m_i
                if f_mid >= f_i
                    m_lower(i) = m_mid;  % 增大 m_i
                else
                    m_upper(i) = m_mid;  % 减小 m_i
                end
            end
            m_vector(i) = (m_upper(i) + m_lower(i)) / 2;  % 更新 m_i
        end

        % 判断吞吐量是否收敛
        if abs(y_new - y_old) < tolerance
            break;
        end
        
        % 更新 y_old 和迭代次数增加
        y_old = y_new;
        k = k + 1;
    else
        % 如果不满足约束条件，调整对应的 m_i 到约束的边界值
        for i = 1:I
            if AOI_array(i) > AOI_th
                m_vector(i) = AOI_th / (0.5 * T_S + T_S / (1 - error_array(i)));  % 满足 AOI 的临界 m
            end
            
            if error_array(i) > error_max
                % 计算 C_r 和 V_r
                C_r = log2(1 + SNR_vector(i));
                V_r = 1 - (1 + SNR_vector(i))^(-2);
                z = erfcinv(error_max);  % 逆 Q 函数

                % 使用新的公式计算边界 m
                discriminant = (2 * C_r * d_vector(i))^2 - 4 * C_r^2 * (d_vector(i)^2 - (z^2 * V_r) / (log(2)^2));

                if discriminant >= 0
                    m_boundary_error_pos = (2 * C_r * d_vector(i) + sqrt(discriminant)) / (2 * C_r^2);
                    m_boundary_error_neg = (2 * C_r * d_vector(i) - sqrt(discriminant)) / (2 * C_r^2);

                    % 选择合理的解
                    m_vector(i) = max(1, min(400, max(m_boundary_error_pos, m_boundary_error_neg)));
                else
                    % 处理判别式为负的情况（无有效解）
                    disp(['警告: 计算得到的 m 超出范围，已调整到合理范围。用户 ', num2str(i)]);
                    m_vector(i) = 400; % 赋一个默认值或相应处理
                end
            end
            
            if B_mi > T_max
                m_vector(i) = (T_max - B_total + m_vector(i) * T_S) / T_S;  % 满足时延的临界 m
            end
        end
        disp('当前迭代未满足约束条件，调整 m_i。');
    end
end

% 最终输出结果时判断约束条件
if all(AOI_array <= AOI_th) && B_mi <= T_max && all(error_array <= error_max)
    disp(['R_max = ', num2str(y_new)]);
    disp('满足所有约束条件');
else
    disp(['R_max = ', num2str(y_new)]);
    disp('未满足所有约束条件');
end

% 返回最优传输块长
disp('Optimal m_vector = ');
disp(m_vector);

% 绘制吞吐量收敛曲线
figure;
plot(1:length(y_values), y_values, 'LineWidth', 2);
xlabel('迭代次数');
ylabel('最大吞吐量');
grid on;

% 定义Q函数的互补误差函数
qfunc = @(x) 0.5 * erfc(x / sqrt(2));

% erro函数计算错误率
function [e] = erro(m_n, SNR_normal, d)
    V_r = 1 - (1 + SNR_normal).^(-2);  % 信道色散
    C_r = log2(1 + SNR_normal);        % 香农容量
    e = qfunc(((m_n / V_r)^(1/2)) * (C_r - d / m_n));  % 错误率公式
    if m_n == 0
        e = 1;
    end
end
