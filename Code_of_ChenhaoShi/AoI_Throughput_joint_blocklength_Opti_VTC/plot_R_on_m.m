clc;
clear;

% 参数初始化
T_S = 1e-5;      % 符号长度
I = 2;           % 传感器数量
T_max = 50e-3;   % 最大时延
error_max = 1e-2; % 最大错误率

% 用户1的参数
SNR_1 = 3;       % 用户1的信噪比
d_1 = 150;        % 用户1的数据量

% 用户2的参数
SNR_2 = 8;       % 用户2的信噪比
d_2 = 130;         % 用户2的数据量

AOI_th = 0.5;      % 设定的平均AOI阈值

% 传输块长 m_i 的取值范围
m_min = 1;       % 最小块长
m_max = 400;     % 最大块长
m_step = 1;      % 块长步长
m_values = m_min:m_step:m_max;  % 块长遍历范围

% 初始化存储系统吞吐量 R
R_values = zeros(length(m_values), length(m_values));  % 存储不同 (m_1, m_2) 的吞吐量

% 遍历不同的 m_1 和 m_2 值
for index1 = 1:length(m_values)
    m_1 = m_values(index1);  % 当前用户1的块长
    for index2 = 1:length(m_values)
        m_2 = m_values(index2);  % 当前用户2的块长
        
        % 计算用户1的错误率
        error_1 = erro(m_1, SNR_1, d_1);
        % 计算用户2的错误率
        error_2 = erro(m_2, SNR_2, d_2);

        % 计算每个用户的平均AOI
        AOI_1 = 0.5 * m_1 * T_S + (m_1 * T_S) / (1 - error_1);
        AOI_2 = 0.5 * m_2 * T_S + (m_2 * T_S) / (1 - error_2);
        
        % 计算 A(m_1) 和 B(m_1)
        A_1 = d_1 * (1 - error_1);  % A(m_1)
        B_1 = m_1 * T_S;             % B(m_1)

        % 计算 A(m_2) 和 B(m_2)
        A_2 = d_2 * (1 - error_2);  % A(m_2)
        B_2 = m_2 * T_S;             % B(m_2)

        % 判断是否满足约束条件
        if AOI_1 <= AOI_th && AOI_2 <= AOI_th && (B_1 + B_2) <= T_max
            R_values(index1, index2) = (A_1 + A_2) / (B_1 + B_2);  % 计算系统吞吐量 R
        else
            R_values(index1, index2) = NaN;  % 不满足条件时标记为 NaN
        end
    end
end

% 绘制 R 关于 m_1 和 m_2 的关系图
figure;
mesh(m_values, m_values, R_values);
xlabel('Blocklength m_1 (symbols)');
ylabel('Blocklength m_2 (symbols)');
zlabel('Effective Sum Throughput R');
%title('系统吞吐量 R 与传输块长 m_1 和 m_2 的关系');
grid on;
