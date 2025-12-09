clc;
clear;

% 参数初始化
T_S = 1e-5;           % 符号长度
SNR = 3;              % 单个用户的信噪比
d = 100;              % 数据量

% m 的取值范围
m_range = 1:400;

% 初始化 AOI 数组
AOI_values = zeros(1, length(m_range));  % 存储 AOI 值

% 遍历 m 的值
for i = 1:length(m_range)
    m = m_range(i);
    
    % 计算该 m 下的错误率
    error_value = erro(m, SNR, d);
    
    % 计算 AOI 值
    AOI = 0.5 * m * T_S + (m * T_S) / (1 - error_value);
    
    % 存储计算得到的 AOI 值
    AOI_values(i) = AOI;
end

% 绘制 AOI 和 m 的关系曲线
figure;
plot(m_range, AOI_values, 'LineWidth', 2);
xlabel('m');
ylabel('AOI');
title('AOI 与 m 的关系');
grid on;

% erro函数定义
function [e] = erro(m_n, SNR_normal, d)
    V_r = 1 - (1 + SNR_normal).^(-2);  % 信道色散
    C_r = log2(1 + SNR_normal);        % 香农容量
    e = qfunc(((m_n / V_r)^(1/2)) * (C_r - d./m_n) * log(2));  % 误差函数
    if m_n == 0
        e = 1;
    end
end
