%% ======================== 基本功能说明 =============================
% 本代码用于分析 Alice → Bob 的短包传输系统中，在有限块长度与 HARQ 重传机制下：
%   1) Bob 解码失败概率 e_tot_B
%   2) Eve 窃听失败概率 e_tot_E
%   3) 泄露失败总概率 e_tot_LF（Leakage Failure）
%   4) Bob–Alice 总能耗 E_tot
%
% 并针对不同的反馈延时 tk 以及重传次数 N，绘制 e_tot_LF 随块长 m 变化曲线。
%
% 关键特征包括：
%   - 有限块长误码模型 erro(m, SNR, d)
%   - (N+1) 次最多重传（含初始传输）
%   - HARQ 中 ACK/NACK 的时延开销 tk
%   - 发送能量、反馈能量与计算能量的综合能耗模型
%   - 满足时延约束、能耗约束与 CPU 频率约束的传输策略筛选

%% ======================== 基本参数定义 =============================

d = 320;                 % 发送数据量（比特）
B = 5e6;                 % 带宽 (Hz)
Ts = 3e-5;               % 单符号持续时间 (s)
T = 60e-3;               % 整个传输帧长 (s)
f_cpu_max = 3.5e9;       % 最大 CPU 主频
PK = 30;                 % Alice 发射功率 (dBm)
k = 10^(-11);            % 计算能耗相关常数

% ------------ 信道参数（Bob 信道更好）------------
Ns_Bob = -164;           % Bob 噪声功率 (dBm/Hz)
Ns_Eve = -163;           % Eve 噪声功率
r_Bob = 135;             % Alice–Bob 距离
r_Eve = 140;             % Alice–Eve 距离
z_Bob = 2.5428;          % Bob 信道增益
z_Eve = 2.4428;          % Eve 信道增益

% 路径损耗模型 PL = 17 + 40·log10(r)
PL_Bob = 17.0 + 40*log10(r_Bob);
PL_Eve = 17.0 + 40*log10(r_Eve);

% Bob/Eve 接收信号功率 (dB)
P_receive_dB_Bob = PK - PL_Bob - (Ns_Bob + 10*log10(B));
P_receive_dB_Eve = PK - PL_Eve - (Ns_Eve + 10*log10(B));

% 转换为线性 SNR = 10^((dBm–30)/10)
SNR_PL_Bob = 10.^((P_receive_dB_Bob - 30)/10);
SNR_PL_Eve = 10.^((P_receive_dB_Eve - 30)/10);

% 平均 SNR（包含信道增益）
SNR_normal_Bob = (abs(z_Bob).^2) .* SNR_PL_Bob;
SNR_normal_Eve = (abs(z_Eve).^2) .* SNR_PL_Eve;

v = 1e-8;                % ACK/NACK 反馈错误概率（极小）
cnt = 1;                 % 记录数组索引

% 不同场景下的 ACK/NACK 延时 tk 与最大重传次数 N
tk_array = [0, 3e-3, 5e-3];
N_value  = [2, 3, 4];

e_tot_min = Inf;         % 记录最小泄露失败率（若需后续优化）
m_min = Inf;             % 最优块长（可扩展用途）
N_min = Inf;             % 最优重传次数

figure;

%% ===================== 主循环：遍历不同 tk 与 N ====================
for idx = 1 : length(N_value)

    tk = tk_array(idx);     % 当前场景下的反馈时延
    N  = N_value(idx);      % 最大重传次数
    
    % 预分配数组（长度 = 最大可行块长 floor(T/Ts)）
    e_tot_B_array = NaN(1, floor(T/Ts));
    e_tot_E_array = NaN(1, floor(T/Ts));
    e_tot_LF_array = NaN(1, floor(T/Ts));
    E_tot_array = NaN(1, floor(T/Ts));
    m_array = 1:floor(T/Ts);

    %% ========== 对不同块长 m = t/Ts 进行遍历 ==========
    for m = 1 : floor(T/Ts)
        t = Ts * m;     % 实际传输时间 t
        
        % 有限块长误码模型
        e_B = erro(m, SNR_normal_Bob, d);   % Bob 解码失败率
        e_E = erro(m, SNR_normal_Eve, d);   % Eve 解码失败率
        e_LF = (e_B .* e_E) + (1 - e_E);    % 泄露失败率（一次传输）

        % ----------------- 计算 Bob 累计解码失败概率 -----------------
        e_tot_B = 0;
        for n = 1:N
            e_tot_B = e_tot_B + (e_B^n) * ((1 - v)^(n-1)) * v;
        end
        % N 次都失败
        e_tot_B = e_tot_B + (e_B^(N+1)) * ((1 - v)^N);

        % ----------------- 计算 Eve 累计窃听失败概率 -----------------
        e_tot_E = e_E;
        for n = 1:N
            e_tot_E = e_E * (e_B * ((1-v) * e_tot_E + v) + (1 - e_B));
        end

        % ----------------- 最终泄露失败概率 -----------------
        e_tot_LF = e_tot_B .* e_tot_E + (1 - e_tot_E);

        % ----------------- 计算发送/反馈/计算能耗 -----------------
        Et0 = t * (10.^((PK-30)/10));   % 对应一次发送能量
        Et = Et0;
        for n = 0:N
            Et = Et + (e_B^n) * ((1 - v)^n) * Et0;
        end

        % 反馈能耗
        Ek0 = tk * (10.^((PK-30)/10));
        Ek = 0;
        for n = 1:N
            Ek = Ek + (e_B^(n+1)) * ((1 - v)^n) * Ek0;
        end

        % 计算能耗（CPU）
        c = 20;
        Ec0 = (k * (c^3)) / ((T - t)^2);
        Ec = 0;
        for n = 1:N
            Ec_n  = (k * c^3) / ((T - (n+1)*t - n*tk)^2);
            Ec_nn = (k * c^3) / ((T -  n   *t - (n-1)*tk)^2);
            Ec = Ec + (e_B^n) * ((1 - v)^(n-1)) * (Ec_n - Ec_nn);
        end
        Ec = Ec0 + Ec;

        % Bob 端总能耗（不考虑 Eve）
        E_tot = Et + Ec + Ek;

        % ------------------- 可行性约束检查 -------------------
        % (1) 计算时延 ≤ T
        % (2) HARQ 设计满足时序递推
        % (3) 能耗 < 0.5（自定义约束）
        if ((20e6 / f_cpu_max + (N+1)*t - N*tk) <= T) && ...
           (N <= floor((T - (20e6/f_cpu_max) - t) / (t + tk))) && ...
           (E_tot < 0.5)

            e_tot_LF_array(cnt) = e_tot_LF;
            e_tot_B_array(cnt) = e_tot_B;
            e_tot_E_array(cnt) = e_tot_E;
            E_tot_array(cnt)   = E_tot;
            m_array(cnt)       = m;
            cnt = cnt + 1;
        end
    end

    % 绘制泄露失败概率 vs 块长 m
    plot(m_array, e_tot_LF_array);
    hold on;
end
