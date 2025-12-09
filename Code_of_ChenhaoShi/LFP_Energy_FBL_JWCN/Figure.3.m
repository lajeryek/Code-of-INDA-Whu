%d = 320; % Alice向服务器Bob传输的数据量为320bit
B = 5e6; % 带宽Hz
Ts = 3e-5; % 符号长度 单位s
T = 60e-3; % 帧长 单位s
f_cpu_max = 3.5e9; % 最大CPU频率
 %PK = 30; % 发射功率dBm
 %d_array=[320,420,520];
%PK_array=[20,25,30];
k = 10^(-11);
% 描述两个接收端的信道质量，认为Bob好于窃听的Eve,分别从传输距离r，噪声大小Ns以及信道增益z考量
Ns_Bob=-165;%Bob的噪声大小
Ns_Eve=-162;%Eve的噪声大小
r_Bob=135;%传输距离
r_Eve=145;%传输距离
 z_Bob=2.5428;
 z_Eve=2.4428;
 
v = 1e-5;%UE解码失误概率
cnt = 1;%数组索引
tk = 3e-3;%发送NACK的时间
N = 2; % Bob最大重传次数
e_tot_max = 1e-4;%可靠性约束阈值
E_tot_min = Inf;
m_min = Inf;
figure; % 创建一个新的图形窗口
  for idx=1:1:5
    if idx==1
        PK=30;
        d=320;
    elseif idx==2
        PK=35;
        d=320;
    elseif idx==3
        PK=40;
        d=320;
    elseif idx==4
        PK=35;
        d=520;
    elseif idx==5
        PK=40;
        d=520;
    end
    PL_Bob=17.0+40*log10(r_Bob);%路径损耗Pathloss_dBm_Bob
    PL_Eve=17.0+40*log10(r_Eve);%路径损耗Pathloss_dBm_Eve
    P_receive_dB_Bob = PK - PL_Bob - (Ns_Bob + 10*log10(B)); %Bob接收功率 
    P_receive_dB_Eve = PK - PL_Eve - (Ns_Eve + 10*log10(B)); %Eve接收功率 

    SNR_PL_Bob=(10.^((P_receive_dB_Bob-30)/10));%信噪比 
    SNR_PL_Eve=(10.^((P_receive_dB_Eve-30)/10));%信噪比 

    SNR_normal_Bob=(abs(z_Bob).^2).*SNR_PL_Bob;%平均信噪比
    SNR_normal_Eve=(abs(z_Eve).^2).*SNR_PL_Eve;%平均信噪比
    e_tot_B_array = NaN(1, floor(T/Ts)); % 预先定义长度，并填充为nan
    e_tot_E_array = NaN(1, floor(T/Ts)); % 预先定义长度，并填充为nan
    e_tot_LF_array = NaN(1, floor(T/Ts)); % 预先定义长度，并填充为nan
    m_array = 1:T/Ts;
    E_tot_array = NaN(1, floor(T/Ts));
    E_tot_min_array=[];
    for m = 1:floor(T/Ts)
        t = Ts * m;
        e_B = erro(m, SNR_normal_Bob,d); % Bob的解码错误概率
        e_E = erro(m, SNR_normal_Eve,d); % Eve窃听失败的概率
        e_LF = (e_B.*e_E) + (1 - e_E); % Eve窃听成功或者Bob解码失败的概率
        % 计算三个部分的e_tot
        e_tot_B = 0;
        e_tot_E = e_E; % 初始传输的错误概率
     
        % 计算Bob_erro
        for n = 1:1:N
            e_tot_B = e_tot_B + (e_B^n) * ((1 - v)^(n-1)) * v; % Bob总的错误率
        end
        e_tot_B = e_tot_B + ((e_B.^(N+1)) .* ((1-v).^N)); % 计算Bob总的错误率
        
        %计算Eve——error
         for n=1:1:N
            e_tot_E= e_E * (e_B * ((1 - v) * e_tot_E + v * 1) + (1 - e_B));
         end
         
        % 计算ELF的总概率
         e_tot_LF=(e_tot_B.*e_tot_E)+(1-e_tot_E);%1-(1-e_B)*e_E

        % 计算能量
        Et0 = t .* (10.^((PK-30)./10));
        Et = Et0;
        for n = 0:1:N
            Et = Et + ((e_B.^n) .* ((1-v)^n)) .* Et0;
        end

        Ek0 = tk * (10.^((PK-30)./10));
        Ek = 0;
        for n = 1:1:N
            Ek = Ek + ((e_B.^(n+1)) * ((1-v)^n)) .* Ek0;
        end

        c = 20;
        Ec0 = ((k * (c^3)) ./ ((T-t).^2));
        Ec = 0;
        for n = 1:1:N
            Ec_n = ((k * (c^3)) ./ ((T-(n+1).*t-n.*(tk)).^2)); % Ec(n)
            Ec_nn = ((k * (c^3)) ./ ((T-n.*t-(n-1).*(tk)).^2)); % Ec(n-1)
            Ec = Ec + ((e_B.^n) .* ((1-v).^(n-1))) .* (Ec_n - Ec_nn);
        end
        Ec = Ec0 + Ec;
        E_tot = (Et + Ec + Ek); %考虑Bob与Alice之间的总能耗，不考虑Eve的能耗

        if  ((20e6 / f_cpu_max + (N+1) * t - N *tk) <= T) &&(N<=floor((T-(20e6/f_cpu_max)-t)/(t+tk)))&&(e_tot_LF<=0.01)%&& (e_tot_B <= e_tot_max)&&(e_tot_E>=0.1)
            e_tot_LF_array(cnt) = e_tot_LF;
            E_tot_array(cnt) = E_tot;
            m_array(cnt) = m;
            cnt = cnt + 1;
        end
    end
     plot(m_array, E_tot_array);
     hold on
  end
 
