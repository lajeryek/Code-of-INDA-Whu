cnt = 1;%数组索引
N_array=[1,2,3,4,5,6,7];
 figure;
set(gcf,'color','w');
set(gca,'FontSize',14);
 for idx=1:1:3
    if idx==1
        tk = 2.5e-3;
    end
    if idx==2
        tk = 7.25e-3;
    end
    if idx==3
       tk = 0;
    end
   % tk=3e-3;
    E_tot_min_array=[];
    m_min_array=[];
    e_tot_min_array=[];
    for N=N_array
        E_tot_min = Inf; 
        e_tot_min=Inf;
        m_min=Inf;
        SNR_normal_Bob=22.1;
        SNR_normal_Eve=17.1;
        e_tot_B_array = NaN(1, floor(T/Ts)); % 预先定义长度，并填充为nan
        e_tot_E_array = NaN(1, floor(T/Ts)); % 预先定义长度，并填充为nan
        e_tot_LF_array = NaN(1, floor(T/Ts)); % 预先定义长度，并填充为nan
        m_array = 1:(T/Ts);
        E_tot_array = NaN(1, floor(T/Ts));
         
        for m = 1:0.1:floor(T/Ts)
            t = Ts .* m;
            e_B = erro(m, SNR_normal_Bob,d); % Bob的解码错误概率
            e_E = erro(m, SNR_normal_Eve,d); % Eve窃听失败的概率
            % 初始传输的错误概率
            e_tot_B = 0;
            e_tot_E = e_E;
            % 计算Bob_erro
            for n = 1:1:N
                e_tot_B = e_tot_B + (e_B.^n) * ((1 - v)^(n-1)) * v; % Bob总的错误率
            end
            e_tot_B = e_tot_B + ((e_B.^(N+1)) .* ((1-v).^N)); % 计算Bob总的错误率

            %计算Eve——error
              for n=1:1:(N-1)
                 e_tot_E=e_tot_E+((e_E^(N-n)).*(e_B^(n)).*(1-v)^(n)).*v;
              end
            e_tot_E=e_tot_E+((e_E.^(N+1)) .* ((1-v).^N).*(e_B^(N)));
       


            % 计算LFP的总概率
              e_tot_LF=(e_tot_B.*e_tot_E)+(1-e_tot_E);%1-(1-e_B)*e_E
           %  e_tot_LF=1-(1-e_tot_B)*e_tot_E;

            % 计算能量energy
            Et0 = t .* (10.^((PK-30)./10));
            Et = 0;
            for n = 0:1:N
                Et = Et + ((e_B.^n) .* ((1-v)^n)) .* Et0;
            end

            Ek0 = tk * (10.^((PK-30)./10));
            Ek =0;
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
            E_tot = Et + Ec + Ek; %考虑Bob与Alice之间的总能耗，不考虑Eve的能耗

            if  ((2e7 / f_cpu_max + (N+1) * t - N *tk) <= T) &&(N<=floor((T-(2e7/f_cpu_max)-t)/(t+tk)))&&(e_tot_LF<=0.01)%&& (e_tot_B <= e_tot_max)&&(e_tot_E>=0.1)
                e_tot_LF_array(cnt) = e_tot_LF;
                E_tot_array(cnt) = E_tot;
                m_array(cnt) = m;
                 if E_tot<=E_tot_min
                      E_tot_min=E_tot;
                      m_min=m;
                      e_tot_min=e_tot_LF;
                 end
                cnt = cnt + 1;
            end
        end
       E_tot_min_array = [E_tot_min_array, E_tot_min];
       m_min_array=[m_min_array,m_min];
       e_tot_min_array=[e_tot_min_array,e_tot_min];
    end
     plot(N_array, E_tot_min_array, line_styles{idx});
     hold on
 end
