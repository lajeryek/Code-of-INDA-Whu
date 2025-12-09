for idx=1:1:3
   if idx==1
    N=2;
    e_tot_max_array=[3e-3,5e-3,7e-3,8.5e-3,1e-2];
   end
   if idx==2
   N=3;
    e_tot_max_array=[1e-8,1e-5,1e-4,0.5e-3,1e-3,1.5e-3,2e-3,3e-3];
   end
   E_tot_min = Inf;
    e_tot_min=Inf;
     m_min=Inf;
     E_tot_min_array=[];
     m_min_array=[];
     for e_tot_max=e_tot_max_array
            SNR_normal_Bob=38.4;
            SNR_normal_Eve=30;
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
                 for n=1:1:N
                    e_tot_E= e_E * (e_B * ((1 - v) * e_tot_E + v * 1) + (1 - e_B));
                 end
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

                if  ((2e7 / f_cpu_max + (N+1) * t - N *tk) <= T) &&(N<=floor((T-(2e7/f_cpu_max)-t)/(t+tk)))&&(e_tot_LF<=e_tot_max)%&& (e_tot_B <= e_tot_max)&&(e_tot_E>=0.1)
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
     end

