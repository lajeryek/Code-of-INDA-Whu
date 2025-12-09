% Function to calculate total energy considering channel errors and decoding failures
function [E_tot_min_array] = calculateEnergy(T_array, d, SNR_normal_Bob, SNR_normal_Eve, v, tk, k, f_cpu_max, PK, Ts)
    E_tot_min_array = [];
    cnt = 1;
    for T = T_array
        E_tot_min = Inf;
        e_tot_B_array = NaN(1, floor(T/Ts)); 
        e_tot_E_array = NaN(1, floor(T/Ts)); 
        e_tot_LF_array = NaN(1, floor(T/Ts)); 
        m_array = 1:T/Ts;
        E_tot_array = NaN(1, floor(T/Ts));
        for N = 0:1:9
            for m = 1:floor(T/Ts)
                t = Ts * m;
                e_B = erro(m, SNR_normal_Bob, d); 
                e_E = erro(m, SNR_normal_Eve, d); 
                e_LF = (e_B .* e_E) + (1 - e_E); 

                e_tot_B = 0;
                e_tot_E = e_E; 

                for n = 0:1:N
                    e_tot_B = e_tot_B + (e_B^n) * ((1 - v)^(n-1)) * v;
                end
                e_tot_B = e_tot_B + ((e_B.^(N+1)) .* ((1-v).^N));

                for n = 0:1:N
                    e_tot_E = e_E * (e_B * ((1 - v) * e_tot_E + v * 1) + (1 - e_B));
                end

                e_tot_LF = (e_tot_B .* e_tot_E) + (1 - e_tot_E);

                Et0 = t .* (10.^((PK-30)./10));
                Et = Et0;
                for n = 0:1:N
                    Et = Et + ((e_B.^n) .* ((1-v)^n)) .* Et0;
                end

                Ek0 = tk * (10.^((PK-30)./10));
                Ek = 0;
                for n = 0:1:N
                    Ek = Ek + ((e_B.^(n+1)) * ((1-v)^n)) .* Ek0;
                end

                c = 350;
                Ec0 = ((k * (c^3)) ./ ((T-t).^2));
                Ec = 0;
                for n = 0:1:N
                    Ec_n = ((k * (c^3)) ./ ((T-(n+1).*t-n.*(tk)).^2)); 
                    Ec_nn = ((k * (c^3)) ./ ((T-n.*t-(n-1).*(tk)).^2)); 
                    Ec = Ec + ((e_B.^n) .* ((1-v).^(n-1))) .* (Ec_n - Ec_nn);
                end
                Ec = Ec0 + Ec;
                E_tot = (Et + Ec + Ek); 

                if  ((20e6 / f_cpu_max + (N+1) * t - N * tk) <= T) && (N <= floor((T-(20e6/f_cpu_max)-t)/(t+tk))) && (e_tot_LF <= 0.0001)
                    e_tot_LF_array(cnt) = e_tot_LF;
                    E_tot_array(cnt) = E_tot;
                    m_array(cnt) = m;
                    if E_tot <= E_tot_min
                        E_tot_min = E_tot;
                    end
                    cnt = cnt + 1;
                end
            end
        end
        E_tot_min_array = [E_tot_min_array, E_tot_min];
    end
end

B = 5e6; 
Ts = 2.5e-5; 
T_array = [40e-3, 50e-3, 60e-3, 70e-3, 80e-3];
f_cpu_max = 3.9e9; 
PK = 30;
k = 10^(-11);
SNR_normal_Bob = 2;
SNR_normal_Eve = 1;
v = 1e-5;
cnt = 1;
tk = 3e-3;
figure;

% Call the function for different values of d
for idx = 1:1:3
    if idx == 1
        d = 240;
    elseif idx == 2
        d = 320;
    elseif idx == 3
        d = 520;
    end
    E_tot_min_array = calculateEnergy(T_array, d, SNR_normal_Bob, SNR_normal_Eve, v, tk, k, f_cpu_max, PK, Ts);
    plot(T_array, E_tot_min_array, '-o');
    hold on
end

% Special cases for d=240 and d=320
d=240;
E_tot_min_array = calculateEnergy(T_array, d, SNR_normal_Bob, SNR_normal_Eve, v, tk, k, f_cpu_max, PK, Ts);
plot(T_array, E_tot_min_array, '-o');
hold on

d=320;
E_tot_min_array = calculateEnergy(T_array, d, SNR_normal_Bob, SNR_normal_Eve, v, tk, k, f_cpu_max, PK, Ts);
plot(T_array, E_tot_min_array, '-o');
hold on

