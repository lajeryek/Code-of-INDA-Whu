%% ========================================================================
function eps = epsilon_fbl_awgn(n, d, gamma)
% 有限块长 AWGN 误码近似
% n : 块长
% d : 比特数，允许 r = d/n > 1
% gamma : SNR

if gamma <= 0 || n <= 0 || d < 0
    eps = 1;
    return;
end

r = d / n;                % 速率
C = C_awgn(gamma);
V = V_awgn(gamma);

arg = sqrt(n ./ V) .* (C - r) * log(2);
eps = qfunc_local(arg);

% 数值截断到 [0,1]
eps = max(min(eps, 1), 0);
end




