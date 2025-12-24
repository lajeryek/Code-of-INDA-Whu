function eps = calc_epsilon(n, d, P, h)
% calc_epsilon(n,d,P,h)
% Compute finite-blocklength BLER: eps = Q( ln2*sqrt(n/V)*(C-d/n) )
% gamma = (h^2 * P) / sigma2, default sigma2 = 1
%
% Input:
%   n: blocklength
%   d: number of bits (can be > n)
%   P: transmit power
%   h: channel amplitude (use z=h^2)
% Output:
%   eps: block error probability in [0,1]

sigma2 = 1;
z = h.^2;
gamma = z .* P ./ sigma2;

eps = epsilon_fbl_awgn(n, d, gamma);

fprintf('[calc_epsilon] n=%g, d=%g, P=%g, h=%g -> gamma=%g, eps=%.6e\n', ...
    n, d, P, h, gamma, eps);

end

%% ===== local helper functions =====
function epsv = epsilon_fbl_awgn(n, d, gamma)
    if gamma <= 0 || n <= 0 || d < 0
        epsv = 1; return;
    end
    C = log2(1 + gamma);
    V = 1 - 1/(1 + gamma)^2;
    X = log(2) * sqrt(n / V) * (C - d/n);
    epsv = qfunc_local(X);
    epsv = min(max(epsv, 0), 1);
end

function y = qfunc_local(x)
    y = 0.5 * erfc(x / sqrt(2));
end
