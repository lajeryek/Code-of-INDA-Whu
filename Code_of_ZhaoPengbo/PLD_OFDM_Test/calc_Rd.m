function Rd = calc_Rd(nM, dM, PM, hB, nK, dK, PK, hE)
% calc_Rd(...)
% Compute four BLERs and deception rate Rd for OFDM (orthogonal) PLD model.
%
% Bob-M: epsilon(nM, dM, gamma=zB*PM/sigma2)
% Bob-K: epsilon(nK, dK, gamma=zB*PK/sigma2)
% Eve-M: epsilon(nM, dM, gamma=zE*PM/sigma2)
% Eve-K: epsilon(nK, dK, gamma=zE*PK/sigma2)
%
% Outputs:
%   eBM, eBK, eEM, eEK: four BLERs
%   Rd: deception rate
%
% Default sigma2 = 1.

sigma2 = 1;

zB = hB.^2;
zE = hE.^2;

gBM = zB .* PM ./ sigma2;
gBK = zB .* PK ./ sigma2;
gEM = zE .* PM ./ sigma2;
gEK = zE .* PK ./ sigma2;

eBM = epsilon_fbl_awgn(nM, dM, gBM);
eBK = epsilon_fbl_awgn(nK, dK, gBK);
eEM = epsilon_fbl_awgn(nM, dM, gEM);
eEK = epsilon_fbl_awgn(nK, dK, gEK);

Rd = Rd_metric(eBM, eBK, eEM, eEK);

fprintf('==== calc_Rd (OFDM) ====\n');
fprintf('Bob: eps=%.6e | hB=%g | M:(n=%g,d=%g,P=%g,gamma=%g) \n', eBM, hB, nM, dM, PM, gBM);
fprintf('Bob: eps=%.6e | hB=%g | K:(n=%g,d=%g,P=%g,gamma=%g) \n', eBK, hB, nK, dK, PK, gBK);
fprintf('Eve: eps=%.6e | hE=%g | M:(n=%g,d=%g,P=%g,gamma=%g) \n', eEM, hE, nM, dM, PM, gEM);
fprintf('Eve: eps=%.6e | hE=%g | K:(n=%g,d=%g,P=%g,gamma=%g) \n', eEK, hE, nK, dK, PK, gEK);
fprintf('Rd = %.12f\n', Rd);
fprintf('========================\n');

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

function rd = Rd_metric(eBM, eBK, eEM, eEK)
    rd = (1 - (1 - eBM)*eBK) * (1 - eEM) * eEK;
end
