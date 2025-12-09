%% ========================================================================
function [Rd_max, dK_opt] = optimize_dK_for_given_resources(params, PM, PK, nM, nK)
% 对固定 (PM,PK,nM,nK)，在 d_K ∈ [0, 2 nK] 上一维搜索 R_d 最大值
% 若在所有 d_K 下都不满足阈值约束，则返回 NaN

sigma2 = params.sigma2;
zB = params.zB; zE = params.zE;

% OFDM 下的 SNR（M/K 正交，无干扰）
gBM = zB * PM / sigma2;
gBK = zB * PK / sigma2;
gEM = zE * PM / sigma2;
gEK = zE * PK / sigma2;

% 与 d_K 无关的误码
eBM = epsilon_fbl_awgn(nM, params.dM, gBM);
eEM = epsilon_fbl_awgn(nM, params.dM, gEM);

% 若这两个本身就违反阈值，则任何 d_K 都不可行，直接返回 NaN
if eBM > params.epsBM_th || eEM > params.epsEM_th
    Rd_max = NaN; dK_opt = NaN;
    return;
end

% d_K 网格（利用凹性，一维搜索可以更粗一点）
dK_vec = linspace(0, 2*nK, params.dK_grid_size);
Rd_vec = nan(size(dK_vec));

for idx = 1:numel(dK_vec)
    dK = dK_vec(idx);

    eBK = epsilon_fbl_awgn(nK, dK, gBK);
    eEK = epsilon_fbl_awgn(nK, dK, gEK);

    % 阈值约束
    if eBK <= params.epsBK_th && ...
       eEM <= params.epsEM_th && ...
       eEK >= params.epsEK_min

        Rd_vec(idx) = Rd_metric(eBM, eBK, eEM, eEK);
    else
        Rd_vec(idx) = NaN;
    end
end

if all(isnan(Rd_vec))
    Rd_max = NaN; dK_opt = NaN;
else
    [Rd_max, idx_best] = max(Rd_vec, [], 'omitnan');
    dK_opt = dK_vec(idx_best);
end
end