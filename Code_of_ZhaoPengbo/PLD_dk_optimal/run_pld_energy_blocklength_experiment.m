%% ========================================================================
function results = run_pld_energy_blocklength_experiment(params, nM_list, PM_list, PK_list)
% 扫描 (n_M, P_M, P_K)，对每个组合在 d_K 上做一维搜索，返回 R_d^max

NnM = numel(nM_list);
NPM = numel(PM_list);
NPK = numel(PK_list);

Rd_best = nan(NnM, NPM, NPK);
dK_best = nan(NnM, NPM, NPK);

best_global.Rd = NaN;
best_global.dK = NaN;
best_global.nM = NaN;
best_global.nK = NaN;
best_global.PM = NaN;
best_global.PK = NaN;

for i = 1:NnM
    nM = nM_list(i);
    nK = params.Ntot - nM;
    if nK <= 0
        continue; % 块长无效
    end

    for j = 1:NPM
        PM = PM_list(j);
        if PM > params.Pmax
            continue;
        end

        for k = 1:NPK
            PK = PK_list(k);
            if PK > params.Pmax
                continue;
            end

            % 约束 1: 能量约束
            if PM * nM + PK * nK > params.Etot
                continue;
            end

            % 对给定 (PM,PK,nM,nK) 在 d_K 上一维搜索
            [Rd_max_here, dK_opt_here] = optimize_dK_for_given_resources( ...
                params, PM, PK, nM, nK);

            Rd_best(i,j,k) = Rd_max_here;
            dK_best(i,j,k) = dK_opt_here;

            % 更新全局最优
            if ~isnan(Rd_max_here) && (isnan(best_global.Rd) || Rd_max_here > best_global.Rd)
                best_global.Rd  = Rd_max_here;
                best_global.dK  = dK_opt_here;
                best_global.nM  = nM;
                best_global.nK  = nK;
                best_global.PM  = PM;
                best_global.PK  = PK;
            end
        end
    end
end

results.Rd_best = Rd_best;
results.dK_best = dK_best;
results.best    = best_global;
results.nM_list = nM_list;
results.PM_list = PM_list;
results.PK_list = PK_list;
end