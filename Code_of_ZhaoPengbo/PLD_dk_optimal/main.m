function main()
% main.m — OFDM-PLD 实验：max R_d(P_M,P_K,n_M,n_K,d_K) under energy & length constraints

close all; clc;

fprintf('==== PLD in OFDM: Joint (P, n, d_K) Experiment ====\n');

%% ---------- 1. 系统参数（可按需要修改） ----------
params.sigma2 = 1;          % 噪声功率
hB = 1.2; params.zB = hB^2; % Bob 信道增益 |h_B|^2
hE = 0.7; params.zE = hE^2; % Eve 信道增益 |h_E|^2

params.dM = 32;             % 明文比特数 d_M
                            % d_K 在 [0, 2 n_K] 内搜索

% 资源预算
params.Ntot   = 200;        % 总块长上界 N_tot
params.Pmax   = 5;          % 单流最大发电功率 P_max
params.Etot   = 600;        % 总能量预算 E_tot (可< Pmax*Ntot)

% 误码阈值（与 OFDM 论文一致，可修改）
params.epsBM_th = 0.5;
params.epsBK_th = 0.5;
params.epsEM_th = 0.5;
params.epsEK_min = 0.5;

% d_K 搜索分辨率
params.dK_grid_size = 400;

%% ---------- 2. 设置扫描网格 ----------
% 块长分配：n_M 网格，n_K = Ntot - n_M
nM_list = 40:5:160;   % 可自行改为更密/更宽
% 功率分配：P_M, P_K 网格
PM_list = linspace(0.5, params.Pmax, 10);
PK_list = linspace(0.5, params.Pmax, 10);

fprintf('Grid sizes: |nM|=%d, |PM|=%d, |PK|=%d\n', ...
    numel(nM_list), numel(PM_list), numel(PK_list));

%% ---------- 3. 扫描 (n_M, P_M, P_K)，内部对 d_K 一维搜索 ----------
results = run_pld_energy_blocklength_experiment(params, nM_list, PM_list, PK_list);

%% ---------- 4. 打印全局最优点 ----------
if isnan(results.best.Rd)
    fprintf('No feasible (P_M,P_K,n_M,n_K,d_K) found under given constraints.\n');
else
    fprintf('\n=== Global best under constraints ===\n');
    fprintf('  n_M   = %d, n_K   = %d\n', results.best.nM, results.best.nK);
    fprintf('  P_M   = %.3f, P_K  = %.3f\n', results.best.PM, results.best.PK);
    fprintf('  d_K^* = %.3f\n', results.best.dK);
    fprintf('  R_d^* = %.6f\n', results.best.Rd);
end

%% ---------- 5. 画一个代表性的 R_d(d_K) 曲线，展示凹性 ----------
if ~isnan(results.best.Rd)
    figure('Color','w'); hold on; grid on;
    [dK_vec, Rd_vec] = sweep_dK_for_fixed_resources(params, ...
        results.best.PM, results.best.PK, ...
        results.best.nM, results.best.nK);
    plot(dK_vec, Rd_vec, 'k-', 'LineWidth', 2);
    xlabel('d_K'); ylabel('R_d');
    title('R_d vs d_K at the globally best (P_M,P_K,n_M,n_K)');
    xline(results.best.dK, 'g--', 'd_K^*');
end

%% ---------- 6. 画一个简单的热力图（对 PK 取 max） ----------
% 聚合：对每个 (n_M, P_M) 取最大 R_d over P_K
Rd_best_2D = nan(numel(nM_list), numel(PM_list));
for i = 1:numel(nM_list)
    for j = 1:numel(PM_list)
        Rd_slice = squeeze(results.Rd_best(i,j,:));
        if all(isnan(Rd_slice)), Rd_best_2D(i,j) = NaN;
        else, Rd_best_2D(i,j) = max(Rd_slice, [], 'omitnan');
        end
    end
end

figure('Color','w');
imagesc(PM_list, nM_list, Rd_best_2D);
set(gca, 'YDir','normal'); colorbar;
xlabel('P_M'); ylabel('n_M');
title('Max R_d over P_K and d_K (in feasible set)');
end


%% ========================================================================
function [dK_vec, Rd_vec] = sweep_dK_for_fixed_resources(params, PM, PK, nM, nK)
% 用于画图：给定 (PM,PK,nM,nK) 时，扫描 d_K 计算 R_d 曲线（不再检查阈值）

sigma2 = params.sigma2;
zB = params.zB; zE = params.zE;

gBM = zB * PM / sigma2;
gBK = zB * PK / sigma2;
gEM = zE * PM / sigma2;
gEK = zE * PK / sigma2;

eBM = epsilon_fbl_awgn(nM, params.dM, gBM);
eEM = epsilon_fbl_awgn(nM, params.dM, gEM);

dK_vec = linspace(0, 2*nK, params.dK_grid_size);
Rd_vec = nan(size(dK_vec));

for idx = 1:numel(dK_vec)
    dK = dK_vec(idx);
    eBK = epsilon_fbl_awgn(nK, dK, gBK);
    eEK = epsilon_fbl_awgn(nK, dK, gEK);
    Rd_vec(idx) = Rd_metric(eBM, eBK, eEM, eEK);
end
end