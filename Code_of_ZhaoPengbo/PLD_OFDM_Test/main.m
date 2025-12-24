function main()
% main.m 
% PLD joint optimization (OFDM, orthogonal):
% Variables: P_M, P_K, n_M, n_K, d_K (integer)
% Goal: maximize R_d
%
%
%
% changes :
% 1) d_K is optimized over integers by enumeration in feasible bounds -> rigorous integer output
% 2) Power block solved by coarse-to-fine grid search (global-ish, avoids gradient flat regions)
% 3) n block solved by integer scan (default enforce nM+nK=Ntot)
% 4) Accept only feasible candidates; keep best with tie-break (optional)

close all; 
% clc;
fprintf('==== Rigorous PLD Joint Optimization (OFDM, integer dK) ====\n');

%% ---------------- Parameters ----------------
params.sigma2 = 1;

hB = 1.2; params.zB = hB^2;
hE = 0.7; params.zE = hE^2;

params.dM = 32;

params.Ntot = 200;
params.Etot = 600;
params.Pmax = 5;

params.epsBM_th  = 0.5;
params.epsBK_th  = 0.5;
params.epsEM_th  = 0.5;
params.epsEK_min = 0.5;

% Outer iteration controls
params.max_outer_iter  = 60;
params.min_outer_iter  = 12;      
params.stop_patience   = 8;       % 连续多轮无显著提升才停
params.outer_abs_tol   = 1e-13;   % 绝对提升阈值（更严格）
params.outer_rel_tol   = 1e-11;   % 相对提升阈值

% Progressive refinement schedules
params.n_window        = 80;
params.n_step_schedule = [8 4 2 1];

params.p_levels_per_outer = 1;           % 每轮只做一层功率网格
params.p_grid_schedule    = [11 15 21 31 41];
params.p_win0             = params.Pmax; % 初始功率窗口
params.p_win_decay        = 0.85;        % 每轮缩小一点（慢）
params.p_win_min          = 0.05*params.Pmax;


% n-search
params.enforce_full_blocklength = true;  % if true: nM+nK=Ntot (recommended)
params.n_min_int = 1;

% power search (coarse-to-fine)
params.p_levels = 3;           % refinement levels
params.p_grid0  = 21;          % grid per dimension at level 1
params.p_shrink = 0.25;        % local window shrink per level

% tie-break among near-equal Rd solutions
params.tie_tol = 1e-12;        % treat Rd equal within this tolerance
params.tie_mode = "min_sumP";  % "min_sumP" or "min_energy" or "none"

%% ---------------- Run BCD search ----------------
res = pld_joint_opt_grid(params);

%% ---------------- Print results ----------------
if ~res.feasible
    fprintf('\nNo feasible solution found.\n');
    return;
end

fprintf('\n==== Final Solution ====\n');
fprintf('n_M = %d, n_K = %d (sum=%d)\n', res.nM, res.nK, res.nM + res.nK);
fprintf('P_M = %.6f, P_K = %.6f\n', res.PM, res.PK);
fprintf('Energy = %.6f (<= %.6f)\n', res.PM*res.nM + res.PK*res.nK, params.Etot);
fprintf('d_K* = %d (range [0, 2n_K]=[0, %d])\n', res.dK, 2*res.nK);
fprintf('R_d* = %.12f\n', res.Rd);

figure('Color','w');
plot(res.Rd_hist, 'o-', 'LineWidth', 2); grid on;
xlabel('Outer iteration');
ylabel('R_d');
title('R_d over iterations');

end


%% ========================================================================
function res = pld_joint_opt_grid(params)

% ---------- Feasible initialization ----------
[PM, PK, nM, nK, dK, Rd, ok] = find_feasible_init(params);
if ~ok
    res.feasible = false; return;
end

best = pack_sol(PM, PK, nM, nK, dK, Rd);
Rd_hist = best.Rd;

fprintf('\n[Init feasible]\n');
fprintf('  n=[%d,%d], P=[%.3f,%.3f], dK=%d, Rd=%.6f\n', best.nM,best.nK,best.PM,best.PK,best.dK,best.Rd);

% progressive states
p_win = params.p_win0;
stall = 0;

for it = 1:params.max_outer_iter
    Rd_prev = best.Rd;

    % -------- choose refinement for this outer iter --------
    n_step = params.n_step_schedule(min(it, numel(params.n_step_schedule)));
    p_grid = params.p_grid_schedule(min(it, numel(params.p_grid_schedule)));

    % -------- (A) n-update: local + coarse-to-fine --------
    [cand_n, okn] = opt_n_scan_local(best.PM, best.PK, best.nM, best.nK, n_step, params);
    if okn
        best = accept_if_better(best, cand_n, params);
    end

    % -------- (B) P-update: one-level grid with progressive p_grid --------
    [cand_p, okp] = opt_p_grid_onelevel(best.nM, best.nK, best.PM, best.PK, p_win, p_grid, params);
    if okp
        best = accept_if_better(best, cand_p, params);
    end

    % update window slowly (so improvements can appear across many iterations)
    p_win = max(params.p_win_min, p_win * params.p_win_decay);

    Rd_hist(end+1) = best.Rd; %#ok<AGROW>

    dRd = best.Rd - Rd_prev;
    fprintf('[Iter %02d] Rd=%.12f | Δ=%.3e | n=[%d,%d] P=[%.4f,%.4f] dK=%d | (nStep=%d, pGrid=%d, pWin=%.3f)\n', ...
        it, best.Rd, dRd, best.nM, best.nK, best.PM, best.PK, best.dK, n_step, p_grid, p_win);

    % -------- stopping: min iter + patience + abs/rel tol --------
    tol_eff = max(params.outer_abs_tol, params.outer_rel_tol * max(1, abs(Rd_prev)));
    if dRd <= tol_eff
        stall = stall + 1;
    else
        stall = 0;
    end

    if it >= params.min_outer_iter && stall >= params.stop_patience
        fprintf('Stop: stall=%d >= %d after minIter=%d (tolEff=%.3e)\n', ...
            stall, params.stop_patience, params.min_outer_iter, tol_eff);
        break;
    end
end

res = best;
res.feasible = true;
res.Rd_hist = Rd_hist;

end

function best = accept_if_better(best, cand, params)
if cand.Rd > best.Rd + params.tie_tol
    best = cand;
elseif abs(cand.Rd - best.Rd) <= params.tie_tol
    best = tie_break(best, cand, params);
end
end


%% ========================================================================
function [cand, ok] = opt_n_scan_local(PM, PK, nM0, nK0, step, params)
ok = false; cand = struct();
best_Rd = -inf; best_sol = [];

N = params.Ntot;
w = params.n_window;

% local candidates around current nM*
nM_min = max(params.n_min_int, nM0 - w);
nM_max = min(N - params.n_min_int, nM0 + w);

for nM = nM_min : step : nM_max
    nK = N - nM; % enforce full length
    [dK, Rd, feas] = eval_best_dK_integer(PM, PK, nM, nK, params);
    if ~feas, continue; end

    s = pack_sol(PM, PK, nM, nK, dK, Rd);
    if Rd > best_Rd + params.tie_tol
        best_Rd = Rd; best_sol = s; ok = true;
    elseif abs(Rd - best_Rd) <= params.tie_tol && ok
        best_sol = tie_break(best_sol, s, params);
    end
end

if ok, cand = best_sol; end
end


%% ========================================================================
function [cand, ok] = opt_n_scan(PM, PK, nM0, nK0, params)
% Optimize integer (nM,nK) with P fixed, by scan.
% Default enforce nM+nK=Ntot, since objective prefers longer blocklength.

ok = false;
cand = struct();

best_local = [];
best_Rd = -inf;

if params.enforce_full_blocklength
    N = params.Ntot;
    for nM = params.n_min_int : (N - params.n_min_int)
        nK = N - nM;
        [dK, Rd, feas] = eval_best_dK_integer(PM, PK, nM, nK, params);
        if ~feas, continue; end

        s = pack_sol(PM, PK, nM, nK, dK, Rd);
        if Rd > best_Rd + params.tie_tol
            best_Rd = Rd; best_local = s; ok = true;
        elseif abs(Rd - best_Rd) <= params.tie_tol && ok
            best_local = tie_break(best_local, s, params);
        end
    end
else
    % allow nM+nK<=Ntot (slower but still okay for N=200)
    for nM = params.n_min_int : params.Ntot
        for nK = params.n_min_int : (params.Ntot - nM)
            [dK, Rd, feas] = eval_best_dK_integer(PM, PK, nM, nK, params);
            if ~feas, continue; end
            s = pack_sol(PM, PK, nM, nK, dK, Rd);
            if Rd > best_Rd + params.tie_tol
                best_Rd = Rd; best_local = s; ok = true;
            elseif abs(Rd - best_Rd) <= params.tie_tol && ok
                best_local = tie_break(best_local, s, params);
            end
        end
    end
end

if ok
    cand = best_local;
end

end



%% ========================================================================
function [cand, ok] = opt_p_grid_onelevel(nM, nK, PM0, PK0, win, gridN, params)
ok = false; cand = struct();
best_Rd = -inf; best_sol = [];

PM_min = max(0, PM0 - win);
PM_max = min(params.Pmax, PM0 + win);
PK_min = max(0, PK0 - win);
PK_max = min(params.Pmax, PK0 + win);

PM_list = linspace(PM_min, PM_max, gridN);
PK_list = linspace(PK_min, PK_max, gridN);

for i = 1:numel(PM_list)
    PM = PM_list(i);
    for j = 1:numel(PK_list)
        PK = PK_list(j);

        if ~check_resource(PM, PK, nM, nK, params), continue; end

        [dK, Rd, feas] = eval_best_dK_integer(PM, PK, nM, nK, params);
        if ~feas, continue; end

        s = pack_sol(PM, PK, nM, nK, dK, Rd);

        if Rd > best_Rd + params.tie_tol
            best_Rd = Rd; best_sol = s; ok = true;
        elseif abs(Rd - best_Rd) <= params.tie_tol && ok
            best_sol = tie_break(best_sol, s, params);
        end
    end
end

if ok, cand = best_sol; end
end

%% ========================================================================
function [cand, ok] = opt_p_grid(nM, nK, PM0, PK0, params)
% Optimize (PM,PK) with n fixed using coarse-to-fine grid search.
% At each grid point, optimize integer dK exactly by enumeration.

ok = false;
cand = struct();
best_Rd = -inf;
best_sol = [];

% initial center
cPM = PM0; cPK = PK0;
win = params.Pmax;  % first level: full [0,Pmax]

for lv = 1:params.p_levels
    if lv == 1
        PM_min = 0; PM_max = params.Pmax;
        PK_min = 0; PK_max = params.Pmax;
        Ng = params.p_grid0;
    else
        % local window around current best
        PM_min = max(0, cPM - win);
        PM_max = min(params.Pmax, cPM + win);
        PK_min = max(0, cPK - win);
        PK_max = min(params.Pmax, cPK + win);
        Ng = max(11, round(params.p_grid0/2));
    end

    PM_list = linspace(PM_min, PM_max, Ng);
    PK_list = linspace(PK_min, PK_max, Ng);

    for i = 1:numel(PM_list)
        PM = PM_list(i);
        for j = 1:numel(PK_list)
            PK = PK_list(j);

            % resource constraints first
            if ~check_resource(PM, PK, nM, nK, params), continue; end

            [dK, Rd, feas] = eval_best_dK_integer(PM, PK, nM, nK, params);
            if ~feas, continue; end

            s = pack_sol(PM, PK, nM, nK, dK, Rd);

            if Rd > best_Rd + params.tie_tol
                best_Rd = Rd; best_sol = s; ok = true;
            elseif abs(Rd - best_Rd) <= params.tie_tol && ok
                best_sol = tie_break(best_sol, s, params);
            end
        end
    end

    if ok
        cPM = best_sol.PM; cPK = best_sol.PK;
        win = win * params.p_shrink;
    else
        break;
    end
end

if ok
    cand = best_sol;
end

end


%% ========================================================================
function [dK_star, Rd_star, feasible] = eval_best_dK_integer(PM, PK, nM, nK, params)
% Given (PM,PK,nM,nK), compute best INTEGER dK in feasible interval:
% Constraints:
%   epsBM<=th, epsEM<=th must hold (independent of dK)
%   epsBK(dK)<=th, epsEK(dK)>=min
%   dK in {0,1,...,2*nK}

feasible = false;
dK_star = NaN;
Rd_star = NaN;

if ~check_resource(PM, PK, nM, nK, params)
    return;
end

sigma2 = params.sigma2;
zB = params.zB; zE = params.zE;

gBM = zB*PM/sigma2;
gBK = zB*PK/sigma2;
gEM = zE*PM/sigma2;
gEK = zE*PK/sigma2;

% M-block eps (fixed wrt dK)
eBM = epsilon_fbl_awgn(nM, params.dM, gBM);
eEM = epsilon_fbl_awgn(nM, params.dM, gEM);

if eBM > params.epsBM_th || eEM > params.epsEM_th
    return;
end

% Exact feasible integer bounds for dK
[dL, dU, feas] = dK_bounds_integer(nK, gBK, gEK, params);
if ~feas
    return;
end

best_Rd = -inf;
best_dK = dL;

for dK = dL:dU
    eBK = epsilon_fbl_awgn(nK, dK, gBK);
    eEK = epsilon_fbl_awgn(nK, dK, gEK);

    if eBK <= params.epsBK_th && eEK >= params.epsEK_min
        Rd = Rd_metric(eBM, eBK, eEM, eEK);
        if Rd > best_Rd
            best_Rd = Rd;
            best_dK = dK;
            feasible = true;
        end
    end
end

if feasible
    dK_star = best_dK;
    Rd_star = best_Rd;
end

end


%% ========================================================================
function [dL, dU, feasible] = dK_bounds_integer(nK, gBK, gEK, params)
% Derive integer bounds from monotonicity using X(d)=A-kd (affine in d):
% e=Q(X)
% Bob: eBK<=th -> XBK >= Q^{-1}(th) -> d <= (A_B - x_th)/k_B
% Eve: eEK>=min -> XEK <= Q^{-1}(min) -> d >= (A_E - x_min)/k_E
% Hard: d in [0, 2nK]

feasible = true;

[C_B, V_B] = CV_awgn(gBK);
[C_E, V_E] = CV_awgn(gEK);

A_B = log(2) * sqrt(nK / V_B) * C_B;
k_B = log(2) * sqrt(1 / (nK * V_B));

A_E = log(2) * sqrt(nK / V_E) * C_E;
k_E = log(2) * sqrt(1 / (nK * V_E));

x_th  = qinv_local(params.epsBK_th);
x_min = qinv_local(params.epsEK_min);

dU_real = (A_B - x_th) / k_B;
dL_real = (A_E - x_min) / k_E;

dL = ceil(max(0, dL_real));
dU = floor(min(2*nK, dU_real));

if ~(isfinite(dL) && isfinite(dU)) || dL > dU
    feasible = false;
end

end


%% ========================================================================
function ok = check_resource(PM, PK, nM, nK, params)
ok = true;
if PM < 0 || PK < 0 || PM > params.Pmax || PK > params.Pmax
    ok = false; return;
end
if nM < params.n_min_int || nK < params.n_min_int
    ok = false; return;
end
if nM + nK > params.Ntot
    ok = false; return;
end
if PM*nM + PK*nK > params.Etot + 1e-12
    ok = false; return;
end
end


%% ========================================================================
function rd = Rd_metric(eBM, eBK, eEM, eEK)
rd = (1 - (1 - eBM)*eBK) * (1 - eEM) * eEK;
end


%% ========================================================================
function epsv = epsilon_fbl_awgn(n, d, gamma)
% eps = Q(X), X = ln2*sqrt(n/V)*(C - d/n)
if gamma <= 0 || n <= 0 || d < 0
    epsv = 1; return;
end
[C, V] = CV_awgn(gamma);
X = log(2) * sqrt(n / V) * (C - d/n);
epsv = qfunc_local(X);
epsv = min(max(epsv, 0), 1);
end


%% ========================================================================
function [C, V] = CV_awgn(gamma)
C = log2(1 + gamma);
V = 1 - 1/(1 + gamma)^2;
end


%% ========================================================================
function y = qfunc_local(x)
y = 0.5 * erfc(x / sqrt(2));
end

function x = qinv_local(p)
% Q(x)=p -> x = sqrt(2)*erfcinv(2p)
p = min(max(p, 1e-15), 1-1e-15);
x = sqrt(2) * erfcinv(2*p);
end


%% ========================================================================
function [PM, PK, nM, nK, dK, Rd, ok] = find_feasible_init(params)
% Deterministic feasible initialization:
% - scan a few n splits
% - choose moderate powers within energy
% - pick best feasible by dK enumeration

ok = false;
PM = NaN; PK = NaN; nM = NaN; nK = NaN; dK = NaN; Rd = NaN;

% splits = linspace(0.3, 0.7, 9);
splits = linspace(0.1, 0.9, 15);
best_Rd = -inf;

for s = splits
    if params.enforce_full_blocklength
        N = params.Ntot;
        nM0 = max(params.n_min_int, round(s*N));
        nM0 = min(nM0, N - params.n_min_int);
        nK0 = N - nM0;
    else
        nM0 = max(params.n_min_int, round(s*params.Ntot));
        nK0 = max(params.n_min_int, params.Ntot - nM0);
    end

    % pick a few power candidates (low/medium/high) but within energy
    Pcand = [0.2, 0.5, 0.8] * params.Pmax;
    for a = 1:numel(Pcand)
        for b = 1:numel(Pcand)
            PM0 = Pcand(a);
            PK0 = Pcand(b);
            if ~check_resource(PM0, PK0, nM0, nK0, params), continue; end

            [dK0, Rd0, feas] = eval_best_dK_integer(PM0, PK0, nM0, nK0, params);
            if ~feas, continue; end

            if Rd0 > best_Rd + params.tie_tol
                best_Rd = Rd0;
                PM = PM0; PK = PK0; nM = nM0; nK = nK0; dK = dK0; Rd = Rd0;
                ok = true;
            elseif abs(Rd0 - best_Rd) <= params.tie_tol && ok
                cur = pack_sol(PM,PK,nM,nK,dK,Rd);
                snew = pack_sol(PM0,PK0,nM0,nK0,dK0,Rd0);
                best = tie_break(cur, snew, params);
                PM=best.PM; PK=best.PK; nM=best.nM; nK=best.nK; dK=best.dK; Rd=best.Rd;
            end
        end
    end
end

end


%% ========================================================================
function s = pack_sol(PM, PK, nM, nK, dK, Rd)
s = struct('PM',PM,'PK',PK,'nM',nM,'nK',nK,'dK',dK,'Rd',Rd);
end


%% ========================================================================
function best = tie_break(a, b, params)
% Choose between two solutions with near-equal Rd
best = a;
if abs(a.Rd - b.Rd) > params.tie_tol
    if b.Rd > a.Rd, best = b; end
    return;
end

switch params.tie_mode
    case "none"
        % keep a
    case "min_sumP"
        if (b.PM + b.PK) < (a.PM + a.PK)
            best = b;
        end
    case "min_energy"
        Ea = a.PM*a.nM + a.PK*a.nK;
        Eb = b.PM*b.nM + b.PK*b.nK;
        if Eb < Ea
            best = b;
        end
    otherwise
        % default
        if (b.PM + b.PK) < (a.PM + a.PK)
            best = b;
        end
end

end