clear; clc; close all;

%% ======================================================================
%  Physical Layer Deception in OFDM Systems - power-domain curvature check
%  Exact formulas follow the paper / Appendix A-B-C.
%  This script is written to avoid the numerical/path-selection issues that
%  arise if one directly sweeps P in [2,5] mW for all curves.
%
%  Main ideas of this version:
%  1) Exact curves (omega, epsilon and their derivatives) are plotted on a
%     LOG power axis over a range that covers the transition epsilon=0.5.
%  2) The lower bound \hat{\epsilon} is plotted locally around the anchor.
%  3) For \hat f^{(t)} and \hat g^{(t)}, only the RAW valid domain implied by
%     the paper formula is plotted, i.e., where the denominator is > 0.
%  4) A plot-clipped version max(\hat eps,0) is also shown ONLY for visual
%     convenience; the raw formula is kept unchanged and is used in checks.
%% ======================================================================

%% ====================== User parameters ================================
par.sigma2   = 1.0;      % [mW]
par.zB_dB    = 0.0;
par.zE_dB    = -10.0;
par.zB       = 10^(par.zB_dB/10);
par.zE       = 10^(par.zE_dB/10);

par.nM       = 64;
par.nK       = 64;
par.dM       = 16;
par.dK       = 16;

par.eps_th_BM = 0.5;
par.eps_th_EM = 0.5;
par.eps_th_BK = 0.5;
par.eps_th_EK = 0.5;
par.Tth       = 0.05;

tol = 1e-12;

%% -------------------- Characteristic powers ----------------------------
% omega_{i,j}(P)=0 <=> C(gamma)=d_j/n_j <=> P = sigma^2/z_i * (2^(d_j/n_j)-1)
P0_BM = power_at_eps_half(par.zB, par.nM, par.dM, par.sigma2);
P0_EM = power_at_eps_half(par.zE, par.nM, par.dM, par.sigma2);
P0_BK = power_at_eps_half(par.zB, par.nK, par.dK, par.sigma2);
P0_EK = power_at_eps_half(par.zE, par.nK, par.dK, par.sigma2);

% Exact-curve sweep: choose a broad range that covers both Bob and Eve
Pmin_exact = min([P0_BM, P0_BK, P0_EM, P0_EK]) / 20;
Pmax_exact = max([P0_BM, P0_BK, P0_EM, P0_EK]) * 3;
P_exact_M  = logspace(log10(Pmin_exact), log10(Pmax_exact), 2501);
P_exact_K  = P_exact_M;

%% -------------------- Exact omega / epsilon curves ---------------------
[w_BM, dw_BM, d2w_BM] = omega_with_derivatives(P_exact_M, par.zB, par.nM, par.dM, par.sigma2);
[w_EM, dw_EM, d2w_EM] = omega_with_derivatives(P_exact_M, par.zE, par.nM, par.dM, par.sigma2);
[w_BK, dw_BK, d2w_BK] = omega_with_derivatives(P_exact_K, par.zB, par.nK, par.dK, par.sigma2);
[w_EK, dw_EK, d2w_EK] = omega_with_derivatives(P_exact_K, par.zE, par.nK, par.dK, par.sigma2);

[eps_BM, deps_BM, d2eps_BM] = epsilon_with_derivatives(w_BM, dw_BM, d2w_BM);
[eps_EM, deps_EM, d2eps_EM] = epsilon_with_derivatives(w_EM, dw_EM, d2w_EM);
[eps_BK, deps_BK, d2eps_BK] = epsilon_with_derivatives(w_BK, dw_BK, d2w_BK);
[eps_EK, deps_EK, d2eps_EK] = epsilon_with_derivatives(w_EK, dw_EK, d2w_EK);

%% -------------------- Local study for lower bounds ---------------------
% Bob-M lower bound: choose a local window around its own transition point.
PM_hat_lb = 1.05 * P0_BM;
PM_lb     = linspace(max(Pmin_exact, 0.35*P0_BM), 2.0*P0_BM, 1401);
eps_BM_lb = epsilon_only(PM_lb, par.zB, par.nM, par.dM, par.sigma2);
epshat_BM_lb_raw  = eps_hat_lower_raw(PM_lb, par.zB, par.nM, par.dM, par.sigma2, PM_hat_lb);
epshat_BM_lb_clip = max(epshat_BM_lb_raw, 0);

% Eve-K lower bound: choose a local window around its transition point.
PK_hat_lb = 0.97 * P0_EK;
PK_lb     = linspace(0.90*P0_EK, 1.06*P0_EK, 1401);
eps_EK_lb = epsilon_only(PK_lb, par.zE, par.nK, par.dK, par.sigma2);
epshat_EK_lb_raw  = eps_hat_lower_raw(PK_lb, par.zE, par.nK, par.dK, par.sigma2, PK_hat_lb);
epshat_EK_lb_clip = max(epshat_EK_lb_raw, 0);

%% -------------------- Local study for \hat f^{(t)} and \hat g^{(t)} ------
% To make \hat f and \hat g visible, the fixed powers must be chosen where
% (1-eps_Eve,M) and raw epshat_Eve,K are both positive and not numerically 0.
PM_fix_fg = 1.08 * P0_EM;
PK_fix_fg = 0.97 * P0_EK;
PM_hat_fg = PM_fix_fg;
PK_hat_fg = PK_fix_fg;

PM_g = linspace(0.96*P0_EM, 1.20*P0_EM, 1501);
PK_f = linspace(0.90*P0_EK, 1.08*P0_EK, 1501);

epshat_BM_t = eps_hat_lower_raw(PM_fix_fg, par.zB, par.nM, par.dM, par.sigma2, PM_hat_fg);
eps_EM_t    = epsilon_only(PM_fix_fg,   par.zE, par.nM, par.dM, par.sigma2);
eps_BK_t    = epsilon_only(PK_fix_fg,   par.zB, par.nK, par.dK, par.sigma2);
epshat_EK_t = eps_hat_lower_raw(PK_fix_fg, par.zE, par.nK, par.dK, par.sigma2, PK_hat_fg);

eps_BK_curve   = epsilon_only(PK_f, par.zB, par.nK, par.dK, par.sigma2);
eps_EM_curve   = epsilon_only(PM_g, par.zE, par.nM, par.dM, par.sigma2);
epshat_BM_curve = eps_hat_lower_raw(PM_g, par.zB, par.nM, par.dM, par.sigma2, PM_hat_fg);
epshat_EK_curve = eps_hat_lower_raw(PK_f, par.zE, par.nK, par.dK, par.sigma2, PK_hat_fg);

A_f = 1 - (1 - epshat_BM_t) .* eps_BK_curve;
B_f = (1 - eps_EM_t) .* epshat_EK_curve;
A_g = (1 - eps_EM_curve) .* epshat_EK_t;
B_g = 1 - (1 - epshat_BM_curve) .* eps_BK_t;

% Eq. (8) / Eq. (11) are only meaningful on the raw-valid domain:
valid_f = (A_f >= 0) & (B_f > 0) & isfinite(A_f) & isfinite(B_f);
valid_g = (A_g >= 0) & (B_g > 0) & isfinite(A_g) & isfinite(B_g);

% Two y choices:
% 1) manual, convenient for observing shape;
% 2) algorithm line in the paper.
y_f_manual = 1.0;
y_g_manual = 1.0;

y_f_alg = sqrt( max(1 - (1 - epshat_BM_t) * eps_BK_t, 0) / max((1 - eps_EM_t) * epshat_EK_t, realmin) );
y_g_alg = sqrt( max((1 - eps_EM_t) * epshat_EK_t, 0) / max(1 - (1 - epshat_BM_t) * eps_BK_t, realmin) );

fh_manual = nan(size(PK_f));
fh_alg    = nan(size(PK_f));
gh_manual = nan(size(PM_g));
gh_alg    = nan(size(PM_g));

fh_manual(valid_f) = 2*y_f_manual*sqrt(A_f(valid_f)) - y_f_manual^2 ./ B_f(valid_f);
fh_alg(valid_f)    = 2*y_f_alg   *sqrt(A_f(valid_f)) - y_f_alg^2    ./ B_f(valid_f);
gh_manual(valid_g) = 2*y_g_manual*sqrt(A_g(valid_g)) - y_g_manual^2 ./ B_g(valid_g);
gh_alg(valid_g)    = 2*y_g_alg   *sqrt(A_g(valid_g)) - y_g_alg^2    ./ B_g(valid_g);

%% -------------------- Feasible masks on exact sweep --------------------
mask_PM = feasible_mask_PM(P_exact_M, PK_fix_fg, par);
mask_PK = feasible_mask_PK(P_exact_K, PM_fix_fg, par);

%% -------------------- Console summary ---------------------------------
fprintf('\n================ characteristic powers (epsilon = 0.5) ================\n');
fprintf('P0_BM = %.12g mW\n', P0_BM);
fprintf('P0_EM = %.12g mW\n', P0_EM);
fprintf('P0_BK = %.12g mW\n', P0_BK);
fprintf('P0_EK = %.12g mW\n', P0_EK);
fprintf('exact-curve sweep: [%.6g, %.6g] mW on log axis\n', Pmin_exact, Pmax_exact);

fprintf('\n================ exact sign / curvature summary =======================\n');
print_sign_summary('BM exact (vs P_M)', dw_BM, d2w_BM, deps_BM, d2eps_BM, mask_PM, tol);
print_sign_summary('EM exact (vs P_M)', dw_EM, d2w_EM, deps_EM, d2eps_EM, mask_PM, tol);
print_sign_summary('BK exact (vs P_K)', dw_BK, d2w_BK, deps_BK, d2eps_BK, mask_PK, tol);
print_sign_summary('EK exact (vs P_K)', dw_EK, d2w_EK, deps_EK, d2eps_EK, mask_PK, tol);

fprintf('\n================ lower-bound checks ===================================\n');
fprintf('Bob-M local tangency at PM_hat_lb? %d\n', abs(eps_hat_lower_raw(PM_hat_lb, par.zB, par.nM, par.dM, par.sigma2, PM_hat_lb) - epsilon_only(PM_hat_lb, par.zB, par.nM, par.dM, par.sigma2)) <= 1e-10);
fprintf('Eve-K local tangency at PK_hat_lb? %d\n', abs(eps_hat_lower_raw(PK_hat_lb, par.zE, par.nK, par.dK, par.sigma2, PK_hat_lb) - epsilon_only(PK_hat_lb, par.zE, par.nK, par.dK, par.sigma2)) <= 1e-10);
fprintf('Bob-M raw lower bound <= exact on local window? %d\n', all(epshat_BM_lb_raw <= eps_BM_lb + 1e-10));
fprintf('Eve-K raw lower bound <= exact on local window? %d\n', all(epshat_EK_lb_raw <= eps_EK_lb + 1e-10));
report_interval('Eve-K raw lower bound positive on PK_lb', PK_lb, epshat_EK_lb_raw > 0);

fprintf('\n================ surrogate-study setup =================================\n');
fprintf('PM_fix_fg = %.12g mW, PK_fix_fg = %.12g mW\n', PM_fix_fg, PK_fix_fg);
fprintf('PM_hat_fg = %.12g mW, PK_hat_fg = %.12g mW\n', PM_hat_fg, PK_hat_fg);
fprintf('epshat_BM_t = %.6e\n', epshat_BM_t);
fprintf('eps_EM_t    = %.6e\n', eps_EM_t);
fprintf('eps_BK_t    = %.6e\n', eps_BK_t);
fprintf('epshat_EK_t = %.6e\n', epshat_EK_t);
fprintf('y_f_manual = %.6g, y_f_alg = %.6g\n', y_f_manual, y_f_alg);
fprintf('y_g_manual = %.6g, y_g_alg = %.6g\n', y_g_manual, y_g_alg);
report_interval('valid raw domain of hat f^(t)', PK_f, valid_f);
report_interval('valid raw domain of hat g^(t)', PM_g, valid_g);

%% ====================== Figure 1: exact M-branch curves ===============
figure('Name', 'Exact curves vs P_M', 'Color', 'w');
subplot(3,2,1);
semilogx(P_exact_M, w_BM, 'LineWidth', 1.35); hold on;
semilogx(P_exact_M, w_EM, '--', 'LineWidth', 1.35);
xline(P0_BM, ':'); xline(P0_EM, ':');
grid on; xlabel('P_M [mW]'); ylabel('\omega_{i,M}(P_M)');
legend('\omega_{Bob,M}', '\omega_{Eve,M}', 'P0_{Bob,M}', 'P0_{Eve,M}', 'Location', 'best');
title('Exact \omega_{i,M}(P_M)');

subplot(3,2,2);
semilogx(P_exact_M, dw_BM, 'LineWidth', 1.35); hold on;
semilogx(P_exact_M, dw_EM, '--', 'LineWidth', 1.35);
grid on; xlabel('P_M [mW]'); ylabel('d\omega_{i,M}/dP_M');
legend('d\omega_{Bob,M}/dP_M', 'd\omega_{Eve,M}/dP_M', 'Location', 'best');
title('Exact first derivative');

subplot(3,2,3);
semilogx(P_exact_M, d2w_BM, 'LineWidth', 1.35); hold on;
semilogx(P_exact_M, d2w_EM, '--', 'LineWidth', 1.35);
grid on; xlabel('P_M [mW]'); ylabel('d^2\omega_{i,M}/dP_M^2');
legend('d^2\omega_{Bob,M}/dP_M^2', 'd^2\omega_{Eve,M}/dP_M^2', 'Location', 'best');
title('Exact second derivative');

subplot(3,2,4);
semilogx(P_exact_M, eps_BM, 'LineWidth', 1.35); hold on;
% semilogx(P_exact_M, eps_EM, '--', 'LineWidth', 1.35);
yline(0.5, ':'); xline(P0_BM, ':'); xline(P0_EM, ':');
grid on; xlabel('P_M [mW]'); ylabel('\epsilon_{i,M}(P_M)');
% legend('\epsilon_{Bob,M}', '\epsilon_{Eve,M}', '\epsilon=0.5', 'P0_{Bob,M}', 'P0_{Eve,M}', 'Location', 'best');
legend('\epsilon_{Bob,M}', '\epsilon=0.5', 'P0_{Bob,M}',  'Location', 'best');
title('Exact \epsilon_{i,M}(P_M)');

subplot(3,2,5);
semilogx(P_exact_M, deps_BM, 'LineWidth', 1.35); hold on;
% semilogx(P_exact_M, deps_EM, '--', 'LineWidth', 1.35);
grid on; xlabel('P_M [mW]'); ylabel('d\epsilon_{i,M}/dP_M');
legend('d\epsilon_{Bob,M}/dP_M', 'Location', 'best');
% legend('d\epsilon_{Bob,M}/dP_M', 'd\epsilon_{Eve,M}/dP_M', 'Location', 'best');
title('Exact first derivative');

subplot(3,2,6);
semilogx(P_exact_M, d2eps_BM, 'LineWidth', 1.35); hold on;
% semilogx(P_exact_M, d2eps_EM, '--', 'LineWidth', 1.35);
grid on; xlabel('P_M [mW]'); ylabel('d^2\epsilon_{i,M}/dP_M^2');
% legend('d^2\epsilon_{Bob,M}/dP_M^2', 'd^2\epsilon_{Eve,M}/dP_M^2', 'Location', 'best');
legend('d^2\epsilon_{Bob,M}/dP_M^2', 'Location', 'best');
title('Exact second derivative');
sgtitle('Exact power-domain curves for M branch');

%% ====================== Figure 2: exact K-branch curves ===============
figure('Name', 'Exact curves vs P_K', 'Color', 'w');
subplot(3,2,1);
semilogx(P_exact_K, w_BK, 'LineWidth', 1.35); hold on;
semilogx(P_exact_K, w_EK, '--', 'LineWidth', 1.35);
xline(P0_BK, ':'); xline(P0_EK, ':');
grid on; xlabel('P_K [mW]'); ylabel('\omega_{i,K}(P_K)');
legend('\omega_{Bob,K}', '\omega_{Eve,K}', 'P0_{Bob,K}', 'P0_{Eve,K}', 'Location', 'best');
title('Exact \omega_{i,K}(P_K)');

subplot(3,2,2);
semilogx(P_exact_K, dw_BK, 'LineWidth', 1.35); hold on;
semilogx(P_exact_K, dw_EK, '--', 'LineWidth', 1.35);
grid on; xlabel('P_K [mW]'); ylabel('d\omega_{i,K}/dP_K');
legend('d\omega_{Bob,K}/dP_K', 'd\omega_{Eve,K}/dP_K', 'Location', 'best');
title('Exact first derivative');

subplot(3,2,3);
semilogx(P_exact_K, d2w_BK, 'LineWidth', 1.35); hold on;
semilogx(P_exact_K, d2w_EK, '--', 'LineWidth', 1.35);
grid on; xlabel('P_K [mW]'); ylabel('d^2\omega_{i,K}/dP_K^2');
legend('d^2\omega_{Bob,K}/dP_K^2', 'd^2\omega_{Eve,K}/dP_K^2', 'Location', 'best');
title('Exact second derivative');

subplot(3,2,4);
semilogx(P_exact_K, eps_BK, 'LineWidth', 1.35); hold on;
% semilogx(P_exact_K, eps_EK, '--', 'LineWidth', 1.35);
yline(0.5, ':'); xline(P0_BK, ':'); xline(P0_EK, ':');
grid on; xlabel('P_K [mW]'); ylabel('\epsilon_{i,K}(P_K)');
% legend('\epsilon_{Bob,K}', '\epsilon_{Eve,K}', '\epsilon=0.5', 'P0_{Bob,K}', 'P0_{Eve,K}', 'Location', 'best');
legend('\epsilon_{Bob,K}', '\epsilon=0.5', 'P0_{Bob,K}', 'Location', 'best');
title('Exact \epsilon_{i,K}(P_K)');

subplot(3,2,5);
semilogx(P_exact_K, deps_BK, 'LineWidth', 1.35); hold on;
% semilogx(P_exact_K, deps_EK, '--', 'LineWidth', 1.35);
grid on; xlabel('P_K [mW]'); ylabel('d\epsilon_{i,K}/dP_K');
% legend('d\epsilon_{Bob,K}/dP_K', 'd\epsilon_{Eve,K}/dP_K', 'Location', 'best');
legend('d\epsilon_{Bob,K}/dP_K', 'Location', 'best');
title('Exact first derivative');

subplot(3,2,6);
semilogx(P_exact_K, d2eps_BK, 'LineWidth', 1.35); hold on;
% semilogx(P_exact_K, d2eps_EK, '--', 'LineWidth', 1.35);
grid on; xlabel('P_K [mW]'); ylabel('d^2\epsilon_{i,K}/dP_K^2');
% legend('d^2\epsilon_{Bob,K}/dP_K^2', 'd^2\epsilon_{Eve,K}/dP_K^2', 'Location', 'best');
legend('d^2\epsilon_{Bob,K}/dP_K^2', 'Location', 'best');
title('Exact second derivative');
sgtitle('Exact power-domain curves for K branch');

%% ====================== Figure 3: local lower bounds ===================
% figure('Name', 'Local lower-bound study', 'Color', 'w');
% subplot(1,2,1);
% plot(PM_lb, eps_BM_lb, 'LineWidth', 1.35); hold on;
% plot(PM_lb, epshat_BM_lb_raw, '--', 'LineWidth', 1.35);
% % plot(PM_lb, epshat_BM_lb_clip, ':', 'LineWidth', 1.35);
% % xline(PM_hat_lb, ':');
% grid on; xlabel('P_M [mW]'); ylabel('Bob-M local lower bound');
% legend('\epsilon_{Bob,M}', '\hat{\epsilon}_{Bob,M} raw', 'max(\hat{\epsilon}_{Bob,M},0) for plot', 'P_M hat', 'Location', 'best');
% title('Local Bob-M lower bound');

% subplot(1,2,2);
subplot(1,1,1);
plot(PK_lb, eps_EK_lb, 'LineWidth', 1.35); hold on;
plot(PK_lb, epshat_EK_lb_raw, '--', 'LineWidth', 1.35);
% plot(PK_lb, epshat_EK_lb_clip, ':', 'LineWidth', 1.35);
% xline(PK_hat_lb, ':');
grid on; xlabel('P_K [mW]'); ylabel('Eve-K local lower bound');
legend('\epsilon_{Eve,K}', '\hat{\epsilon}_{Eve,K} raw', 'max(\hat{\epsilon}_{Eve,K},0) for plot', 'P_K hat', 'Location', 'best');
title('Local Eve-K lower bound');
sgtitle('Local study of the Appendix-A lower bound');

%% ====================== Figure 4: \hat f and \hat g ====================
figure('Name', 'hat f and hat g on raw-valid domain', 'Color', 'w');
subplot(1,2,1);
plot(PK_f(valid_f), fh_manual(valid_f), 'LineWidth', 1.45); hold on;
plot(PK_f(valid_f), fh_alg(valid_f), '--', 'LineWidth', 1.45);
xline(PK_fix_fg, ':');
grid on; xlabel('P_K [mW]'); ylabel('\hat f^{(t)}');
legend(sprintf('\hat f^{(t)} manual y=%.3g', y_f_manual), sprintf('\hat f^{(t)} alg y=%.3g', y_f_alg), 'P_K^{(t)}', 'Location', 'best');
title('Eq. (8) on raw-valid domain');

subplot(1,2,2);
plot(PM_g(valid_g), gh_manual(valid_g), 'LineWidth', 1.45); hold on;
plot(PM_g(valid_g), gh_alg(valid_g), '--', 'LineWidth', 1.45);
xline(PM_fix_fg, ':');
grid on; xlabel('P_M [mW]'); ylabel('\hat g^{(t)}');
legend(sprintf('\hat g^{(t)} manual y=%.3g', y_g_manual), sprintf('\hat g^{(t)} alg y=%.3g', y_g_alg), 'P_M^{(t)}', 'Location', 'best');
title('Eq. (11) on raw-valid domain');
sgtitle('Surrogate functions from the paper: raw-valid plotting only');

%% ====================== local functions =================================
function P0 = power_at_eps_half(z, n, d, sigma2)
    P0 = sigma2 / z * (2^(d/n) - 1);
end

function [omega, domega_dP, d2omega_dP2] = omega_with_derivatives(P, z, n, d, sigma2)
    alpha = z / sigma2;
    gamma = alpha .* P;
    r = d / n;

    V  = 1 - (1 + gamma).^(-2);
    A  = log(1 + gamma) - r * log(2);  % = [C(gamma) - d/n] ln 2
    B  = V.^(-1/2);

    V1 = 2 * (1 + gamma).^(-3);
    V2 = -6 * (1 + gamma).^(-4);
    A1 = (1 + gamma).^(-1);
    A2 = -(1 + gamma).^(-2);
    B1 = -0.5 * V.^(-3/2) .* V1;
    B2 = 0.75 * V.^(-5/2) .* (V1.^2) - 0.5 * V.^(-3/2) .* V2;

    sqn = sqrt(n);
    omega       = sqn .* A .* B;
    domega_dg   = sqn .* (A1 .* B + A .* B1);
    d2omega_dg2 = sqn .* (A2 .* B + 2 .* A1 .* B1 + A .* B2);

    domega_dP   = alpha .* domega_dg;
    d2omega_dP2 = (alpha.^2) .* d2omega_dg2;
end

function [epsv, deps_dP, d2eps_dP2] = epsilon_with_derivatives(omega, domega_dP, d2omega_dP2)
    phi = exp(-0.5 .* omega.^2 - 0.5*log(2*pi));
    epsv = 0.5 .* erfc(omega ./ sqrt(2));
    deps_dP = -phi .* domega_dP;
    d2eps_dP2 = phi .* (omega .* (domega_dP.^2) - d2omega_dP2);
end

function epsv = epsilon_only(P, z, n, d, sigma2)
    [omega, ~, ~] = omega_with_derivatives(P, z, n, d, sigma2);
    epsv = 0.5 .* erfc(omega ./ sqrt(2));
end

% function epshat = eps_hat_lower_raw(P, z, n, d, sigma2, P_hat)
%     [omega_hat, ~, ~] = omega_with_derivatives(P_hat, z, n, d, sigma2);
%     [a_m, b_m, c_m]   = abc_qbound(-omega_hat);
%     [omega, ~, ~]     = omega_with_derivatives(P, z, n, d, sigma2);
%     epshat = 1 - b_m .* exp(-a_m .* omega) - c_m;
% end

function epshat = eps_hat_lower_raw(P, z, n, d, sigma2, P_hat)
    [omega_hat, ~, ~] = omega_with_derivatives(P_hat, z, n, d, sigma2);
    [omega, ~, ~]     = omega_with_derivatives(P,     z, n, d, sigma2);

    if omega_hat >= 0
        % Q(omega) with omega >= 0
        [a_m, b_m, c_m] = abc_qbound(omega_hat);
        epshat = b_m .* exp(-a_m .* omega) + c_m;
    else
        % Q(omega) = 1 - Q(-omega), with -omega >= 0
        [a_m, b_m, c_m] = abc_qbound(-omega_hat);
        epshat = 1 - b_m .* exp(a_m .* omega) - c_m;
    end
end

function [a, b, c] = abc_qbound(w_hat)
    Qhat = 0.5 .* erfc(w_hat ./ sqrt(2));
    a = max( exp(-0.5 .* w_hat.^2) ./ (sqrt(2*pi) .* Qhat), w_hat );
    b = exp(a .* w_hat - 0.5 .* w_hat.^2) ./ (sqrt(2*pi) .* a);
    c = Qhat - b .* exp(-a .* w_hat);
end

% function epsLF = eps_lf(PM, PK, par)
%     epsBM = epsilon_only(PM, par.zB, par.nM, par.dM, par.sigma2);
%     epsEM = epsilon_only(PM, par.zE, par.nM, par.dM, par.sigma2);
%     epsBK = epsilon_only(PK, par.zB, par.nK, par.dK, par.sigma2);
%     epsEK = epsilon_only(PK, par.zE, par.nK, par.dK, par.sigma2);
%     epsLF = 1 - (1 - epsBM) .* (1 - epsBK) .* (1 - (1 - epsEM) .* (1 - epsEK));
% end

% function T = throughput(PM, PK, par)
%     epsLF = eps_lf(PM, PK, par);
%     T = (1 - epsLF) .* (par.dM / (par.nM + par.nK));
% end

function mask = feasible_mask_PM(PM_grid, PK_fix, par)
    epsBM = epsilon_only(PM_grid, par.zB, par.nM, par.dM, par.sigma2);
    epsEM = epsilon_only(PM_grid, par.zE, par.nM, par.dM, par.sigma2);
    epsBK = epsilon_only(PK_fix,  par.zB, par.nK, par.dK, par.sigma2);
    epsEK = epsilon_only(PK_fix,  par.zE, par.nK, par.dK, par.sigma2);
    % T     = throughput(PM_grid, PK_fix, par);
    mask = (epsBM <= par.eps_th_BM) & ...
           (epsEM <= par.eps_th_EM) & ...
           (epsBK <= par.eps_th_BK) & ...
           (epsEK >= par.eps_th_EK) ;
        %    (T >= par.Tth);
end

function mask = feasible_mask_PK(PK_grid, PM_fix, par)
    epsBK = epsilon_only(PK_grid, par.zB, par.nK, par.dK, par.sigma2);
    epsEK = epsilon_only(PK_grid, par.zE, par.nK, par.dK, par.sigma2);
    epsBM = epsilon_only(PM_fix,  par.zB, par.nM, par.dM, par.sigma2);
    epsEM = epsilon_only(PM_fix,  par.zE, par.nM, par.dM, par.sigma2);
    % T     = throughput(PM_fix, PK_grid, par);
    mask = (epsBM <= par.eps_th_BM) & ...
           (epsEM <= par.eps_th_EM) & ...
           (epsBK <= par.eps_th_BK) & ...
           (epsEK >= par.eps_th_EK) ;
        %    (T >= par.Tth);
end

function print_sign_summary(name, dw, d2w, deps, d2eps, mask, tol)
    fprintf('%s:\n', name);
    fprintf('  d omega / dP  >= 0 on full sweep?      %d\n', all(dw >= -tol));
    fprintf('  d2 omega / dP2 <= 0 on full sweep?     %d\n', all(d2w <= tol));
    fprintf('  d eps / dP    <= 0 on full sweep?      %d\n', all(deps <= tol));
    fprintf('  d2 eps / dP2 >= 0 on full sweep?       %d\n', all(d2eps >= -tol));
    if any(mask)
        fprintf('  d omega / dP  >= 0 on feasible sweep?  %d\n', all(dw(mask) >= -tol));
        fprintf('  d2 omega / dP2 <= 0 on feasible sweep? %d\n', all(d2w(mask) <= tol));
        fprintf('  d eps / dP    <= 0 on feasible sweep?  %d\n', all(deps(mask) <= tol));
        fprintf('  d2 eps / dP2 >= 0 on feasible sweep?   %d\n', all(d2eps(mask) >= -tol));
    else
        fprintf('  feasible sweep is empty under current setup.\n');
    end
end

function report_interval(name, x, mask)
    if any(mask)
        idx = find(mask);
        fprintf('%s: [%.12g, %.12g]\n', name, x(idx(1)), x(idx(end)));
    else
        fprintf('%s: empty\n', name);
    end
end
