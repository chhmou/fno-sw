% RESIDUAL_FNO  Governing-equation residuals for the FNO results.
% Set EXPERIMENT = 1 (free-fall limit) or 2 (linear limit). Reads the .mat
% files exported by the notebooks into ../results.
% Uses check_wave_residual.m (same folder or on the path).
%
% For every test solution in a results file, prints the relative residual of
% the truth fields (the floor of the check at 101 timeslices, ~1.5e-5) and
% of the FNO prediction. For the exploratory eps_try predictions, where no
% period is known, the period is extrapolated quadratically in eps over the
% same solutions used to construct the input profiles.

EXPERIMENT = 1;                       % 1 or 2

DATA_DIR = getenv('STANDING_WAVE_DATA');
if isempty(DATA_DIR), DATA_DIR = fullfile('..', 'data'); end
MAT_DIR = fullfile('..', 'results');
if EXPERIMENT == 1
    MAT_FILE = fullfile(MAT_DIR, 'fno_wave_results.mat');
    TRY_FILE = fullfile(MAT_DIR, 'fno_wave_results_eps_try.mat');                 % optional
else
    MAT_FILE = fullfile(MAT_DIR, 'fno_wave_results_ood_smallnu.mat');
    TRY_FILE = fullfile(MAT_DIR, 'fno_wave_results_ood_smallnu_eps_try.mat');     % optional
end

P = load(fullfile(DATA_DIR, 'parameters'), '-ascii');

%% ---------------- main results: truth floor vs FNO prediction -------------
R = load(MAT_FILE);
idx = double(R.test_indices(:)); ood = logical(R.is_ood(:));
Yt = double(R.Y_true); Yp = double(R.Y_pred);

fprintf('%-9s %-9s %-5s %-14s %-14s\n', 'solution', 'eps', 'role', 'resid(truth)', 'resid(FNO)');
for k = 1:numel(idx)
    T  = P(P(:, 1) == idx(k), 2);
    ep = 1 - P(P(:, 1) == idx(k), 4);
    rt = check_wave_residual(squeeze(Yt(k, :, 1, :)), squeeze(Yt(k, :, 2, :)), T);
    rp = check_wave_residual(squeeze(Yp(k, :, 1, :)), squeeze(Yp(k, :, 2, :)), T);
    if ood(k), role = 'OOD'; else, role = 'ID'; end
    fprintf('V%-8d %-9.4f %-5s %-14.3e %-14.3e\n', idx(k), ep, role, rt, rp);
    rel_fno(k) = rp; rel_tru(k) = rt;                       %#ok<SAGROW>
    % accumulate for the POOLED residual (sums before the ratio, matching
    % the energy-weighted pooling of the error tables)
    [~, Rk] = check_wave_residual(squeeze(Yp(k, :, 1, :)), squeeze(Yp(k, :, 2, :)), T);
    Mt = size(Yp, 2); dtk = (T/4) / (Mt - 1);
    pt = ( squeeze(Yp(k, 3:Mt, 2, :)) - squeeze(Yp(k, 1:Mt-2, 2, :)) ) / (2*dtk);  % crude phi_t scale
    num_R(k) = sum(Rk(:).^2);
    den_R(k) = sum(pt(:).^2) * (Mt/(Mt-2)) + sum(reshape(squeeze(Yp(k, :, 1, :)), [], 1).^2);
end
fprintf('\nPooled residuals (sums over the group before the ratio):\n');
for gg = 1:2
    if gg == 1, msk = ~ood; lbl = 'ID '; else, msk = ood; lbl = 'OOD'; end
    pooled = sqrt(sum(num_R(msk)) / sum(den_R(msk)));
    fprintf('  %s: pooled %.3e   (per-solution mean %.3e, min %.3e, max %.3e)\n', ...
            lbl, pooled, mean(rel_fno(msk)), min(rel_fno(msk)), max(rel_fno(msk)));
end

%% ---------------- exploratory predictions: fit the period -----------------
if exist(TRY_FILE, 'file')
    T5 = load(TRY_FILE);
    Ytry = double(T5.Y_pred); etry = double(T5.eps_try(:));
    T240 = P(P(:, 1) == 240, 2);          % period of the last data member
    % family-trend period at the hypothetical eps: quadratic fit of T(eps)
    % over the same solutions used to extrapolate the profiles (V231..240)
    fit_ids = double(T5.fit_indices(:));
    efit = zeros(size(fit_ids)); Tfam = zeros(size(fit_ids));
    for q = 1:numel(fit_ids)
        efit(q) = 1 - P(P(:, 1) == fit_ids(q), 4);
        Tfam(q) = P(P(:, 1) == fit_ids(q), 2);
    end
    cT = polyfit(efit, Tfam, 2);
    fprintf('\nExploratory predictions (NO truth): residual at the extrapolated period\n');
    fprintf('(T extrapolated in eps exactly as the initial profiles were)\n');
    fprintf('%-8s %-14s %-14s\n', 'eps', 'T', 'residual');
    for j = 1:numel(etry)
        e = squeeze(Ytry(j, :, 1, :)); f = squeeze(Ytry(j, :, 2, :));
        Tj = polyval(cT, etry(j));               % period from the branch trend
        fprintf('%-8.3f %-14.5f %-14.3e\n', etry(j), Tj, check_wave_residual(e, f, Tj));
    end
    fprintf('(residuals at the level of the tested FNO predictions above support\n');
    fprintf(' the extrapolation; percent-level changes of T move r only slightly)\n');
end
