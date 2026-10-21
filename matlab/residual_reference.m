% RESIDUAL_REFERENCE  Verify that Wilkening's stored solutions satisfy the
% governing equations, using check_wave_residual.m (which must be on the
% path or in the same folder).
%
% Usage: set CHECK_IDS below, then run. For each solution the
% script prints the relative residual of the DNO-eliminated dynamic
% equation; values of order 1e-6 (at the full 201 timeslices) mean the
% stored fields satisfy the equations to the accuracy of the time
% differencing. With TIME_STRIDE = 2 the check runs on the same 101 time
% instances as the FNO exports, so truth and predictions share one baseline.
clear;
CHECK_IDS   = 101:240;               % ALL solutions (any subset works)
TIME_STRIDE = 2;                     % 2 -> the same 101 time instances as the
                                     % operator exports (set 1 for all 201)
DATA_DIR = getenv('STANDING_WAVE_DATA');
if isempty(DATA_DIR), DATA_DIR = fullfile('..', 'data'); end
P = load(fullfile(DATA_DIR, 'parameters'), '-ascii');
% columns: idx | T | h | A_c | max|eta_x| | ||eta_x||^2

fprintf('%-8s %-10s %-12s\n', 'solution', 'eps', 'rel residual');
nC   = numel(CHECK_IDS);
epsA = zeros(nC, 1);  relA = zeros(nC, 1);  cnt = 0;
for i = CHECK_IDS
    % --- read the quarter-solution text file: 201 x 257 x (eta, phi, mu) ---
    vals = sscanf(fileread(fullfile(DATA_DIR, sprintf('V%d', i))), '%f');
    B    = reshape(vals, 3, 257, 201);
    eta  = squeeze(B(1, :, :))';               % (201, 257) on [0, pi]
    phi  = squeeze(B(2, :, :))';

    % --- even mirror onto the 512-point grid (-pi, pi] ---
    eta = [eta(:, end-1:-1:2), eta];
    phi = [phi(:, end-1:-1:2), phi];

    % --- restrict to the same time instances as the operator exports ---
    eta = eta(1:TIME_STRIDE:end, :);
    phi = phi(1:TIME_STRIDE:end, :);

    % --- period and epsilon from the parameters table ---
    T   = P(P(:, 1) == i, 2);
    eps = 1 - P(P(:, 1) == i, 4);

    rel = check_wave_residual(eta, phi, T);
    fprintf('V%-7d %-10.4f %-12.3e\n', i, eps, rel);
    cnt = cnt + 1;  epsA(cnt) = eps;  relA(cnt) = rel;
end

%% ---------------- summary ------------------------------------------------
fprintf('\nSummary over %d solutions:\n', cnt);
fprintf('  mean %.3e   median %.3e   min %.3e   max %.3e\n', ...
        mean(relA), median(relA), min(relA), max(relA));
