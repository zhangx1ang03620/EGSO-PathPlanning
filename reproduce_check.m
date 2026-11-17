function ok = reproduce_check()
% REPRODUCE_CHECK  Recompute key numbers of the paper and compare them with
% the published values (about one minute).
%
%   >> setup_egso; reproduce_check
%
% Checks
%   1. Local-minima census (Sec. 4.2): the Eikonal field has a 0% trap rate
%      and no spurious sink on the five benchmark maps.
%   2. Regression: EGSO with fixed seeds returns exactly the values produced
%      by the implementation used for the paper.
%   3. Stored results (Sec. 4.3, 4.8): EGSO average rank 1.2 and the
%      Friedman statistic 79.29, recomputed from results/main.

root = setup_egso();
old = pwd; cleanup = onCleanup(@() cd(old));
cd(root);
ok = true;
[Maps, ~] = GenerateEGSOMaps();
close all;

%% 1. Local-minima census of the Eikonal field
fprintf('\n[1] Local-minima census of T_nav (paper: 0%% trap rate, 0 sinks on every map)\n');
for m = 1:5
    G = Maps{m}; M = size(G, 1);
    goalMask = false(size(G)); goalMask(M, M) = true;
    T = FastMarchingField(G, goalMask, 'navigate');
    [tr, ns] = GreedyDescentCensus(T, G, [M, M]);
    pass = tr == 0 && ns == 0;
    ok = ok && pass;
    fprintf('    Map %d: trap rate %.2f%%, sinks %d  %s\n', m, tr, ns, tag(pass));
end

%% 2. Bit-level regression of EGSO
fprintf('\n[2] EGSO regression (N = 50, T = 100, rng(2026))\n');
ref = [1 77.621255606113; 3 86.901364976489; 5 150.602230345583];
for k = 1:size(ref, 1)
    m = ref(k, 1); G = Maps{m}; M = size(G, 1);
    rng(2026);
    b = EGSO(50, 100, 1, M, M - 2, @(x) fitness(x, G), G);
    pass = abs(b - ref(k, 2)) < 1e-9;
    ok = ok && pass;
    fprintf('    Map %d: %.12f (reference %.12f)  %s\n', m, b, ref(k, 2), tag(pass));
end

%% 3. Stored results: ranks and Friedman statistic
fprintf('\n[3] Stored main comparison (results/main/EGSO_Comparison_Results.mat)\n');
R = load(fullfile('results', 'main', 'EGSO_Comparison_Results.mat'));
names = R.Algorithm_Names; k = numel(names);
ranks = zeros(5, k);
for m = 1:5
    v = inf(1, k);
    for a = 1:k
        x = R.Results(m).(names{a}).Mean;
        if isfinite(x), v(a) = x; end
    end
    [~, order] = sort(v);
    ranks(m, order) = 1:k;
end
egso_rank = mean(ranks(:, strcmp(names, 'EGSO')));
[chi2, ~, p] = ManualFriedman(ranks);
pass1 = abs(egso_rank - 1.2) < 1e-12;
pass2 = abs(chi2 - 79.29) < 0.005;
ok = ok && pass1 && pass2;
fprintf('    EGSO average rank %.1f (paper 1.2)  %s\n', egso_rank, tag(pass1));
fprintf('    Friedman chi-sq %.2f, p = %.3g (paper 79.29, 2.12e-08)  %s\n', chi2, p, tag(pass2));

fprintf('\n%s\n', ternary(ok, 'All checks passed.', 'Some checks FAILED.'));
if nargout == 0, clear ok; end
end

function s = tag(pass)
s = ternary(pass, '[ok]', '[FAIL]');
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
