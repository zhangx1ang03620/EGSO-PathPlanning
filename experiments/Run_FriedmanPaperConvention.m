% Run_FriedmanPaperConvention
% Friedman test exactly as reported in the paper (chi-sq = 79.29, df = 22,
% p = 2.12e-08) for the fixed-iteration protocol (N = 50, T = 100).
%
% Ranks are taken per map on the valid mean path length (V.Mean). Algorithms
% with a 0% success rate on a map (undefined V.Mean) are placed last, and ties
% among them are broken by the listing order of the algorithm table (the order
% of Algorithm_Names). These are the ranks of the ranking table in the paper,
% so the statistic can be reproduced from that table alone.
%
% Run_SignificanceTests_v2.m gives these failed algorithms the average of the
% tied ranks instead, which yields chi-sq = 74.30 on the same data.
%
% Run from the repository root after setup_egso.

R = load(fullfile('results', 'main', 'EGSO_Comparison_Results.mat'));
names = R.Algorithm_Names;
k = numel(names);
ranks = zeros(5, k);
for m = 1:5
    v = inf(1, k);
    for a = 1:k
        x = R.Results(m).(names{a}).Mean;
        if isfinite(x), v(a) = x; end
    end
    [~, order] = sort(v);            % stable: ties keep the listing order
    ranks(m, order) = 1:k;
end
[chi2, df, p] = ManualFriedman(ranks);
fprintf('Ranks (maps x algorithms):\n');
for a = 1:k
    fprintf('  %-5s %s  avg %.1f\n', names{a}, sprintf('%3d', ranks(:, a)), mean(ranks(:, a)));
end
fprintf('Friedman chi-sq = %.2f, df = %d, p = %.3g\n', chi2, df, p);
