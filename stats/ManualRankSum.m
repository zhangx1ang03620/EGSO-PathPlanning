function [p, W, z] = ManualRankSum(x, y)
% ManualRankSum  Self-contained two-sided Wilcoxon rank-sum (Mann-Whitney)
% test with tie-corrected normal approximation, avoiding a dependency on
% the Statistics and Machine Learning Toolbox's ranksum (not available in
% this MATLAB installation). Standard formulation, e.g. Lehmann (1975).
%
%   x, y : independent samples (vectors)
%   Returns the two-sided p-value, the rank sum of x (W), and the
%   (tie-corrected) z-statistic.

    n1 = numel(x); n2 = numel(y);
    combined = [x(:); y(:)];
    [sortedVals, sortIdx] = sort(combined);
    ranks = zeros(size(combined));

    i = 1;
    n = numel(combined);
    while i <= n
        j = i;
        while j < n && sortedVals(j+1) == sortedVals(i)
            j = j + 1;
        end
        avgRank = (i + j) / 2;
        ranks(sortIdx(i:j)) = avgRank;
        i = j + 1;
    end

    W = sum(ranks(1:n1));

    % tie correction term
    i = 1; tieSum = 0;
    while i <= n
        j = i;
        while j < n && sortedVals(j+1) == sortedVals(i)
            j = j + 1;
        end
        t = j - i + 1;
        tieSum = tieSum + (t^3 - t);
        i = j + 1;
    end

    mu_W = n1 * (n1 + n2 + 1) / 2;
    var_W = (n1 * n2 / 12) * ((n1 + n2 + 1) - tieSum / ((n1+n2) * (n1+n2-1)));
    var_W = max(var_W, eps);

    if W > mu_W
        z = (W - mu_W - 0.5) / sqrt(var_W);
    else
        z = (W - mu_W + 0.5) / sqrt(var_W);
    end

    p = 2 * (1 - 0.5*(1 + erf(abs(z)/sqrt(2))));
    p = min(max(p, 0), 1);
end
