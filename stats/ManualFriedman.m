function [chi2, df, p] = ManualFriedman(X)
% ManualFriedman  Self-contained Friedman test (rows = blocks/maps, columns
% = treatments/algorithms), avoiding a dependency on the Statistics and
% Machine Learning Toolbox's friedman (license entitled but not actually
% installed in this MATLAB installation -- see ManualRankSum.m for the
% same situation with ranksum). Standard formulation with average-rank tie
% handling (Lehmann 1975); chi-square upper-tail p-value via gammainc,
% which is part of core MATLAB (elfun/specfun) and requires no toolbox.
%
%   X : n-by-k matrix of scores, n blocks (maps) x k treatments (algorithms)
%   Returns chi2 statistic, df = k-1, and the p-value.

[n, k] = size(X);
R = zeros(n, k);
for i = 1:n
    row = X(i, :);
    [sortedVals, sortIdx] = sort(row);
    ranks = zeros(1, k);
    j = 1;
    while j <= k
        jj = j;
        while jj < k && sortedVals(jj+1) == sortedVals(j)
            jj = jj + 1;
        end
        avgRank = (j + jj) / 2;
        ranks(sortIdx(j:jj)) = avgRank;
        j = jj + 1;
    end
    R(i, :) = ranks;
end

Rj = sum(R, 1);
chi2 = (12 / (n*k*(k+1))) * sum(Rj.^2) - 3*n*(k+1);
df = k - 1;
p = gammainc(chi2/2, df/2, 'upper');

end
