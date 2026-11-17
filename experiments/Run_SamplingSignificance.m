% =========================================================================
% Run_SamplingSignificance.m
% Wilcoxon rank-sum tests (own implementation, ManualRankSum.m; no Statistics Toolbox):
% EGSO vs RRT* / Informed-RRT*, from the saved raw data of 20 independent runs,
% collision-free runs only, alpha = 0.05.
% =========================================================================
clc; clear;

R = load('results/main/EGSO_Comparison_Results.mat');
S = load('results/sampling/Sampling_Baselines.mat');
collision_threshold = 1e4;
alpha = 0.05;

fid = fopen('results/sampling/Sampling_Significance.txt', 'w');
fprintf(fid, 'Wilcoxon Rank-Sum Test (manual implementation): EGSO vs RRT* / Informed-RRT* (alpha=0.05)\n');
fprintf(fid, '=========================================================================\n\n');

names = {'RRTStar', 'InformedRRTStar'};
labels = {'RRT*', 'Informed-RRT*'};
WTL = zeros(2,3);

for map_idx = 1:5
    egso_all = R.Results(map_idx).EGSO.Best;
    egso_valid = egso_all(egso_all < collision_threshold);
    line = sprintf('===== Map %d =====', map_idx);
    fprintf('%s\n', line); fprintf(fid, '%s\n', line);

    for k = 1:2
        lens = S.Sampling(map_idx).(names{k}).AllLens;
        other_valid = lens(isfinite(lens));

        [p, W, z] = ManualRankSum(egso_valid, other_valid);
        med_egso = median(egso_valid);
        med_other = median(other_valid);
        if p < alpha
            if med_egso < med_other
                verdict = '+ (EGSO better)'; WTL(k,1) = WTL(k,1)+1;
            else
                verdict = '- (EGSO worse)'; WTL(k,3) = WTL(k,3)+1;
            end
        else
            verdict = '~ (no sig. diff.)'; WTL(k,2) = WTL(k,2)+1;
        end
        line = sprintf('  EGSO vs %-16s: median EGSO=%.2f, median %s=%.2f, z=%.3f, p=%.4g  %s', ...
                labels{k}, med_egso, labels{k}, med_other, z, p, verdict);
        fprintf('%s\n', line); fprintf(fid, '%s\n', line);
    end
end

fprintf(fid, '\nSummary (Win/Tie/Loss across 5 maps):\n');
fprintf('\nSummary (Win/Tie/Loss across 5 maps):\n');
for k = 1:2
    line = sprintf('  EGSO vs %-16s: Win=%d Tie=%d Loss=%d', labels{k}, WTL(k,1), WTL(k,2), WTL(k,3));
    fprintf('%s\n', line); fprintf(fid, '%s\n', line);
end
fclose(fid);
fprintf('\nSaved: results/sampling/Sampling_Significance.txt\n');
