% =========================================================================
% Run_SignificanceTests_v2.m
% Significance tests without the Statistics Toolbox (the values reported in the paper):
%
%   1. Uses the toolbox-free ManualRankSum.m / ManualFriedman.m, so ranksum/friedman
%      from the Statistics Toolbox are not required.
%
%
%   2. Friedman statistic: recomputed from the same VMean_matrix as
%      Run_SignificanceTests.m, and cross-checked against the value obtained from the
%      published rank matrix alone.
%   3. Additional analysis including failed runs: besides the collision-free-only
%      Wilcoxon test, the same test on all 20 raw runs (collision penalties included),
%      accounting for feasibility and solution quality together. No new runs are
%      needed: Results(map_idx).(Algo).Best stores all 20 raw values.
%
% =========================================================================
clc; clear;

if ~exist('results/main/EGSO_Comparison_Results.mat', 'file')
    error('Run_SignificanceTests_v2:Missing', 'run Main_Comparison_Experiment / Aggregate_Comparison_Results first to create results/main/EGSO_Comparison_Results.mat');
end

R = load('results/main/EGSO_Comparison_Results.mat');
Results = R.Results;
Algorithm_Names = R.Algorithm_Names;
Map_Info = R.Map_Info;
num_algorithms = numel(Algorithm_Names);
collision_threshold = 1e4;
alpha = 0.05;

egso_idx = find(strcmp(Algorithm_Names, 'EGSO'));
other_idx = setdiff(1:num_algorithms, egso_idx);

%% ===== 1. Wilcoxon (valid runs only; recomputed with ManualRankSum) =====
WTL_valid = zeros(num_algorithms, 3);
%% ===== 1b. Wilcoxon (all 20 raw runs, collision penalties included) =====
WTL_all   = zeros(num_algorithms, 3);

fid = fopen('results/main/Significance_v2_Wilcoxon.txt', 'w');
fprintf(fid, 'Wilcoxon Rank-Sum Test (recomputed with ManualRankSum.m, no toolbox dependency)\n');
fprintf(fid, 'Two variants reported side by side for each map x baseline:\n');
fprintf(fid, '  [Valid-only]  collision-free runs only (as in Run_SignificanceTests.m)\n');
fprintf(fid, '  [All-runs]    all 20 raw runs per algorithm, including the large\n');
fprintf(fid, '                collision-penalty fitness values for failed runs -- jointly reflects\n');
fprintf(fid, '                feasibility (success rate) and solution quality in one test\n');
fprintf(fid, '=========================================================================\n\n');

for map_idx = 1:5
    fprintf('\n===== Map %d: %s =====\n', map_idx, Map_Info(map_idx).size);
    fprintf(fid, '\n===== Map %d: %s =====\n', map_idx, Map_Info(map_idx).size);

    egso_all_raw = Results(map_idx).EGSO.Best;
    egso_valid   = egso_all_raw(egso_all_raw < collision_threshold);

    header = sprintf('%-10s | %-28s | %-28s', 'Algorithm', 'Valid-only (orig. protocol)', 'All-runs (incl. failures, NEW)');
    fprintf('%s\n', header); fprintf(fid, '%s\n', header);

    for k = other_idx
        opp_name = Algorithm_Names{k};
        opp_all_raw = Results(map_idx).(opp_name).Best;
        opp_valid   = opp_all_raw(opp_all_raw < collision_threshold);

        % --- Valid-only (collision-free runs) ---
        if isempty(opp_valid)
            res_valid = '+ (0%SR)'; WTL_valid(k,1) = WTL_valid(k,1) + 1;
        elseif isempty(egso_valid)
            res_valid = 'N/A';
        else
            p_val = ManualRankSum(egso_valid, opp_valid);
            if p_val < alpha
                if mean(egso_valid) < mean(opp_valid)
                    res_valid = '+'; WTL_valid(k,1) = WTL_valid(k,1) + 1;
                else
                    res_valid = '-'; WTL_valid(k,3) = WTL_valid(k,3) + 1;
                end
            else
                res_valid = '~'; WTL_valid(k,2) = WTL_valid(k,2) + 1;
            end
        end

        % --- All-runs, including collision-penalty values (NEW) ---
        p_val_all = ManualRankSum(egso_all_raw, opp_all_raw);
        if p_val_all < alpha
            if mean(egso_all_raw) < mean(opp_all_raw)
                res_all = '+'; WTL_all(k,1) = WTL_all(k,1) + 1;
            else
                res_all = '-'; WTL_all(k,3) = WTL_all(k,3) + 1;
            end
        else
            res_all = '~'; WTL_all(k,2) = WTL_all(k,2) + 1;
        end

        line_str = sprintf('%-10s | %-28s | %-28s', opp_name, res_valid, sprintf('%s (p=%.3g)', res_all, p_val_all));
        fprintf('%s\n', line_str); fprintf(fid, '%s\n', line_str);
    end
end

fprintf('\n\n===== W/T/L Summary =====\n');
fprintf(fid, '\n\n===== W/T/L Summary =====\n');
totW1=0;totT1=0;totL1=0; totW2=0;totT2=0;totL2=0;
header2 = sprintf('%-10s | %-16s | %-16s', 'Algorithm', 'Valid-only W/T/L', 'All-runs W/T/L');
fprintf('%s\n', header2); fprintf(fid, '%s\n', header2);
for k = other_idx
    line_str = sprintf('%-10s | %2d/%2d/%2d          | %2d/%2d/%2d', Algorithm_Names{k}, ...
        WTL_valid(k,1),WTL_valid(k,2),WTL_valid(k,3), WTL_all(k,1),WTL_all(k,2),WTL_all(k,3));
    fprintf('%s\n', line_str); fprintf(fid, '%s\n', line_str);
    totW1=totW1+WTL_valid(k,1); totT1=totT1+WTL_valid(k,2); totL1=totL1+WTL_valid(k,3);
    totW2=totW2+WTL_all(k,1);   totT2=totT2+WTL_all(k,2);   totL2=totL2+WTL_all(k,3);
end
n_comp = numel(other_idx)*5;
summary_str = sprintf('\nOverall (%d comparisons): Valid-only Win=%d(%.1f%%) Tie=%d(%.1f%%) Loss=%d(%.1f%%)  |  All-runs Win=%d(%.1f%%) Tie=%d(%.1f%%) Loss=%d(%.1f%%)', ...
    n_comp, totW1,100*totW1/n_comp, totT1,100*totT1/n_comp, totL1,100*totL1/n_comp, ...
    totW2,100*totW2/n_comp, totT2,100*totT2/n_comp, totL2,100*totL2/n_comp);
fprintf('%s\n', summary_str); fprintf(fid, '%s\n', summary_str);
fclose(fid);
fprintf('\nSaved: results/main/Significance_v2_Wilcoxon.txt\n');

%% ===== 2. Friedman (ManualFriedman, no toolbox) =====
VMean_matrix = zeros(5, num_algorithms);
for map_idx = 1:5
    for k = 1:num_algorithms
        m = Results(map_idx).(Algorithm_Names{k}).Mean;
        if isinf(m) || isnan(m)
            VMean_matrix(map_idx, k) = NaN;
        else
            VMean_matrix(map_idx, k) = m;
        end
    end
end
for map_idx = 1:5
    row = VMean_matrix(map_idx, :);
    finite_max = max(row(isfinite(row)));
    row(isnan(row)) = finite_max * 1.5;
    VMean_matrix(map_idx, :) = row;
end

[chi2, df, p_friedman] = ManualFriedman(VMean_matrix);
fprintf('\n\n===== Friedman Test (ManualFriedman, %d algorithms x 5 maps) =====\n', num_algorithms);
fprintf('Chi-sq = %.4f, df = %d, p-value = %.4g\n', chi2, df, p_friedman);

fid2 = fopen('results/main/Significance_v2_Friedman.txt', 'w');
fprintf(fid2, 'Friedman Test (recomputed with ManualFriedman.m, no toolbox dependency)\n');
fprintf(fid2, '%d algorithms across 5 maps (blocked by map, V.Mean as score)\n', num_algorithms);
fprintf(fid2, 'Chi-sq = %.4f, df = %d, p-value = %.4g\n', chi2, df, p_friedman);
fprintf(fid2, '\nNote: here algorithms with 0%% SR on a map share the average of the tied ranks.\n');
fprintf(fid2, 'The paper breaks these ties by the listing order of the algorithm table, which gives\n');
fprintf(fid2, 'chi-sq = 79.29 (p = 2.12e-08); see experiments/Run_FriedmanPaperConvention.m.\n');
if p_friedman < alpha
    fprintf(fid2, 'Conclusion: significant overall difference among algorithms (p < 0.05).\n');
end
fclose(fid2);
fprintf('Saved: results/main/Significance_v2_Friedman.txt\n');

save('results/main/Significance_v2_Data.mat', 'WTL_valid', 'WTL_all', 'VMean_matrix', 'chi2', 'df', 'p_friedman', 'Algorithm_Names');

fprintf('\n========================================\n');
fprintf('Significance tests (v2) finished\n');
fprintf('========================================\n');
