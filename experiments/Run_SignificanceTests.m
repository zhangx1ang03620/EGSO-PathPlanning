% =========================================================================
% Run_SignificanceTests.m
% Statistical significance tests on the saved raw fitness of 20 independent runs:
%   1. Wilcoxon rank-sum test: on every map, EGSO vs each of the other 22 algorithms,
%      on the fitness of collision-free (valid) runs only (alpha = 0.05)
%   2. Friedman test on the 5 maps x 23 algorithms V.Mean matrix
% =========================================================================
clc; clear;

if ~exist('results/main/EGSO_Comparison_Results.mat', 'file')
    error('Run_SignificanceTests:Missing', 'run Main_Comparison_Experiment / Aggregate_Comparison_Results first to create results/main/EGSO_Comparison_Results.mat');
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

%% ===== 1. Wilcoxon rank-sum: EGSO vs the other 22 algorithms on every map =====
fid = fopen('results/main/Significance_Wilcoxon.txt', 'w');
fprintf(fid, 'Wilcoxon Rank-Sum Test: EGSO vs each baseline (alpha=0.05)\n');
fprintf(fid, 'Test conducted on collision-free (valid) runs only.\n');
fprintf(fid, '"+" = EGSO significantly better, "-" = EGSO significantly worse, "~" = no significant difference\n');
fprintf(fid, '"+ (0%%SR)" = opponent found zero collision-free paths (trivial win, no ranksum possible)\n');
fprintf(fid, '=========================================================================\n\n');

% W/T/L counts (win/tie/loss of every competitor over the five maps)
WTL = zeros(num_algorithms, 3);  % columns: [Win(EGSO better), Tie, Loss(EGSO worse)]

for map_idx = 1:5
    fprintf('\n===== Map %d: %s =====\n', map_idx, Map_Info(map_idx).size);
    fprintf(fid, '\n===== Map %d: %s =====\n', map_idx, Map_Info(map_idx).size);

    egso_all = Results(map_idx).EGSO.Best;
    egso_valid = egso_all(egso_all < collision_threshold);

    header = sprintf('%-10s %10s %10s %12s %8s', 'Algorithm', 'EGSO.Mean', 'Opp.Mean', 'p-value', 'Result');
    fprintf('%s\n', header); fprintf(fid, '%s\n', header);

    for k = other_idx
        opp_name = Algorithm_Names{k};
        opp_all = Results(map_idx).(opp_name).Best;
        opp_valid = opp_all(opp_all < collision_threshold);

        if isempty(opp_valid)
            result_str = '+ (0%SR)';
            p_val = NaN;
            opp_mean = NaN;
            WTL(k, 1) = WTL(k, 1) + 1;
        elseif isempty(egso_valid)
            result_str = 'N/A';
            p_val = NaN;
            opp_mean = mean(opp_valid);
        else
            p_val = ranksum(egso_valid, opp_valid);
            opp_mean = mean(opp_valid);
            egso_mean = mean(egso_valid);
            if p_val < alpha
                if egso_mean < opp_mean
                    result_str = '+';
                    WTL(k, 1) = WTL(k, 1) + 1;
                else
                    result_str = '-';
                    WTL(k, 3) = WTL(k, 3) + 1;
                end
            else
                result_str = '~';
                WTL(k, 2) = WTL(k, 2) + 1;
            end
        end

        if isnan(opp_mean)
            line_str = sprintf('%-10s %10.2f %10s %12s %8s', opp_name, mean(egso_valid), '---', '---', result_str);
        elseif isnan(p_val)
            line_str = sprintf('%-10s %10.2f %10.2f %12s %8s', opp_name, mean(egso_valid), opp_mean, '---', result_str);
        else
            line_str = sprintf('%-10s %10.2f %10.2f %12.4g %8s', opp_name, mean(egso_valid), opp_mean, p_val, result_str);
        end
        fprintf('%s\n', line_str); fprintf(fid, '%s\n', line_str);
    end
end

%% ===== W/T/L summary table =====
fprintf('\n\n===== Win/Tie/Loss Summary (EGSO vs each baseline, across 5 maps) =====\n');
fprintf(fid, '\n\n===== Win/Tie/Loss Summary (EGSO vs each baseline, across 5 maps) =====\n');
header2 = sprintf('%-10s %6s %6s %6s', 'Algorithm', 'Win', 'Tie', 'Loss');
fprintf('%s\n', header2); fprintf(fid, '%s\n', header2);
total_W = 0; total_T = 0; total_L = 0;
for k = other_idx
    line_str = sprintf('%-10s %6d %6d %6d', Algorithm_Names{k}, WTL(k,1), WTL(k,2), WTL(k,3));
    fprintf('%s\n', line_str); fprintf(fid, '%s\n', line_str);
    total_W = total_W + WTL(k,1); total_T = total_T + WTL(k,2); total_L = total_L + WTL(k,3);
end
overall_str = sprintf('\nOverall across all %d baselines x 5 maps = %d comparisons: EGSO Win=%d (%.1f%%), Tie=%d (%.1f%%), Loss=%d (%.1f%%)', ...
    numel(other_idx), numel(other_idx)*5, total_W, 100*total_W/(numel(other_idx)*5), ...
    total_T, 100*total_T/(numel(other_idx)*5), total_L, 100*total_L/(numel(other_idx)*5));
fprintf('%s\n', overall_str); fprintf(fid, '%s\n', overall_str);
fclose(fid);
fprintf('\nSaved: results/main/Significance_Wilcoxon.txt\n');

%% ===== 2. Friedman test: 5 maps x 23 algorithms V.Mean matrix =====
fprintf('\n\n===== Friedman Test (5 maps x %d algorithms, blocked by map) =====\n', num_algorithms);
VMean_matrix = zeros(5, num_algorithms);
for map_idx = 1:5
    for k = 1:num_algorithms
        m = Results(map_idx).(Algorithm_Names{k}).Mean;
        if isinf(m)
            % Inf is replaced by the largest valid V.Mean on that map plus a penalty margin,
            % which keeps the matrix usable and still ranks the failure last
            m = NaN;
        end
        VMean_matrix(map_idx, k) = m;
    end
end
% NaN (complete failure) is replaced by the largest finite value of its row (map), i.e. worst on that map
for map_idx = 1:5
    row = VMean_matrix(map_idx, :);
    finite_max = max(row(isfinite(row)));
    row(isnan(row)) = finite_max * 1.5;
    VMean_matrix(map_idx, :) = row;
end

[p_friedman, tbl, stats] = friedman(VMean_matrix, 1, 'off');
fprintf('Friedman chi-sq statistic = %.4f, df = %d, p-value = %.4g\n', tbl{2,5}, tbl{2,3}, p_friedman);

fid2 = fopen('results/main/Significance_Friedman.txt', 'w');
fprintf(fid2, 'Friedman Test: %d algorithms across 5 maps (blocked by map, V.Mean as score)\n', num_algorithms);
fprintf(fid2, 'Chi-sq = %.4f, df = %d, p-value = %.4g\n', tbl{2,5}, tbl{2,3}, p_friedman);
if p_friedman < alpha
    fprintf(fid2, 'Conclusion: significant overall difference among algorithms (p < 0.05).\n');
else
    fprintf(fid2, 'Conclusion: no significant overall difference detected (p >= 0.05).\n');
end
fclose(fid2);
fprintf('Saved: results/main/Significance_Friedman.txt\n');

save('results/main/Significance_Data.mat', 'WTL', 'VMean_matrix', 'p_friedman', 'tbl', 'stats', 'Algorithm_Names');

fprintf('\n========================================\n');
fprintf('Significance tests finished\n');
fprintf('========================================\n');
