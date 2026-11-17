% =========================================================================
% EGSO ablation study
% Contribution of each of the four core components:
%   1. EGSO (Full)          - full algorithm
%   2. EGSO w/o EVF          - hand-crafted APF instead of the Eikonal field
%   3. EGSO w/o GSI          - uniform random feasible initialisation
%   4. EGSO w/o CALS         - plain Laplacian smoothing without clearance
%   5. EGSO w/o AntiStag     - no anti-stagnation restart
%
% Maps: Map 1 (decoy dead-end trap) / Map 3 (Voronoi rooms) / Map 5 (warehouse shelves),
% covering trap-sensitive, generic and large dense topologies
% =========================================================================
clc; clear; close all;

if ~exist('results/ablation', 'dir')
    mkdir('results/ablation');
end

fprintf('========================================\n');
fprintf('EGSO ablation study\n');
fprintf('========================================\n\n');

[Maps, Map_Info] = GenerateEGSOMaps();
close(gcf);

%% Experiment parameters (paper setting)
max_gen        = 100;
num_population = 50;
num_runs       = 10;
collision_threshold = 1e4;

Ablation_Names = {'EGSO', 'w/o EVF', 'w/o GSI', 'w/o CALS', 'w/o AntiStag'};
num_variants = length(Ablation_Names);

abl_colors = [
    0.90 0.10 0.10;   % EGSO Full     - red
    0.00 0.45 0.74;   % w/o EVF       - blue
    0.00 0.70 0.30;   % w/o GSI       - green
    1.00 0.55 0.00;   % w/o CALS      - orange
    0.55 0.00 0.55;   % w/o AntiStag  - purple
];
abl_markers = {'o', 's', '^', 'd', 'v'};
abl_lines   = {'-', '--', '-.', ':', '-'};

fprintf('Configuration:\n');
fprintf('  variants: %d\n', num_variants);
fprintf('  maps: 3 (Map 1/3/5)\n');
fprintf('  runs: %d\n', num_runs);
fprintf('  iterations: %d\n', max_gen);
fprintf('  population: %d\n', num_population);

%% Result storage
AblResults = struct();
map_ids = [1, 3, 5];
for map_idx = map_ids
    for v = 1:num_variants
        AblResults(map_idx).Fitness{v} = zeros(num_runs, 1);
        AblResults(map_idx).Curves{v}  = zeros(num_runs, max_gen);
        AblResults(map_idx).Time{v}    = zeros(num_runs, 1);
    end
end

%% Main loop
for map_idx = map_ids
    fprintf('\n========== Map %d/5: %s ==========\n', map_idx, Map_Info(map_idx).size);

    G = Maps{map_idx};
    [Xmax, ~] = size(G);
    X_min = 1;
    dimensions = Xmax - 2;
    fobj = @(x) fitness(x, G);

    for run = 1:num_runs
        fprintf('  Run %d/%d ...', run, num_runs);

        for v = 1:num_variants
            tic;
            switch v
                case 1
                    [Best_score, ~, Curve] = EGSO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                case 2
                    [Best_score, ~, Curve] = EGSO_noEVF(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                case 3
                    [Best_score, ~, Curve] = EGSO_noGSI(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                case 4
                    [Best_score, ~, Curve] = EGSO_noCALS(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                case 5
                    [Best_score, ~, Curve] = EGSO_noAntiStag(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
            end
            elapsed = toc;

            AblResults(map_idx).Fitness{v}(run)   = Best_score;
            AblResults(map_idx).Curves{v}(run, :) = Curve;
            AblResults(map_idx).Time{v}(run)      = elapsed;
        end

        fprintf(' done\n');
    end
end

%% Statistics and output
fprintf('\n\n========================================\n');
fprintf('Ablation results\n');
fprintf('========================================\n');

fid = fopen('results/ablation/Ablation_Results.txt', 'w');
fprintf(fid, 'EGSO Ablation Study Results\n');
fprintf(fid, 'Config: %d runs, %d iterations, pop=%d\n', num_runs, max_gen, num_population);
fprintf(fid, '=========================================================\n\n');

Summary_SR    = zeros(num_variants, 5);
Summary_VMean = zeros(num_variants, 5);
Summary_VStd  = zeros(num_variants, 5);
Summary_VBest = zeros(num_variants, 5);

for map_idx = map_ids
    fprintf('\n===== Map %d: %s =====\n', map_idx, Map_Info(map_idx).size);
    fprintf(fid, '===== Map %d: %s =====\n', map_idx, Map_Info(map_idx).size);

    header = sprintf('%-18s %7s %10s %10s %10s %10s %8s', ...
                     'Variant', 'SR(%)', 'V.Best', 'V.Mean', 'V.Worst', 'V.Std', 'Time(s)');
    fprintf('%s\n', header);
    fprintf(fid, '%s\n', header);
    fprintf('%s\n', repmat('-', 1, 80));
    fprintf(fid, '%s\n', repmat('-', 1, 80));

    for v = 1:num_variants
        all_fit = AblResults(map_idx).Fitness{v};
        valid_mask = all_fit < collision_threshold;
        valid_fit = all_fit(valid_mask);
        sr = sum(valid_mask) / num_runs * 100;
        avg_time = mean(AblResults(map_idx).Time{v});

        Summary_SR(v, map_idx) = sr;

        if ~isempty(valid_fit)
            v_best = min(valid_fit); v_mean = mean(valid_fit);
            v_worst = max(valid_fit); v_std = std(valid_fit);
            Summary_VMean(v, map_idx) = v_mean;
            Summary_VStd(v, map_idx)  = v_std;
            Summary_VBest(v, map_idx) = v_best;
            line_str = sprintf('%-18s %6.1f%% %10.2f %10.2f %10.2f %10.2f %8.2f', ...
                              Ablation_Names{v}, sr, v_best, v_mean, v_worst, v_std, avg_time);
        else
            Summary_VMean(v, map_idx) = inf;
            Summary_VStd(v, map_idx)  = 0;
            Summary_VBest(v, map_idx) = inf;
            line_str = sprintf('%-18s %6.1f%% %10s %10s %10s %10s %8.2f', ...
                              Ablation_Names{v}, sr, '---', '---', '---', '---', avg_time);
        end
        fprintf('%s\n', line_str);
        fprintf(fid, '%s\n', line_str);
    end
    fprintf(fid, '\n');
end

%% Degradation table
fprintf('\n\n===== Performance Degradation Analysis (V.Mean increase %%) =====\n');
fprintf(fid, '\n===== Performance Degradation Analysis (V.Mean increase %%) =====\n');

deg_header = sprintf('%-18s', 'Variant');
for m = map_ids
    deg_header = [deg_header, sprintf(' %8s', sprintf('Map%d', m))];
end
deg_header = [deg_header, sprintf(' %8s', 'Avg')];
fprintf('%s\n', deg_header);
fprintf(fid, '%s\n', deg_header);
fprintf('%s\n', repmat('-', 1, 60));
fprintf(fid, '%s\n', repmat('-', 1, 60));

for v = 2:num_variants
    line_str = sprintf('%-18s', Ablation_Names{v});
    deg_vals = [];
    for m = map_ids
        if Summary_VMean(1, m) > 0 && ~isinf(Summary_VMean(v, m))
            deg = (Summary_VMean(v, m) - Summary_VMean(1, m)) / Summary_VMean(1, m) * 100;
            line_str = [line_str, sprintf(' %+7.1f%%', deg)];
            deg_vals = [deg_vals, deg];
        else
            line_str = [line_str, sprintf(' %8s', 'FAIL')];
        end
    end
    if ~isempty(deg_vals)
        line_str = [line_str, sprintf(' %+7.1f%%', mean(deg_vals))];
    end
    fprintf('%s\n', line_str);
    fprintf(fid, '%s\n', line_str);
end

%% SR degradation
fprintf('\n\n===== SR Degradation Analysis =====\n');
fprintf(fid, '\n===== SR Degradation Analysis =====\n');

sr_header = sprintf('%-18s', 'Variant');
for m = map_ids
    sr_header = [sr_header, sprintf(' %8s', sprintf('Map%d', m))];
end
fprintf('%s\n', sr_header);
fprintf(fid, '%s\n', sr_header);
fprintf('%s\n', repmat('-', 1, 50));
fprintf(fid, '%s\n', repmat('-', 1, 50));

for v = 1:num_variants
    line_str = sprintf('%-18s', Ablation_Names{v});
    for m = map_ids
        line_str = [line_str, sprintf(' %7.1f%%', Summary_SR(v, m))];
    end
    fprintf('%s\n', line_str);
    fprintf(fid, '%s\n', line_str);
end

fclose(fid);
fprintf('\nSaved: results/ablation/Ablation_Results.txt\n');

%% Convergence curves
for map_idx = map_ids
    fig = figure('Position', [100, 100, 900, 550], 'Visible', 'off');
    hold on;
    plot_handles = [];
    legend_entries = {};

    for v = 1:num_variants
        valid_mask = AblResults(map_idx).Fitness{v} < collision_threshold;
        if sum(valid_mask) > 0
            avg_curve = mean(AblResults(map_idx).Curves{v}(valid_mask, :), 1);
        else
            avg_curve = mean(AblResults(map_idx).Curves{v}, 1);
        end

        if v == 1
            lw = 3.0;
        else
            lw = 1.8;
        end

        h = plot(avg_curve, 'Color', abl_colors(v, :), 'LineWidth', lw, ...
                 'LineStyle', abl_lines{v}, 'Marker', abl_markers{v}, ...
                 'MarkerIndices', 1:max(1,round(max_gen/10)):max_gen, 'MarkerSize', 7);
        plot_handles = [plot_handles, h];
        legend_entries{end+1} = Ablation_Names{v};
    end

    xlabel('Iteration', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Fitness Value', 'FontSize', 13, 'FontWeight', 'bold');
    title(sprintf('Ablation Study - Map %d (%s)', map_idx, Map_Info(map_idx).size), ...
          'FontSize', 14, 'FontWeight', 'bold');
    legend(plot_handles, legend_entries, 'Location', 'northeast', 'FontSize', 10);
    grid on;
    set(gca, 'YScale', 'log');
    set(gcf, 'Color', 'w');
    hold off;

    print(fig, sprintf('results/ablation/Ablation_Convergence_Map%d', map_idx), '-dpng', '-r300');
    close(fig);
end

%% SR bar chart
fig_sr = figure('Position', [100, 100, 1000, 500], 'Visible', 'off');
bar_data = Summary_SR(:, map_ids)';
b = bar(bar_data, 'grouped');
for v = 1:num_variants
    b(v).FaceColor = abl_colors(v, :);
end
set(gca, 'XTickLabel', {'Map1', 'Map3', 'Map5'}, 'FontSize', 12);
ylabel('Success Rate (%)', 'FontSize', 13, 'FontWeight', 'bold');
title('Ablation Study - Collision-Free Success Rate', 'FontSize', 14, 'FontWeight', 'bold');
legend(Ablation_Names, 'Location', 'southwest', 'FontSize', 10);
ylim([0 105]);
grid on;
set(gcf, 'Color', 'w');
print(fig_sr, 'results/ablation/Ablation_SR_BarChart', '-dpng', '-r300');
close(fig_sr);

%% V.Mean bar chart
fig_vm = figure('Position', [100, 100, 1000, 500], 'Visible', 'off');
vm_data = Summary_VMean(:, map_ids)';
vm_data(isinf(vm_data)) = 0;
b2 = bar(vm_data, 'grouped');
for v = 1:num_variants
    b2(v).FaceColor = abl_colors(v, :);
end
set(gca, 'XTickLabel', {'Map1', 'Map3', 'Map5'}, 'FontSize', 12);
ylabel('Valid Mean Fitness', 'FontSize', 13, 'FontWeight', 'bold');
title('Ablation Study - Path Quality Comparison (V.Mean)', 'FontSize', 14, 'FontWeight', 'bold');
legend(Ablation_Names, 'Location', 'northwest', 'FontSize', 10);
grid on;
set(gcf, 'Color', 'w');
print(fig_vm, 'results/ablation/Ablation_VMean_BarChart', '-dpng', '-r300');
close(fig_vm);

%% Save data
save('results/ablation/Ablation_Data.mat', 'AblResults', 'Map_Info', ...
     'Ablation_Names', 'Summary_SR', 'Summary_VMean', 'Summary_VStd', ...
     'num_runs', 'max_gen', 'num_population');
fprintf('\nSaved: results/ablation/Ablation_Data.mat\n');

%% LaTeX table
fid_latex = fopen('results/ablation/Ablation_Table_LaTeX.txt', 'w');
fprintf(fid_latex, '%% EGSO Ablation Study Results - LaTeX Table\n\n');
fprintf(fid_latex, '\\begin{table}[htbp]\n\\centering\n');
fprintf(fid_latex, '\\caption{Ablation study results of EGSO across 3 representative maps (%d runs)}\n', num_runs);
fprintf(fid_latex, '\\begin{tabular}{l');
for m = map_ids
    fprintf(fid_latex, 'rr');
end
fprintf(fid_latex, '}\n\\hline\n');
fprintf(fid_latex, 'Variant');
for m = map_ids
    fprintf(fid_latex, ' & \\multicolumn{2}{c}{Map %d}', m);
end
fprintf(fid_latex, ' \\\\\n');
fprintf(fid_latex, ' ');
for m = map_ids
    fprintf(fid_latex, ' & SR(\\%%) & V.Mean');
end
fprintf(fid_latex, ' \\\\\n\\hline\n');

for v = 1:num_variants
    if v == 1
        fprintf(fid_latex, '\\textbf{%s}', Ablation_Names{v});
    else
        fprintf(fid_latex, '%s', Ablation_Names{v});
    end
    for m = map_ids
        if isinf(Summary_VMean(v,m))
            fprintf(fid_latex, ' & %.1f & ---', Summary_SR(v,m));
        else
            fprintf(fid_latex, ' & %.1f & %.2f', Summary_SR(v,m), Summary_VMean(v,m));
        end
    end
    fprintf(fid_latex, ' \\\\\n');
end
fprintf(fid_latex, '\\hline\n\\end{tabular}\n\\end{table}\n');
fclose(fid_latex);
fprintf('Saved: results/ablation/Ablation_Table_LaTeX.txt\n');

fprintf('\n========================================\n');
fprintf('Ablation study finished. All results saved in results/ablation/\n');
fprintf('========================================\n');
