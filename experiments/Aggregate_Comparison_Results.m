% =========================================================================
% Aggregate_Comparison_Results.m
% Aggregates results/main/partial_map<1..5>.mat produced by five independent processes
% (Run_Comparison_SingleMap), rebuilds the Results structure of Main_Comparison_Experiment.m,
% and writes the same statistics tables, convergence curves, path figures and LaTeX tables.
% =========================================================================
clc; clear; close all;

if ~exist('results/main', 'dir')
    mkdir('results/main');
end

collision_threshold = 1e4;

[Maps, Map_Info] = GenerateEGSOMaps();
print(gcf, 'results/main/EGSO_Test_Maps', '-dpng', '-r300');
close(gcf);

Results = struct();
Algorithm_Names = {};
max_gen = 0; num_population = 0; num_runs = 0;

for map_idx = 1:5
    fname = sprintf('results/main/partial_map%d.mat', map_idx);
    if ~exist(fname, 'file')
        error('Aggregate:MissingFile', 'Missing %s -- run Run_Comparison_SingleMap(%d) first.', fname, map_idx);
    end
    S = load(fname);
    Algorithm_Names   = S.Algorithm_Names;
    max_gen           = S.max_gen;
    num_population    = S.num_population;
    num_runs          = S.num_runs;

    Results(map_idx).Map_Name = sprintf('Map%d', map_idx);
    Results(map_idx).Map_Size = Map_Info(map_idx).size;

    for algo_idx = 1:numel(Algorithm_Names)
        algo_name = Algorithm_Names{algo_idx};
        Results(map_idx).(algo_name).Best   = S.MapResult.(algo_name).Best;
        Results(map_idx).(algo_name).Curves = S.MapResult.(algo_name).Curves;
        Results(map_idx).(algo_name).Paths  = S.MapResult.(algo_name).Paths;

        % ===== Robustness fix: filter wall-clock outliers caused by SIGSTOP/SIGCONT pauses =====
        % The original runs were paused overnight with SIGSTOP and resumed with SIGCONT; the
        % tic/toc of the run executing at that moment includes the pause (>> 100 s), whereas
        % every genuine single run takes < 60 s. Only the time column is affected, not the
        % fitness or paths. Outliers are replaced by the median time of the remaining runs of
        % the same algorithm on the same map.
        raw_time = S.MapResult.(algo_name).Time;
        time_outlier_threshold = 100;  % seconds; every genuine single run is far below this
        is_outlier = raw_time > time_outlier_threshold;
        if any(is_outlier) && ~all(is_outlier)
            raw_time(is_outlier) = median(raw_time(~is_outlier));
        end
        Results(map_idx).(algo_name).Time = raw_time;

        all_fitness = Results(map_idx).(algo_name).Best;
        valid_mask  = all_fitness < collision_threshold;
        valid_fitness = all_fitness(valid_mask);

        if ~isempty(valid_fitness)
            Results(map_idx).(algo_name).Mean  = mean(valid_fitness);
            Results(map_idx).(algo_name).Worst = max(valid_fitness);
            Results(map_idx).(algo_name).Std   = std(valid_fitness);
            Results(map_idx).(algo_name).SR    = sum(valid_mask) / num_runs * 100;
        else
            Results(map_idx).(algo_name).Mean  = inf;
            Results(map_idx).(algo_name).Worst = inf;
            Results(map_idx).(algo_name).Std   = 0;
            Results(map_idx).(algo_name).SR    = 0;
        end
    end
end

num_algorithms = numel(Algorithm_Names);
algo_colors = [1.00 0.05 0.05; hsv(num_algorithms-1)*0.85];
markers_pool = {'o', 's', '^', 'd', 'v', 'p', 'h', '<', '>', 'x', ...
                '+', '*', 'pentagram', 'hexagram'};
num_markers = length(markers_pool);

fprintf('Loaded all 5 partial result files. Config: %d algos, %d runs, %d iters, pop=%d\n', ...
        num_algorithms, num_runs, max_gen, num_population);

%% Save the experiment data
save('results/main/EGSO_Comparison_Results.mat', 'Results', 'Map_Info', ...
     'Algorithm_Names', 'num_runs', 'max_gen', 'num_population');
fprintf('Saved: results/main/EGSO_Comparison_Results.mat\n');

%% Convergence curves / path figures / overlays (per map)
for map_idx = 1:5
    fig_conv = figure('Name', sprintf('Map%d Convergence', map_idx), ...
                      'Position', [50, 50, 1200, 700], 'Visible', 'off');
    hold on;
    legend_entries = {};
    plot_handles = [];

    for algo_idx = 1:num_algorithms
        algo_name = Algorithm_Names{algo_idx};
        valid_mask = Results(map_idx).(algo_name).Best < collision_threshold;

        if sum(valid_mask) > 0
            avg_curve = mean(Results(map_idx).(algo_name).Curves(valid_mask, :), 1);
        else
            avg_curve = mean(Results(map_idx).(algo_name).Curves, 1);
        end

        if strcmp(algo_name, 'EGSO')
            h = plot(avg_curve, 'Color', algo_colors(1, :), 'LineWidth', 3);
        else
            marker_idx = mod(algo_idx - 2, num_markers) + 1;
            h = plot(avg_curve, 'Color', algo_colors(algo_idx, :), 'LineWidth', 1.2, ...
                     'Marker', markers_pool{marker_idx}, ...
                     'MarkerIndices', 1:max(1,round(max_gen/10)):max_gen, ...
                     'MarkerSize', 5);
        end
        plot_handles = [plot_handles, h];
        legend_entries{end+1} = algo_name;
    end

    xlabel('Iteration', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Fitness Value (log scale)', 'FontSize', 13, 'FontWeight', 'bold');
    title(sprintf('Map %d (%s) - Convergence Curves', map_idx, Map_Info(map_idx).size), ...
          'FontSize', 14, 'FontWeight', 'bold');
    legend(plot_handles, legend_entries, 'Location', 'northeast', 'NumColumns', 4, 'FontSize', 7);
    grid on;
    set(gca, 'YScale', 'log');
    set(gcf, 'Color', 'w');
    hold off;
    print(fig_conv, sprintf('results/main/Convergence_Map%d', map_idx), '-dpng', '-r300');
    close(fig_conv);

    fig_path = figure('Name', sprintf('Path Map%d', map_idx), ...
                      'Position', [50, 50, 1800, 1500], 'Visible', 'off');
    for algo_idx = 1:num_algorithms
        algo_name = Algorithm_Names{algo_idx};
        subplot(5, 5, algo_idx);
        imagesc(Maps{map_idx});
        colormap(gca, [1 1 1; 0.7 0.5 0.8]);
        axis equal tight;
        hold on;

        valid_runs_idx = find(Results(map_idx).(algo_name).Best < collision_threshold);
        if ~isempty(valid_runs_idx)
            [~, best_run_local] = min(Results(map_idx).(algo_name).Best(valid_runs_idx));
            best_run = valid_runs_idx(best_run_local);
            path = Results(map_idx).(algo_name).Paths{best_run};
            if ~isempty(path)
                plot(path(:, 2), path(:, 1), '-', 'Color', algo_colors(algo_idx, :), 'LineWidth', 2);
                plot(path(1, 2), path(1, 1), 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g');
                plot(path(end, 2), path(end, 1), 'rs', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
            end
            best_val = min(Results(map_idx).(algo_name).Best(valid_runs_idx));
            title_str = sprintf('%s (%.1f)', algo_name, best_val);
        else
            title_str = sprintf('%s (FAIL)', algo_name);
        end
        title(title_str, 'FontSize', 9, 'FontWeight', 'bold');
        set(gca, 'YDir', 'normal', 'XTick', [], 'YTick', []);
        hold off;
    end
    sgtitle(sprintf('Map %d (%s) - Best Path Comparison (%d Algorithms)', ...
            map_idx, Map_Info(map_idx).size, num_algorithms), 'FontSize', 14, 'FontWeight', 'bold');
    set(gcf, 'Color', 'w');
    print(fig_path, sprintf('results/main/Paths_Map%d', map_idx), '-dpng', '-r300');
    close(fig_path);

    fig_overlay = figure('Name', sprintf('Overlay Map%d', map_idx), ...
                         'Position', [50, 50, 900, 800], 'Visible', 'off');
    imagesc(Maps{map_idx});
    colormap(gca, [1 1 1; 0.85 0.85 0.85]);
    axis equal tight;
    hold on;
    overlay_handles = []; overlay_names = {};
    for algo_idx = 1:num_algorithms
        algo_name = Algorithm_Names{algo_idx};
        valid_runs_idx = find(Results(map_idx).(algo_name).Best < collision_threshold);
        if ~isempty(valid_runs_idx)
            [~, best_run_local] = min(Results(map_idx).(algo_name).Best(valid_runs_idx));
            best_run = valid_runs_idx(best_run_local);
            path = Results(map_idx).(algo_name).Paths{best_run};
            if ~isempty(path)
                if strcmp(algo_name, 'EGSO'), lw = 3.0; else, lw = 1.0; end
                h = plot(path(:, 2), path(:, 1), '-', 'Color', [algo_colors(algo_idx, :), 0.7], 'LineWidth', lw);
                overlay_handles = [overlay_handles, h];
                overlay_names{end+1} = algo_name;
            end
        end
    end
    map_size = size(Maps{map_idx}, 1);
    plot(1, 1, 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'LineWidth', 2);
    plot(map_size, map_size, 'rs', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'LineWidth', 2);
    title(sprintf('Map %d (%s) - All Paths Overlay', map_idx, Map_Info(map_idx).size), 'FontSize', 14, 'FontWeight', 'bold');
    if ~isempty(overlay_handles)
        legend(overlay_handles, overlay_names, 'Location', 'best', 'NumColumns', 4, 'FontSize', 7);
    end
    set(gca, 'YDir', 'normal');
    set(gcf, 'Color', 'w');
    hold off;
    print(fig_overlay, sprintf('results/main/PathOverlay_Map%d', map_idx), '-dpng', '-r300');
    close(fig_overlay);
end

%% Statistics table
fid = fopen('results/main/Statistics_Summary.txt', 'w');
fprintf(fid, 'EGSO Path Planning Comparison Experiment Results\n');
fprintf(fid, 'Date: %s\n', datestr(now));
fprintf(fid, 'Config: %d algorithms, 5 maps, %d runs, %d iterations, pop=%d\n', ...
        num_algorithms, num_runs, max_gen, num_population);
fprintf(fid, 'Collision threshold: %.0e\n', collision_threshold);
fprintf(fid, '========================================================================\n\n');

ranking_data = zeros(num_algorithms, 5);

for map_idx = 1:5
    fprintf('\n===== Map %d: %s =====\n', map_idx, Map_Info(map_idx).size);
    fprintf(fid, '===== Map %d: %s =====\n', map_idx, Map_Info(map_idx).size);
    header = sprintf('%-10s %7s %10s %10s %10s %10s %8s', ...
                     'Algorithm', 'SR(%)', 'V.Best', 'V.Mean', 'V.Worst', 'V.Std', 'Time(s)');
    fprintf('%s\n', header); fprintf(fid, '%s\n', header);
    fprintf('%s\n', repmat('-', 1, 75)); fprintf(fid, '%s\n', repmat('-', 1, 75));

    mean_vals = inf(num_algorithms, 1);
    for algo_idx = 1:num_algorithms
        algo_name = Algorithm_Names{algo_idx};
        sr = Results(map_idx).(algo_name).SR;
        valid_fitness = Results(map_idx).(algo_name).Best(Results(map_idx).(algo_name).Best < collision_threshold);
        avg_time = mean(Results(map_idx).(algo_name).Time);
        if ~isempty(valid_fitness)
            v_best = min(valid_fitness); v_mean = Results(map_idx).(algo_name).Mean;
            v_worst = Results(map_idx).(algo_name).Worst; v_std = Results(map_idx).(algo_name).Std;
            mean_vals(algo_idx) = v_mean;
            line_str = sprintf('%-10s %6.1f%% %10.2f %10.2f %10.2f %10.2f %8.2f', ...
                              algo_name, sr, v_best, v_mean, v_worst, v_std, avg_time);
        else
            line_str = sprintf('%-10s %6.1f%% %10s %10s %10s %10s %8.2f', ...
                              algo_name, sr, '---', '---', '---', '---', avg_time);
        end
        fprintf('%s\n', line_str); fprintf(fid, '%s\n', line_str);
    end

    [~, rank_idx] = sort(mean_vals);
    for r = 1:num_algorithms
        ranking_data(rank_idx(r), map_idx) = r;
    end
    [best_mean, best_algo_idx] = min(mean_vals);
    if ~isinf(best_mean)
        best_note = sprintf('  >> Best Valid Mean: %s (%.2f)', Algorithm_Names{best_algo_idx}, best_mean);
        fprintf('%s\n', best_note); fprintf(fid, '%s\n', best_note);
    end
    fprintf(fid, '\n');
end

fprintf('\n\n===== Algorithm Ranking Summary (by Valid Mean) =====\n');
fprintf(fid, '\n===== Algorithm Ranking Summary (by Valid Mean) =====\n');
rank_header = sprintf('%-10s', 'Algorithm');
for m = 1:5, rank_header = [rank_header, sprintf(' %6s', sprintf('Map%d', m))]; end
rank_header = [rank_header, sprintf(' %8s', 'AvgRank')];
fprintf('%s\n', rank_header); fprintf(fid, '%s\n', rank_header);
fprintf('%s\n', repmat('-', 1, 60)); fprintf(fid, '%s\n', repmat('-', 1, 60));

avg_ranks = mean(ranking_data, 2);
[~, overall_rank] = sort(avg_ranks);
for r = 1:num_algorithms
    algo_idx = overall_rank(r);
    algo_name = Algorithm_Names{algo_idx};
    line_str = sprintf('%-10s', algo_name);
    for m = 1:5
        if ranking_data(algo_idx, m) == 0 || isinf(ranking_data(algo_idx, m))
            line_str = [line_str, sprintf(' %6s', '---')];
        else
            line_str = [line_str, sprintf(' %6d', ranking_data(algo_idx, m))];
        end
    end
    line_str = [line_str, sprintf(' %8.1f', avg_ranks(algo_idx))];
    fprintf('%s\n', line_str); fprintf(fid, '%s\n', line_str);
end
fclose(fid);
fprintf('\nSaved: results/main/Statistics_Summary.txt\n');

%% LaTeX table
fid_latex = fopen('results/main/Table_LaTeX.txt', 'w');
fprintf(fid_latex, '%% EGSO Comparison Results - LaTeX Table\n\n');
for map_idx = 1:5
    fprintf(fid_latex, '\\begin{table}[htbp]\n\\centering\n');
    fprintf(fid_latex, '\\caption{Map %d (%s) - Statistical comparison (%d runs)}\n', map_idx, Map_Info(map_idx).size, num_runs);
    fprintf(fid_latex, '\\begin{tabular}{lrrrrrl}\n\\hline\n');
    fprintf(fid_latex, 'Algorithm & SR(\\%%) & V.Best & V.Mean & V.Worst & V.Std & Time(s) \\\\\n\\hline\n');
    for algo_idx = 1:num_algorithms
        algo_name = Algorithm_Names{algo_idx};
        sr = Results(map_idx).(algo_name).SR;
        valid_fitness = Results(map_idx).(algo_name).Best(Results(map_idx).(algo_name).Best < collision_threshold);
        avg_time = mean(Results(map_idx).(algo_name).Time);
        if ~isempty(valid_fitness)
            v_best = min(valid_fitness); v_mean = Results(map_idx).(algo_name).Mean;
            v_worst = Results(map_idx).(algo_name).Worst; v_std = Results(map_idx).(algo_name).Std;
            if strcmp(algo_name, 'EGSO')
                fprintf(fid_latex, '\\textbf{%s} & %.1f & %.2f & \\textbf{%.2f} & %.2f & %.2f & %.2f \\\\\n', ...
                        algo_name, sr, v_best, v_mean, v_worst, v_std, avg_time);
            else
                fprintf(fid_latex, '%s & %.1f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n', ...
                        algo_name, sr, v_best, v_mean, v_worst, v_std, avg_time);
            end
        else
            fprintf(fid_latex, '%s & %.1f & --- & --- & --- & --- & %.2f \\\\\n', algo_name, sr, avg_time);
        end
    end
    fprintf(fid_latex, '\\hline\n\\end{tabular}\n\\end{table}\n\n');
end
fclose(fid_latex);
fprintf('Saved: results/main/Table_LaTeX.txt\n');

fprintf('\n========================================\n');
fprintf('Aggregation finished. All results saved in results/main/\n');
fprintf('========================================\n');
