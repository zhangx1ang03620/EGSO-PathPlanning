% =========================================================================
% EGSO main comparison experiment
% EGSO against 22 baseline algorithms on the five benchmark maps
% All results are saved in results/main/ (figures at 300 DPI)
%
% Smoke test: set num_runs=1, max_gen=20 to check that every algorithm runs
% Paper setting: num_runs=20, max_gen=100, num_population=50
% =========================================================================
clc; clear; close all;

%% 0. Output folder
if ~exist('results/main', 'dir')
    mkdir('results/main');
end

%% 1. Benchmark maps
fprintf('========================================\n');
fprintf('Generating the benchmark maps...\n');
fprintf('========================================\n\n');

[Maps, Map_Info] = GenerateEGSOMaps();
print(gcf, 'results/main/EGSO_Test_Maps', '-dpng', '-r300');
close(gcf);
fprintf('Map figure saved: results/main/EGSO_Test_Maps.png (300 DPI)\n');

%% 2. Experiment parameters =================================================
%  Smoke-test parameters (switch to the paper setting after a successful run)
max_gen        = 20;    % paper setting: 100
num_population = 30;    % paper setting: 50
num_runs       = 1;     % paper setting: 20
% ==========================================================================

% ===== Algorithms: EGSO + 19 classical/swarm baselines + 3 recent baselines (23 in total) =====
Algorithm_Names = {'EGSO', 'GWO', 'PSO', 'GA', 'DE', 'SSA', ...
                   'WOA', 'MFO', 'ALO', 'HHO', 'SCA', ...
                   'ABC', 'ACO', 'BBO', 'CS', 'FA', ...
                   'TLBO', 'BA', 'FPA', 'SMA', 'INFO', 'SBOA', 'ARO'};
num_algorithms = length(Algorithm_Names);

% ===== Collision threshold =====
collision_threshold = 1e4;

% ===== Colours (EGSO highlighted in red, the others generated programmatically) =====
algo_colors = [1.00 0.05 0.05; hsv(num_algorithms-1)*0.85];

markers_pool = {'o', 's', '^', 'd', 'v', 'p', 'h', '<', '>', 'x', ...
                '+', '*', 'pentagram', 'hexagram'};
num_markers = length(markers_pool);

fprintf('\nConfiguration:\n');
fprintf('  algorithms: %d\n', num_algorithms);
fprintf('  maps: 5\n');
fprintf('  runs: %d\n', num_runs);
fprintf('  iterations: %d\n', max_gen);
fprintf('  population: %d\n', num_population);
fprintf('  collision threshold: %.0e\n', collision_threshold);

%% 3. Result storage
Results = struct();
for map_idx = 1:5
    Results(map_idx).Map_Name = sprintf('Map%d', map_idx);
    Results(map_idx).Map_Size = Map_Info(map_idx).size;

    for algo_idx = 1:num_algorithms
        algo_name = Algorithm_Names{algo_idx};
        Results(map_idx).(algo_name).Best   = zeros(num_runs, 1);
        Results(map_idx).(algo_name).Mean   = 0;
        Results(map_idx).(algo_name).Worst  = 0;
        Results(map_idx).(algo_name).Std    = 0;
        Results(map_idx).(algo_name).SR     = 0;
        Results(map_idx).(algo_name).Curves = zeros(num_runs, max_gen);
        Results(map_idx).(algo_name).Paths  = cell(num_runs, 1);
        Results(map_idx).(algo_name).Time   = zeros(num_runs, 1);
    end
end

%% 4. Main loop
fprintf('\n========================================\n');
fprintf('Comparison (%d algorithms x %d maps x %d runs)\n', ...
        num_algorithms, 5, num_runs);
fprintf('========================================\n\n');

for map_idx = 1:5
    fprintf('\n========== Map %d/%d: %s ==========\n', ...
            map_idx, 5, Map_Info(map_idx).size);

    G = Maps{map_idx};
    [Xmax, ~] = size(G);
    X_min = 1;
    dimensions = Xmax - 2;
    S = [1, 1];
    E = [Xmax, Xmax];

    fobj = @(x) fitness(x, G);

    for run = 1:num_runs
        fprintf('  run %d/%d ...\n', run, num_runs);

        for algo_idx = 1:num_algorithms
            algo_name = Algorithm_Names{algo_idx};

            try
                tic;
                switch algo_name
                    case 'EGSO'
                        [Best_score, Best_pos, Curve] = EGSO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                    case 'GWO'
                        [Best_score, Best_pos, Curve] = GWO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                    case 'PSO'
                        [Best_score, Best_pos, Curve] = PSO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                    case 'GA'
                        [Best_score, Best_pos, Curve] = GA(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                    case 'DE'
                        [Best_score, Best_pos, Curve] = DE(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                    case 'SSA'
                        [Best_score, Best_pos, Curve] = SSA(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                    case 'WOA'
                        [Best_score, Best_pos, Curve] = WOA(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                    case 'MFO'
                        [Best_score, Best_pos, Curve] = MFO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                    case 'ALO'
                        [Best_score, Best_pos, Curve] = ALO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                    case 'HHO'
                        [Best_score, Best_pos, Curve] = HHO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                    case 'SCA'
                        [Best_score, Best_pos, Curve] = SCA(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                    case 'ABC'
                        [Best_score, Best_pos, Curve] = ABC(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                    case 'ACO'
                        [Best_score, Best_pos, Curve] = ACO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                    case 'BBO'
                        [Best_score, Best_pos, Curve] = BBO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                    case 'CS'
                        [Best_score, Best_pos, Curve] = CS(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                    case 'FA'
                        [Best_score, Best_pos, Curve] = FA(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                    case 'TLBO'
                        [Best_score, Best_pos, Curve] = TLBO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                    case 'BA'
                        [Best_score, Best_pos, Curve] = BA(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                    case 'FPA'
                        [Best_score, Best_pos, Curve] = FPA(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                    case 'SMA'
                        [Best_score, Best_pos, Curve] = SMA(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                    case 'INFO'
                        [Best_score, Best_pos, Curve] = INFO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                    case 'SBOA'
                        [Best_score, Best_pos, Curve] = SBOA(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                    case 'ARO'
                        [Best_score, Best_pos, Curve] = ARO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G);
                    otherwise
                        Best_score = inf;
                        Best_pos = zeros(1, dimensions);
                        Curve = inf(1, max_gen);
                end
                elapsed_time = toc;

                % Generate the path
                Best_pos = round(Best_pos);
                Best_pos = max(Best_pos, 1);
                Best_pos = min(Best_pos, Xmax);
                route = [S(1) Best_pos E(1)];
                path = generateContinuousRoute(route, G);
                if ~isempty(path)
                    path = GenerateSmoothPath(path, G);
                    path = GenerateSmoothPath(path, G);
                end

                % Store the results
                Results(map_idx).(algo_name).Best(run)   = Best_score;
                Results(map_idx).(algo_name).Curves(run, :) = Curve;
                Results(map_idx).(algo_name).Paths{run}  = path;
                Results(map_idx).(algo_name).Time(run)   = elapsed_time;

                if Best_score < collision_threshold
                    status = 'collision-free';
                else
                    status = 'collision';
                end
                fprintf('    [%s] %s done: fitness=%.4f, time=%.2fs\n', ...
                        status, algo_name, Best_score, elapsed_time);

            catch ME
                fprintf('    [error] %s failed: %s\n', algo_name, ME.message);
                Results(map_idx).(algo_name).Best(run) = inf;
            end
        end
    end

    %% Statistics (collision-free valid paths only)
    for algo_idx = 1:num_algorithms
        algo_name = Algorithm_Names{algo_idx};

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

    %% ===== Convergence curves (300 DPI) =====
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
    legend(plot_handles, legend_entries, 'Location', 'northeast', ...
           'NumColumns', 4, 'FontSize', 7);
    grid on;
    set(gca, 'YScale', 'log');
    set(gcf, 'Color', 'w');
    hold off;

    print(fig_conv, sprintf('results/main/Convergence_Map%d', map_idx), '-dpng', '-r300');
    close(fig_conv);

    %% ===== Path grid (5x5, 300 DPI) =====
    fig_path = figure('Name', sprintf('Path Map%d', map_idx), ...
                      'Position', [50, 50, 1800, 1500], 'Visible', 'off');

    for algo_idx = 1:num_algorithms
        algo_name = Algorithm_Names{algo_idx};
        subplot(5, 5, algo_idx);
        imagesc(G);
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
            map_idx, Map_Info(map_idx).size, num_algorithms), ...
            'FontSize', 14, 'FontWeight', 'bold');
    set(gcf, 'Color', 'w');

    print(fig_path, sprintf('results/main/Paths_Map%d', map_idx), '-dpng', '-r300');
    close(fig_path);

    %% ===== Overlay of the best paths =====
    fig_overlay = figure('Name', sprintf('Overlay Map%d', map_idx), ...
                         'Position', [50, 50, 900, 800], 'Visible', 'off');
    imagesc(G);
    colormap(gca, [1 1 1; 0.85 0.85 0.85]);
    axis equal tight;
    hold on;

    overlay_handles = [];
    overlay_names = {};

    for algo_idx = 1:num_algorithms
        algo_name = Algorithm_Names{algo_idx};
        valid_runs_idx = find(Results(map_idx).(algo_name).Best < collision_threshold);

        if ~isempty(valid_runs_idx)
            [~, best_run_local] = min(Results(map_idx).(algo_name).Best(valid_runs_idx));
            best_run = valid_runs_idx(best_run_local);
            path = Results(map_idx).(algo_name).Paths{best_run};

            if ~isempty(path)
                if strcmp(algo_name, 'EGSO')
                    lw = 3.0;
                else
                    lw = 1.0;
                end
                h = plot(path(:, 2), path(:, 1), '-', ...
                         'Color', [algo_colors(algo_idx, :), 0.7], ...
                         'LineWidth', lw);
                overlay_handles = [overlay_handles, h];
                overlay_names{end+1} = algo_name;
            end
        end
    end

    map_size = size(G, 1);
    plot(1, 1, 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'LineWidth', 2);
    plot(map_size, map_size, 'rs', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'LineWidth', 2);

    title(sprintf('Map %d (%s) - All Paths Overlay', map_idx, Map_Info(map_idx).size), ...
          'FontSize', 14, 'FontWeight', 'bold');
    if ~isempty(overlay_handles)
        legend(overlay_handles, overlay_names, 'Location', 'best', ...
               'NumColumns', 4, 'FontSize', 7);
    end
    set(gca, 'YDir', 'normal');
    set(gcf, 'Color', 'w');
    hold off;

    print(fig_overlay, sprintf('results/main/PathOverlay_Map%d', map_idx), '-dpng', '-r300');
    close(fig_overlay);
end

%% 5. Save the experiment data
save('results/main/EGSO_Comparison_Results.mat', 'Results', 'Map_Info', ...
     'Algorithm_Names', 'num_runs', 'max_gen', 'num_population');
fprintf('\nSaved: results/main/EGSO_Comparison_Results.mat\n');

%% 6. Summary statistics table (paper format)
fprintf('\n========================================\n');
fprintf('Summary statistics table (paper format)\n');
fprintf('========================================\n');

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
    fprintf('%s\n', header);
    fprintf(fid, '%s\n', header);

    sep_line = repmat('-', 1, 75);
    fprintf('%s\n', sep_line);
    fprintf(fid, '%s\n', sep_line);

    mean_vals = inf(num_algorithms, 1);

    for algo_idx = 1:num_algorithms
        algo_name = Algorithm_Names{algo_idx};
        sr = Results(map_idx).(algo_name).SR;

        valid_fitness = Results(map_idx).(algo_name).Best(Results(map_idx).(algo_name).Best < collision_threshold);

        if ~isempty(valid_fitness)
            v_best  = min(valid_fitness);
            v_mean  = Results(map_idx).(algo_name).Mean;
            v_worst = Results(map_idx).(algo_name).Worst;
            v_std   = Results(map_idx).(algo_name).Std;
            mean_vals(algo_idx) = v_mean;
        else
            v_best = NaN; v_mean = NaN; v_worst = NaN; v_std = NaN;
        end

        avg_time = mean(Results(map_idx).(algo_name).Time);

        if isnan(v_best)
            line_str = sprintf('%-10s %6.1f%% %10s %10s %10s %10s %8.2f', ...
                              algo_name, sr, '---', '---', '---', '---', avg_time);
        else
            line_str = sprintf('%-10s %6.1f%% %10.2f %10.2f %10.2f %10.2f %8.2f', ...
                              algo_name, sr, v_best, v_mean, v_worst, v_std, avg_time);
        end
        fprintf('%s\n', line_str);
        fprintf(fid, '%s\n', line_str);
    end

    [~, rank_idx] = sort(mean_vals);
    for r = 1:num_algorithms
        ranking_data(rank_idx(r), map_idx) = r;
    end

    [best_mean, best_algo_idx] = min(mean_vals);
    if ~isinf(best_mean)
        best_note = sprintf('  >> Best Valid Mean: %s (%.2f)', Algorithm_Names{best_algo_idx}, best_mean);
        fprintf('%s\n', best_note);
        fprintf(fid, '%s\n', best_note);
    end

    fprintf(fid, '\n');
end

fprintf('\n\n===== Algorithm Ranking Summary (by Valid Mean) =====\n');
fprintf(fid, '\n===== Algorithm Ranking Summary (by Valid Mean) =====\n');

rank_header = sprintf('%-10s', 'Algorithm');
for m = 1:5
    rank_header = [rank_header, sprintf(' %6s', sprintf('Map%d', m))];
end
rank_header = [rank_header, sprintf(' %8s', 'AvgRank')];
fprintf('%s\n', rank_header);
fprintf(fid, '%s\n', rank_header);
fprintf('%s\n', repmat('-', 1, 60));
fprintf(fid, '%s\n', repmat('-', 1, 60));

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
    fprintf('%s\n', line_str);
    fprintf(fid, '%s\n', line_str);
end

fclose(fid);
fprintf('\nSaved: results/main/Statistics_Summary.txt\n');

%% 7. LaTeX table
fid_latex = fopen('results/main/Table_LaTeX.txt', 'w');
fprintf(fid_latex, '%% EGSO Comparison Results - LaTeX Table\n\n');

for map_idx = 1:5
    fprintf(fid_latex, '\\begin{table}[htbp]\n');
    fprintf(fid_latex, '\\centering\n');
    fprintf(fid_latex, '\\caption{Map %d (%s) - Statistical comparison (%d runs)}\n', ...
            map_idx, Map_Info(map_idx).size, num_runs);
    fprintf(fid_latex, '\\begin{tabular}{lrrrrrl}\n');
    fprintf(fid_latex, '\\hline\n');
    fprintf(fid_latex, 'Algorithm & SR(\\%%) & V.Best & V.Mean & V.Worst & V.Std & Time(s) \\\\\n');
    fprintf(fid_latex, '\\hline\n');

    for algo_idx = 1:num_algorithms
        algo_name = Algorithm_Names{algo_idx};
        sr = Results(map_idx).(algo_name).SR;
        valid_fitness = Results(map_idx).(algo_name).Best(Results(map_idx).(algo_name).Best < collision_threshold);
        avg_time = mean(Results(map_idx).(algo_name).Time);

        if ~isempty(valid_fitness)
            v_best  = min(valid_fitness);
            v_mean  = Results(map_idx).(algo_name).Mean;
            v_worst = Results(map_idx).(algo_name).Worst;
            v_std   = Results(map_idx).(algo_name).Std;

            if strcmp(algo_name, 'EGSO')
                fprintf(fid_latex, '\\textbf{%s} & %.1f & %.2f & \\textbf{%.2f} & %.2f & %.2f & %.2f \\\\\n', ...
                        algo_name, sr, v_best, v_mean, v_worst, v_std, avg_time);
            else
                fprintf(fid_latex, '%s & %.1f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n', ...
                        algo_name, sr, v_best, v_mean, v_worst, v_std, avg_time);
            end
        else
            fprintf(fid_latex, '%s & %.1f & --- & --- & --- & --- & %.2f \\\\\n', ...
                    algo_name, sr, avg_time);
        end
    end

    fprintf(fid_latex, '\\hline\n');
    fprintf(fid_latex, '\\end{tabular}\n');
    fprintf(fid_latex, '\\end{table}\n\n');
end

fclose(fid_latex);
fprintf('Saved: results/main/Table_LaTeX.txt\n');

fprintf('\n========================================\n');
fprintf('Finished. All results saved in results/main/\n');
fprintf('========================================\n');
