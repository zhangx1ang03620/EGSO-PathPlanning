% =========================================================================
% Run_SamplingBaselines.m
%
% Sampling-based baselines (the RRT* family). Run_ClassicalBaselines.m covers graph
% search (A*) and pure field descent (FMM-descent); this script adds the
% sampling-based paradigm.
%
% RRT* and Informed-RRT* are stochastic (sampling depends on rng), so like the 23
% swarm algorithms they are run num_runs times and reported as SR/V.Best/V.Mean/V.Worst/
% V.Std/Time, not as single deterministic references like A*/FMM-descent.
% =========================================================================
clc; clear; close all;

if ~exist('results/sampling', 'dir')
    mkdir('results/sampling');
end

fprintf('========================================\n');
fprintf('Sampling-based baselines (RRT* / Informed-RRT*)\n');
fprintf('========================================\n\n');

[Maps, Map_Info] = GenerateEGSOMaps();
close(gcf);

num_runs = 20;
max_iter_rrt = 6000;
collision_threshold = 1e4;   % RRT* returns Inf on failure; valid runs are judged by this threshold

Sampling = struct();
for map_idx = 1:5
    G = Maps{map_idx};
    M = size(G,1);
    S = [1,1]; E = [M,M];
    fprintf('===== Map %d (%s) =====\n', map_idx, Map_Info(map_idx).size);

    for variant = 1:2
        is_informed = (variant == 2);
        name = {'RRTStar', 'InformedRRTStar'};
        name = name{variant};

        lens  = zeros(num_runs, 1);
        times = zeros(num_runs, 1);
        paths = cell(num_runs, 1);
        for run = 1:num_runs
            tic;
            [path, len] = RRTStarPath(G, S, E, max_iter_rrt, is_informed);
            times(run) = toc;
            lens(run) = len;
            paths{run} = path;
        end

        valid = isfinite(lens);
        sr = 100 * sum(valid) / num_runs;
        if any(valid)
            v_best  = min(lens(valid));
            v_mean  = mean(lens(valid));
            v_worst = max(lens(valid));
            v_std   = std(lens(valid));
            [~, best_run] = min(lens + (~valid)*1e9);
        else
            v_best = Inf; v_mean = Inf; v_worst = Inf; v_std = 0; best_run = 1;
        end
        avg_time = mean(times);

        Sampling(map_idx).(name).SR      = sr;
        Sampling(map_idx).(name).Best    = v_best;
        Sampling(map_idx).(name).Mean    = v_mean;
        Sampling(map_idx).(name).Worst   = v_worst;
        Sampling(map_idx).(name).Std     = v_std;
        Sampling(map_idx).(name).Time    = avg_time;
        Sampling(map_idx).(name).AllLens = lens;
        Sampling(map_idx).(name).BestPath = paths{best_run};

        fprintf('  %-16s: SR=%6.1f%%  V.Best=%9.2f  V.Mean=%9.2f  V.Worst=%9.2f  V.Std=%7.2f  Time=%.3fs\n', ...
                name, sr, v_best, v_mean, v_worst, v_std, avg_time);
    end
    fprintf('\n');
end

save('results/sampling/Sampling_Baselines.mat', 'Sampling', 'Map_Info', 'num_runs', 'max_iter_rrt');

%% -- Ranking against EGSO / A* / swarm algorithms ----------------------------
fid = fopen('results/sampling/Sampling_vs_EGSO_Summary.txt', 'w');
fprintf(fid, 'Sampling-based References (RRT*, Informed-RRT*) vs EGSO vs A*\n');
fprintf(fid, '========================================================================\n\n');

haveClassical = exist('results/main/Classical_Baselines.mat', 'file');
haveComparison = exist('results/main/EGSO_Comparison_Results.mat', 'file');
if haveClassical, C = load('results/main/Classical_Baselines.mat'); end
if haveComparison, R = load('results/main/EGSO_Comparison_Results.mat'); end

header = sprintf('%-6s %10s %10s %10s %10s %10s %10s', 'Map', 'A*(len)', 'RRT*.Mean', 'RRT*.SR%%', 'IRRT*.Mean', 'IRRT*.SR%%', 'EGSO.Mean');
fprintf('%s\n', header); fprintf(fid, '%s\n', header);
for map_idx = 1:5
    a_len = NaN; egso_mean = NaN;
    if haveClassical, a_len = C.Classical(map_idx).AStar.Length; end
    if haveComparison, egso_mean = R.Results(map_idx).EGSO.Mean; end
    line_str = sprintf('Map%-3d %10.2f %10.2f %9.1f%% %10.2f %9.1f%% %10.2f', map_idx, a_len, ...
        Sampling(map_idx).RRTStar.Mean, Sampling(map_idx).RRTStar.SR, ...
        Sampling(map_idx).InformedRRTStar.Mean, Sampling(map_idx).InformedRRTStar.SR, egso_mean);
    fprintf('%s\n', line_str); fprintf(fid, '%s\n', line_str);
end
fclose(fid);
fprintf('\nSaved: results/sampling/Sampling_vs_EGSO_Summary.txt\n');

%% -- LaTeX table ------------------------------------------------------------
fid = fopen('results/sampling/Sampling_LaTeX.txt', 'w');
fprintf(fid, '%% Table: RRT* / Informed-RRT* sampling-based baselines (%d independent runs)\n', num_runs);
fprintf(fid, '\\begin{table*}[h]\n\\caption{Comparison with sampling-based references (RRT*, Informed-RRT*, %d independent runs) across the five benchmark maps.}\\label{tab:sampling}\n', num_runs);
fprintf(fid, '\\begin{tabular*}{\\textwidth}{@{\\extracolsep\\fill}lccccccc@{}}\n\\toprule\n');
fprintf(fid, 'Map & Algorithm & SR(\\%%) & V.Best & V.Mean & V.Worst & V.Std & Time(s)\\\\\n\\midrule\n');
for map_idx = 1:5
    fprintf(fid, 'Map %d & RRT* & %.1f & %.2f & %.2f & %.2f & %.2f & %.3f\\\\\n', map_idx, ...
        Sampling(map_idx).RRTStar.SR, Sampling(map_idx).RRTStar.Best, Sampling(map_idx).RRTStar.Mean, ...
        Sampling(map_idx).RRTStar.Worst, Sampling(map_idx).RRTStar.Std, Sampling(map_idx).RRTStar.Time);
    fprintf(fid, ' & Informed-RRT* & %.1f & %.2f & %.2f & %.2f & %.2f & %.3f\\\\\n', ...
        Sampling(map_idx).InformedRRTStar.SR, Sampling(map_idx).InformedRRTStar.Best, Sampling(map_idx).InformedRRTStar.Mean, ...
        Sampling(map_idx).InformedRRTStar.Worst, Sampling(map_idx).InformedRRTStar.Std, Sampling(map_idx).InformedRRTStar.Time);
end
fprintf(fid, '\\botrule\n\\end{tabular*}\n\\end{table*}\n');
fclose(fid);
fprintf('Saved: results/sampling/Sampling_LaTeX.txt\n');

%% -- Figure: best EGSO / A* / RRT* / Informed-RRT* paths on the five maps -----
if haveComparison
    for map_idx = 1:5
        G = Maps{map_idx};
        M = size(G,1);
        fig = figure('Position', [100,100,750,700], 'Visible', 'off', 'Color', 'w');
        imagesc(G); colormap(gca, [1 1 1; 0.85 0.85 0.85]); axis equal tight; hold on;

        h = []; names = {};
        if haveClassical
            p = C.Classical(map_idx).AStar.Path;
            if ~isempty(p)
                hh = plot(p(:,2), p(:,1), '-', 'Color', [0.1 0.1 0.1], 'LineWidth', 2.2);
                h=[h hh]; names{end+1} = sprintf('A* (%.1f)', C.Classical(map_idx).AStar.Length);
            end
        end
        p = Sampling(map_idx).RRTStar.BestPath;
        if ~isempty(p)
            hh = plot(p(:,2), p(:,1), '--', 'Color', [0.85 0.55 0.0], 'LineWidth', 1.8);
            h=[h hh]; names{end+1} = sprintf('RRT* (%.1f)', Sampling(map_idx).RRTStar.Best);
        end
        p = Sampling(map_idx).InformedRRTStar.BestPath;
        if ~isempty(p)
            hh = plot(p(:,2), p(:,1), '-.', 'Color', [0.55 0.0 0.65], 'LineWidth', 1.8);
            h=[h hh]; names{end+1} = sprintf('Informed-RRT* (%.1f)', Sampling(map_idx).InformedRRTStar.Best);
        end
        valid_runs_idx = find(R.Results(map_idx).EGSO.Best < collision_threshold);
        if ~isempty(valid_runs_idx)
            [~, best_local] = min(R.Results(map_idx).EGSO.Best(valid_runs_idx));
            best_run = valid_runs_idx(best_local);
            p = R.Results(map_idx).EGSO.Paths{best_run};
            hh = plot(p(:,2), p(:,1), '-', 'Color', [0.90 0.10 0.10], 'LineWidth', 2.5);
            h=[h hh]; names{end+1} = sprintf('EGSO (%.1f)', min(R.Results(map_idx).EGSO.Best(valid_runs_idx)));
        end

        plot(1,1,'go','MarkerSize',10,'MarkerFaceColor','g');
        plot(M,M,'rs','MarkerSize',10,'MarkerFaceColor','r');
        legend(h, names, 'Location', 'best', 'FontSize', 9);
        title(sprintf('Map %d - A* vs RRT* vs Informed-RRT* vs EGSO', map_idx), 'FontSize', 12, 'FontWeight', 'bold');
        set(gca, 'YDir', 'normal'); set(gcf, 'Color', 'w'); hold off;
        print(fig, sprintf('results/sampling/SamplingCompare_Map%d', map_idx), '-dpng', '-r300');
        close(fig);
    end
    fprintf('Saved: results/sampling/SamplingCompare_Map1-5.png\n');
end

fprintf('\n========================================\n');
fprintf('Sampling-based baselines finished\n');
fprintf('========================================\n');
