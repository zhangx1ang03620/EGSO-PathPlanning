% =========================================================================
% Run_ClassicalBaselines.m
% Two deterministic classical references on the five benchmark maps (no repeated runs needed):
%   1. A*  (8-connected, octile heuristic): how far are the swarm methods from a
%      classical graph-search solution?
%   2. FMM-descent (4-connected greedy descent on the same T field EGSO uses, without
%      swarm optimisation or CALS): what does the swarm layer add on top of the
%      guidance field itself?
% Both use exactly the same GenerateSmoothPath post-processing as the other 23
% algorithms, so path lengths are comparable.
% =========================================================================
clc; clear;

if ~exist('results/main', 'dir')
    mkdir('results/main');
end

[Maps, Map_Info] = GenerateEGSOMaps();
close(gcf);

Classical = struct();
for map_idx = 1:5
    G = Maps{map_idx};
    M = size(G, 1);
    S = [1, 1];
    E = [M, M];

    fprintf('===== Map %d (%s) =====\n', map_idx, Map_Info(map_idx).size);

    % ---- A* ----
    tic;
    [pathA, ~] = AStarPath(G, S, E);
    tA = toc;
    if ~isempty(pathA)
        pathA_s = GenerateSmoothPath(pathA, G);
        pathA_s = GenerateSmoothPath(pathA_s, G);
        lenA = 0;
        for i = 1:size(pathA_s,1)-1
            lenA = lenA + norm(pathA_s(i+1,:) - pathA_s(i,:));
        end
    else
        pathA_s = []; lenA = Inf;
    end
    Classical(map_idx).AStar.Path = pathA_s;
    Classical(map_idx).AStar.Length = lenA;
    Classical(map_idx).AStar.Time = tA;
    fprintf('  A*          : length=%.2f, time=%.3fs\n', lenA, tA);

    % ---- FMM-descent ----
    tic;
    [pathF, ~] = FMMDescentPath(G, S, E);
    tF = toc;
    if ~isempty(pathF)
        pathF_s = GenerateSmoothPath(pathF, G);
        pathF_s = GenerateSmoothPath(pathF_s, G);
        lenF = 0;
        for i = 1:size(pathF_s,1)-1
            lenF = lenF + norm(pathF_s(i+1,:) - pathF_s(i,:));
        end
    else
        pathF_s = []; lenF = Inf;
    end
    Classical(map_idx).FMMDescent.Path = pathF_s;
    Classical(map_idx).FMMDescent.Length = lenF;
    Classical(map_idx).FMMDescent.Time = tF;
    fprintf('  FMM-descent : length=%.2f, time=%.3fs\n', lenF, tF);
end

save('results/main/Classical_Baselines.mat', 'Classical', 'Map_Info');
fprintf('\nSaved: results/main/Classical_Baselines.mat\n');

%% Side-by-side table with EGSO and the swarm algorithms
if exist('results/main/EGSO_Comparison_Results.mat', 'file')
    R = load('results/main/EGSO_Comparison_Results.mat');
    fid = fopen('results/main/Classical_vs_EGSO_Summary.txt', 'w');
    fprintf(fid, 'Classical Deterministic References vs EGSO / Best Swarm Baseline\n');
    fprintf(fid, '=========================================================================\n\n');
    header = sprintf('%-8s %10s %10s %14s %14s %10s', 'Map', 'A*(len)', 'FMM-desc', 'EGSO(V.Mean)', 'EGSO(V.Best)', 'EGSO SR%%');
    fprintf('%s\n', header); fprintf(fid, '%s\n', header);
    for map_idx = 1:5
        egso_mean = R.Results(map_idx).EGSO.Mean;
        egso_best = min(R.Results(map_idx).EGSO.Best(R.Results(map_idx).EGSO.Best < 1e4));
        egso_sr   = R.Results(map_idx).EGSO.SR;
        line_str = sprintf('Map%-4d %10.2f %10.2f %14.2f %14.2f %9.1f%%', ...
            map_idx, Classical(map_idx).AStar.Length, Classical(map_idx).FMMDescent.Length, ...
            egso_mean, egso_best, egso_sr);
        fprintf('%s\n', line_str); fprintf(fid, '%s\n', line_str);
    end
    fclose(fid);
    fprintf('\nSaved: results/main/Classical_vs_EGSO_Summary.txt\n');
else
    fprintf('\n[note] results/main/EGSO_Comparison_Results.mat not found; only the classical baselines were saved.\n');
end

%% Path comparison figure (per map: A* vs FMM-descent vs best EGSO)
if exist('results/main/EGSO_Comparison_Results.mat', 'file')
    R = load('results/main/EGSO_Comparison_Results.mat');
    for map_idx = 1:5
        fig = figure('Position', [100, 100, 700, 650], 'Visible', 'off');
        imagesc(Maps{map_idx});
        colormap(gca, [1 1 1; 0.85 0.85 0.85]);
        axis equal tight;
        hold on;

        h = []; names = {};
        if ~isempty(Classical(map_idx).AStar.Path)
            p = Classical(map_idx).AStar.Path;
            hh = plot(p(:,2), p(:,1), '-', 'Color', [0.1 0.1 0.1], 'LineWidth', 2.5);
            h = [h hh]; names{end+1} = sprintf('A* (%.1f)', Classical(map_idx).AStar.Length);
        end
        if ~isempty(Classical(map_idx).FMMDescent.Path)
            p = Classical(map_idx).FMMDescent.Path;
            hh = plot(p(:,2), p(:,1), '--', 'Color', [0.0 0.5 0.8], 'LineWidth', 2);
            h = [h hh]; names{end+1} = sprintf('FMM-descent (%.1f)', Classical(map_idx).FMMDescent.Length);
        end
        valid_runs_idx = find(R.Results(map_idx).EGSO.Best < 1e4);
        if ~isempty(valid_runs_idx)
            [~, best_local] = min(R.Results(map_idx).EGSO.Best(valid_runs_idx));
            best_run = valid_runs_idx(best_local);
            p = R.Results(map_idx).EGSO.Paths{best_run};
            hh = plot(p(:,2), p(:,1), '-', 'Color', [0.90 0.10 0.10], 'LineWidth', 2.5);
            h = [h hh]; names{end+1} = sprintf('EGSO (%.1f)', min(R.Results(map_idx).EGSO.Best(valid_runs_idx)));
        end

        plot(1,1,'go','MarkerSize',10,'MarkerFaceColor','g');
        plot(size(Maps{map_idx},1), size(Maps{map_idx},1), 'rs', 'MarkerSize',10,'MarkerFaceColor','r');
        legend(h, names, 'Location', 'best', 'FontSize', 10);
        title(sprintf('Map %d - A* vs FMM-descent vs EGSO', map_idx), 'FontSize', 13, 'FontWeight', 'bold');
        set(gca, 'YDir', 'normal');
        set(gcf, 'Color', 'w');
        hold off;
        print(fig, sprintf('results/main/ClassicalCompare_Map%d', map_idx), '-dpng', '-r300');
        close(fig);
    end
    fprintf('Saved: results/main/ClassicalCompare_Map1-5.png\n');
end

fprintf('\n========================================\n');
fprintf('Classical baselines finished\n');
fprintf('========================================\n');
