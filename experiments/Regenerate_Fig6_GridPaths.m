% =========================================================================
% Regenerate_Fig6_GridPaths.m
% Regenerates the grid-path figure (all 23 algorithms on Map 1 / Map 5) at a
% larger size: larger canvas (4200x3600), larger panel titles (9 -> 13), thicker
% path lines, 450 DPI, and failed algorithms marked with a red "FAIL" title.
%
%
%
%
% =========================================================================
clc; clear;

if ~exist('results/main/EGSO_Comparison_Results.mat', 'file')
    error('Missing results/main/EGSO_Comparison_Results.mat');
end
R = load('results/main/EGSO_Comparison_Results.mat');
Results = R.Results;
Algorithm_Names = R.Algorithm_Names;
Map_Info = R.Map_Info;
num_algorithms = numel(Algorithm_Names);
collision_threshold = 1e4;

[Maps, ~] = GenerateEGSOMaps();
close(gcf);

algo_colors = lines(num_algorithms);

for map_idx = [1, 5]
    fig = figure('Position', [50, 50, 2600, 2200], 'Visible', 'off', 'Color', 'w');
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
                plot(path(:, 2), path(:, 1), '-', 'Color', algo_colors(algo_idx, :), 'LineWidth', 3.5);
                plot(path(1, 2), path(1, 1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
                plot(path(end, 2), path(end, 1), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
            end
            best_val = min(Results(map_idx).(algo_name).Best(valid_runs_idx));
            title_str = sprintf('%s (%.1f)', algo_name, best_val);
            title_color = [0 0 0];
        else
            title_str = sprintf('%s (FAIL)', algo_name);
            title_color = [0.75 0 0];
        end
        title(title_str, 'FontSize', 17, 'FontWeight', 'bold', 'Color', title_color);
        set(gca, 'YDir', 'normal', 'XTick', [], 'YTick', [], 'LineWidth', 1.2);
        hold off;
    end
    sgtitle(sprintf('Map %d (%s): Best Path of Each Algorithm (%d Algorithms)', ...
            map_idx, Map_Info(map_idx).size, num_algorithms), 'FontSize', 22, 'FontWeight', 'bold');
    set(gcf, 'Color', 'w');
    exportgraphics(fig, sprintf('results/main/Fig6_GridPaths_Map%d_hires.png', map_idx), 'Resolution', 400);
    close(fig);
    fprintf('Saved: results/main/Fig6_GridPaths_Map%d_hires.png\n', map_idx);
end

fprintf('\nDone -- the two figures (Map 1, Map 5) form panels (a)/(b)\n');
fprintf('of the grid-path figure\n');
