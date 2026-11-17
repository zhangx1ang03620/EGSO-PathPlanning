% =========================================================================
% Regenerate_Fig8_Heatmap.m
% Success-rate heatmap (23 algorithms x 5 maps), built directly from the stored SR
% field in EGSO_Comparison_Results.mat, using the same 23-algorithm
% ordering as the ranking table (by overall average rank) for visual consistency
% with the rest of the paper, large fonts, and high DPI.
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

% Ranking-table ordering (by overall average rank)
table4_order = {'EGSO','SSA','DE','ARO','GA','GWO','SMA','HHO','SBOA','WOA', ...
                'PSO','INFO','TLBO','ALO','BBO','SCA','MFO','ABC','CS','ACO','FA','FPA','BA'};
assert(numel(table4_order) == num_algorithms, 'algorithm count mismatch');

SR_matrix = zeros(5, num_algorithms);
for map_idx = 1:5
    for k = 1:num_algorithms
        SR_matrix(map_idx, k) = Results(map_idx).(table4_order{k}).SR;
    end
end

fig = figure('Position', [50, 50, 2600, 900], 'Visible', 'off', 'Color', 'w');
imagesc(SR_matrix);
colormap(flipud(autumn(256)));
cb = colorbar;
cb.Label.String = 'Success Rate (%)';
cb.Label.FontSize = 16;
cb.Label.FontWeight = 'bold';
cb.FontSize = 13;
caxis([0 100]);

set(gca, 'XTick', 1:num_algorithms, 'XTickLabel', table4_order, 'XTickLabelRotation', 45, ...
         'YTick', 1:5, 'YTickLabel', arrayfun(@(i) sprintf('Map%d', i), 1:5, 'UniformOutput', false), ...
         'FontSize', 14, 'FontWeight', 'bold', 'TickLength', [0 0]);
xlabel('Algorithm (ordered by overall average rank)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Benchmark Map', 'FontSize', 16, 'FontWeight', 'bold');
title('Success Rate Heat Map: 23 Algorithms x 5 Benchmark Maps', 'FontSize', 19, 'FontWeight', 'bold');

for map_idx = 1:5
    for k = 1:num_algorithms
        v = SR_matrix(map_idx, k);
        txt_color = 'k';
        if v < 50
            txt_color = 'w';
        end
        text(k, map_idx, sprintf('%.0f', v), 'HorizontalAlignment', 'center', ...
             'FontSize', 10, 'Color', txt_color, 'FontWeight', 'bold');
    end
end

set(gcf, 'Color', 'w');
print(fig, 'results/main/Fig8_SuccessRate_Heatmap_hires', '-dpng', '-r450');
close(fig);
fprintf('Saved: results/main/Fig8_SuccessRate_Heatmap_hires.png\n');
fprintf('(success-rate heatmap figure)\n');
