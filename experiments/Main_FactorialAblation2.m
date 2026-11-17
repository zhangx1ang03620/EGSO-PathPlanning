% =========================================================================
% Main_FactorialAblation2.m
% Full crossed factorial ablation: 3 (initialisation weighting field) x 3 (elite guidance field)
% + 3 (random initialisation, per elite field) = 12 cells, with initialisation and elite
% guidance varied independently (the 6-cell Main_FactorialAblation.m shares one field for both).
%
% Map 1/3/5, N=50, T=100, 10 independent runs, same protocol as the first factorial
% ablation and the ablation table, so the results can be cross-checked.
% =========================================================================
clc; clear; close all;

if ~exist('results/factorial2', 'dir')
    mkdir('results/factorial2');
end

fprintf('========================================\n');
fprintf('EGSO crossed factorial ablation (init_field x elite_field, 12 cells)\n');
fprintf('========================================\n\n');

[Maps, Map_Info] = GenerateEGSOMaps();
close(gcf);

map_indices = [1, 3, 5];
max_gen        = 100;
num_population = 50;
num_runs       = 10;
collision_threshold = 1e4;

fields = {'eikonal', 'apf', 'dijkstra'};

% 12 cells: {init_mode, init_field, elite_field, label}
Cells = {};
for i = 1:3
    for j = 1:3
        Cells{end+1} = {'gsi', fields{i}, fields{j}, sprintf('GSI(%s)+Elite(%s)', fields{i}, fields{j})}; %#ok<SAGROW>
    end
end
for j = 1:3
    Cells{end+1} = {'random', 'eikonal', fields{j}, sprintf('Random+Elite(%s)', fields{j})}; %#ok<SAGROW>
end
num_cells = numel(Cells);

SR_data   = zeros(num_cells, numel(map_indices));
Mean_data = nan(num_cells, numel(map_indices));
Std_data  = nan(num_cells, numel(map_indices));

t_start_all = tic;
for mi = 1:numel(map_indices)
    map_idx = map_indices(mi);
    G = Maps{map_idx};
    M = size(G, 1);
    Xmax = M; X_min = 1; dim = Xmax - 2;
    fobj = @(x) fitness(x, G);

    fprintf('\n===== Map %d (%s) =====\n', map_idx, Map_Info(map_idx).size);

    for ci = 1:num_cells
        init_mode   = Cells{ci}{1};
        init_field  = Cells{ci}{2};
        elite_field = Cells{ci}{3};
        label       = Cells{ci}{4};

        fits = zeros(num_runs, 1);
        t0 = tic;
        for run = 1:num_runs
            fits(run) = EGSO_Variant2(num_population, max_gen, X_min, Xmax, dim, fobj, G, init_mode, init_field, elite_field);
        end
        t_elapsed = toc(t0);

        valid = fits < collision_threshold;
        SR_data(ci, mi) = sum(valid) / num_runs * 100;
        if any(valid)
            Mean_data(ci, mi) = mean(fits(valid));
            Std_data(ci, mi)  = std(fits(valid));
        end
        fprintf('  [%-28s] SR=%6.1f%%, V.Mean=%8.2f, V.Std=%6.2f  (%.1fs)\n', ...
            label, SR_data(ci, mi), Mean_data(ci, mi), Std_data(ci, mi), t_elapsed);
    end
end
fprintf('\nTotal time: %.1f min\n', toc(t_start_all)/60);

save('results/factorial2/FactorialAblation2_Data.mat', 'Cells', 'map_indices', 'SR_data', 'Mean_data', 'Std_data', ...
     'num_runs', 'max_gen', 'num_population', 'Map_Info');

%% Text summary
fid = fopen('results/factorial2/FactorialAblation2_Summary.txt', 'w');
fprintf(fid, 'Fully-Crossed Factorial Ablation: 3(init field) x 3(elite field) + 3(random init), Map1/3/5\n');
fprintf(fid, '=========================================================================\n\n');
for ci = 1:num_cells
    fprintf(fid, '%-30s', Cells{ci}{4});
    for mi = 1:numel(map_indices)
        fprintf(fid, ' | Map%d SR=%6.1f%% V.Mean=%8.2f', map_indices(mi), SR_data(ci,mi), Mean_data(ci,mi));
    end
    fprintf(fid, '\n');
end

%% Main effects: elite field (GSI initialisation, averaged over init_field) vs. initialisation
fprintf(fid, '\n--- Main effect: elite_field (averaged over the 3 GSI init_field choices), Map avg over 1/3/5 ---\n');
for j = 1:3
    idxs = find(cellfun(@(c) strcmp(c{1},'gsi') && strcmp(c{3}, fields{j}), Cells));
    vals = Mean_data(idxs, :);
    fprintf(fid, '  elite=%-8s: overall mean V.Mean = %.2f\n', fields{j}, mean(vals(:), 'omitnan'));
end
fprintf(fid, '\n--- Main effect: init_field (averaged over the 3 elite_field choices), GSI only ---\n');
for i = 1:3
    idxs = find(cellfun(@(c) strcmp(c{1},'gsi') && strcmp(c{2}, fields{i}), Cells));
    vals = Mean_data(idxs, :);
    fprintf(fid, '  init=%-8s: overall mean V.Mean = %.2f\n', fields{i}, mean(vals(:), 'omitnan'));
end
fprintf(fid, '\n--- Main effect: init_mode (gsi vs random), averaged over elite_field, using elite=matching field for gsi diag---\n');
gsi_idxs = find(cellfun(@(c) strcmp(c{1},'gsi'), Cells));
rand_idxs = find(cellfun(@(c) strcmp(c{1},'random'), Cells));
fprintf(fid, '  gsi (all 9 cells)   : overall mean V.Mean = %.2f\n', mean(Mean_data(gsi_idxs,:), 'all', 'omitnan'));
fprintf(fid, '  random (3 cells)    : overall mean V.Mean = %.2f\n', mean(Mean_data(rand_idxs,:), 'all', 'omitnan'));
fclose(fid);
fprintf('Saved: results/factorial2/FactorialAblation2_Data.mat / Summary.txt\n');

%% Internal consistency: the three diagonal cells should be close to the earlier ablations
fprintf('\n===== Internal consistency check (diagonal cells vs. ablation table / first factorial ablation) =====\n');
diag_labels = {'GSI(eikonal)+Elite(eikonal)', 'GSI(apf)+Elite(apf)', 'GSI(dijkstra)+Elite(dijkstra)'};
for d = 1:numel(diag_labels)
    idx = find(strcmp(cellfun(@(c) c{4}, Cells, 'UniformOutput', false), diag_labels{d}));
    fprintf('  %s: Map1=%.2f Map3=%.2f Map5=%.2f\n', diag_labels{d}, Mean_data(idx,1), Mean_data(idx,2), Mean_data(idx,3));
end

fprintf('\n========================================\nCrossed factorial ablation finished\n========================================\n');
