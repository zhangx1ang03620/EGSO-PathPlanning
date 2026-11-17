% =========================================================================
% Main_FactorialAblation.m
% Full 2 (initialisation) x 3 (guidance field) factorial ablation: varies
% (a) the initialisation (GSI weighting vs. uniform random feasible sampling) and
% (b) the elite-layer guidance field (Eikonal / APF / Dijkstra-BFS) in one pipeline (EGSO_Variant.m).
%
%   (GSI,    Eikonal)  == Full EGSO
%   (GSI,    APF)      == the w/o EVF ablation
%   (GSI,    Dijkstra) == new cell
%   (Random, Eikonal)  == new cell
%   (Random, APF)      == new cell
%   (Random, Dijkstra) == new cell
%
% Maps and parameters as in Main_Ablation_Experiment.m: Map 1/3/5, N=50,
% T=100, 10 independent runs, so the GSI+Eikonal and GSI+APF cells can be
% checked against the Full EGSO / w/o EVF rows of the ablation table
% (internal consistency check).
% =========================================================================
clc; clear; close all;

if ~exist('results/factorial', 'dir')
    mkdir('results/factorial');
end

fprintf('========================================\n');
fprintf('EGSO factorial ablation (2 init x 3 field)\n');
fprintf('========================================\n\n');

[Maps, Map_Info] = GenerateEGSOMaps();
close(gcf);

map_indices = [1, 3, 5];
max_gen        = 100;
num_population = 50;
num_runs       = 10;
collision_threshold = 1e4;

Combos = {
    'gsi',    'eikonal',  'GSI + Eikonal (= Full EGSO)';
    'gsi',    'apf',      'GSI + APF (= w/o EVF)';
    'gsi',    'dijkstra', 'GSI + Dijkstra/BFS';
    'random', 'eikonal',  'Random + Eikonal';
    'random', 'apf',      'Random + APF';
    'random', 'dijkstra', 'Random + Dijkstra/BFS';
};
num_combos = size(Combos, 1);

SR_data   = zeros(num_combos, numel(map_indices));
Mean_data = nan(num_combos, numel(map_indices));
Std_data  = nan(num_combos, numel(map_indices));
AllFits   = cell(num_combos, numel(map_indices));

t_start_all = tic;
for mi = 1:numel(map_indices)
    map_idx = map_indices(mi);
    G = Maps{map_idx};
    M = size(G, 1);
    Xmax = M; X_min = 1; dim = Xmax - 2;
    fobj = @(x) fitness(x, G);

    fprintf('\n===== Map %d (%s) =====\n', map_idx, Map_Info(map_idx).size);

    for ci = 1:num_combos
        init_mode  = Combos{ci, 1};
        field_mode = Combos{ci, 2};
        label      = Combos{ci, 3};

        fits = zeros(num_runs, 1);
        t0 = tic;
        for run = 1:num_runs
            [Best_score, ~, ~] = EGSO_Variant(num_population, max_gen, X_min, Xmax, dim, fobj, G, init_mode, field_mode);
            fits(run) = Best_score;
        end
        t_elapsed = toc(t0);

        AllFits{ci, mi} = fits;
        valid = fits < collision_threshold;
        SR_data(ci, mi) = sum(valid) / num_runs * 100;
        if any(valid)
            Mean_data(ci, mi) = mean(fits(valid));
            Std_data(ci, mi)  = std(fits(valid));
        end
        fprintf('  [%-24s] SR=%6.1f%%, V.Mean=%8.2f, V.Std=%6.2f  (%.1fs for %d runs)\n', ...
            label, SR_data(ci, mi), Mean_data(ci, mi), Std_data(ci, mi), t_elapsed, num_runs);
    end
end
fprintf('\nTotal time: %.1f min\n', toc(t_start_all)/60);

save('results/factorial/FactorialAblation_Data.mat', 'Combos', 'map_indices', 'SR_data', 'Mean_data', 'Std_data', ...
     'AllFits', 'num_runs', 'max_gen', 'num_population', 'Map_Info');

%% Text summary
fid = fopen('results/factorial/FactorialAblation_Summary.txt', 'w');
fprintf(fid, 'Factorial Ablation: 2 (init) x 3 (guidance field), Map1/3/5, %d independent runs each\n', num_runs);
fprintf(fid, '=========================================================================\n\n');
fprintf(fid, '%-32s', 'Combination');
for mi = 1:numel(map_indices)
    fprintf(fid, ' | Map%d SR(%%)  V.Mean   V.Std ', map_indices(mi));
end
fprintf(fid, '\n');
for ci = 1:num_combos
    fprintf(fid, '%-32s', Combos{ci,3});
    for mi = 1:numel(map_indices)
        fprintf(fid, ' | %8.1f  %8.2f  %6.2f', SR_data(ci,mi), Mean_data(ci,mi), Std_data(ci,mi));
    end
    fprintf(fid, '\n');
end
fclose(fid);
fprintf('Saved: results/factorial/FactorialAblation_Data.mat / Summary.txt\n');

%% Consistency with the ablation table (GSI+Eikonal ~ Full EGSO, GSI+APF ~ w/o EVF)
if exist('results/ablation/Ablation_Data.mat', 'file')
    A = load('results/ablation/Ablation_Data.mat');
    fprintf('\n===== Internal consistency check (vs. ablation table) =====\n');
    fprintf('(The ablation-table rows come from EGSO.m/EGSO_noEVF.m; GSI+Eikonal/GSI+APF here\n');
    fprintf(' come from EGSO_Variant.m with the same mechanisms, so values should agree within run-to-run variation)\n\n');
    disp(A);
else
    fprintf('\n[note] results/ablation/Ablation_Data.mat not found; skipping the consistency check\n');
end

fprintf('\n========================================\n');
fprintf('Factorial ablation finished\n');
fprintf('========================================\n');
