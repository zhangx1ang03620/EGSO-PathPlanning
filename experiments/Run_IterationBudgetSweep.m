% =========================================================================
% Run_IterationBudgetSweep.m
% Iteration-budget sweep: tests whether EGSO's gap to A* on Map 4 and Map 5
% (the two largest gaps, +1.8% and +9.4%) narrows as the iteration budget grows,
% sweeping T in {50,100,200,400} with N = 50 fixed.
%
% =========================================================================
clc; clear;

[Maps, Map_Info] = GenerateEGSOMaps();
close(gcf);

if ~exist('results/main/Classical_Baselines.mat', 'file')
    error('Run_IterationBudgetSweep:Missing', 'results/main/Classical_Baselines.mat (A* reference) is required');
end
CB = load('results/main/Classical_Baselines.mat');

map_indices = [4, 5];
T_values = [50, 100, 200, 400];
num_population = 50;
num_runs = 10;
collision_threshold = 1e4;

GapData = nan(numel(map_indices), numel(T_values));
VMeanData = nan(numel(map_indices), numel(T_values));

fprintf('========================================\n');
fprintf('Iteration-budget sweep: T in %s, Map 4/Map 5\n', mat2str(T_values));
fprintf('========================================\n\n');

for mi = 1:numel(map_indices)
    map_idx = map_indices(mi);
    G = Maps{map_idx};
    M = size(G,1);
    Xmax = M; X_min = 1; dim = Xmax - 2;
    fobj = @(x) fitness(x, G);
    a_star_len = CB.Classical(map_idx).AStar.Length;

    fprintf('===== Map %d (A* length = %.2f) =====\n', map_idx, a_star_len);

    for ti = 1:numel(T_values)
        T = T_values(ti);
        fits = zeros(num_runs,1);
        t0 = tic;
        for run = 1:num_runs
            [Best_score, ~, ~] = EGSO(num_population, T, X_min, Xmax, dim, fobj, G);
            fits(run) = Best_score;
        end
        t_elapsed = toc(t0);
        valid = fits < collision_threshold;
        vmean = mean(fits(valid));
        gap_pct = 100*(vmean - a_star_len)/a_star_len;
        VMeanData(mi, ti) = vmean;
        GapData(mi, ti) = gap_pct;
        fprintf('  T=%4d: V.Mean=%8.2f, Gap-to-A*=%+6.2f%%  (%.1fs for %d runs)\n', T, vmean, gap_pct, t_elapsed, num_runs);
    end
    fprintf('\n');
end

save('results/main/IterationBudgetSweep_Data.mat', 'map_indices', 'T_values', 'GapData', 'VMeanData', 'Map_Info');

fid = fopen('results/main/IterationBudgetSweep_Summary.txt', 'w');
fprintf(fid, 'Iteration-Budget Sweep: EGSO Gap-to-A%% vs T, on Map4/Map5\n');
fprintf(fid, '=========================================================================\n\n');
for mi = 1:numel(map_indices)
    fprintf(fid, 'Map %d:\n', map_indices(mi));
    fprintf(fid, '  T       : '); fprintf(fid, '%8d', T_values); fprintf(fid, '\n');
    fprintf(fid, '  V.Mean  : '); fprintf(fid, '%8.2f', VMeanData(mi,:)); fprintf(fid, '\n');
    fprintf(fid, '  Gap(%%)  : '); fprintf(fid, '%+8.2f', GapData(mi,:)); fprintf(fid, '\n\n');
end
fclose(fid);
fprintf('Saved: results/main/IterationBudgetSweep_Data.mat / Summary.txt\n');
fprintf('\nDone\n');
