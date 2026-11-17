% =========================================================================
% Run_SensitivityAnalysis.m
% One-at-a-time sensitivity analysis (set 1), scanned independently on Map 1
% (decoy trap) and Map 5 (large scale), with all other parameters at their defaults.
% See EGSO_Sensitivity.m for the parameters exposed; the values reported in the
% paper are in results/main/SensitivityAnalysis_Summary.txt.
% =========================================================================
clc; clear;

[Maps, Map_Info] = GenerateEGSOMaps();
close(gcf);

map_indices = [1, 5];
max_gen = 100; num_population = 50; num_runs = 10;
collision_threshold = 1e4;

kappa_values      = [0.05, 0.1, 0.2, 0.4];
elite_frac_values = [0.1, 0.2, 0.3, 0.4];
pool_factor_values = [5, 10, 15, 20];

default_kappa = 0.1; default_elite = 0.2; default_pool = 10;

Results = struct();

fprintf('========================================\n');
fprintf('One-at-a-time sensitivity analysis (set 1)\n');
fprintf('========================================\n\n');

for mi = 1:numel(map_indices)
    map_idx = map_indices(mi);
    G = Maps{map_idx};
    M = size(G,1);
    Xmax = M; X_min = 1; dim = Xmax - 2;
    fobj = @(x) fitness(x, G);

    fprintf('===== Map %d (%s) =====\n', map_idx, Map_Info(map_idx).size);

    % --- kappa sweep ---
    fprintf('  -- kappa sweep (elite_frac=%.1f, pool_factor=%d) --\n', default_elite, default_pool);
    kappa_vmean = nan(1,numel(kappa_values)); kappa_std = nan(1,numel(kappa_values));
    for ki = 1:numel(kappa_values)
        fits = zeros(num_runs,1);
        for run = 1:num_runs
            fits(run) = EGSO_Sensitivity(num_population, max_gen, X_min, Xmax, dim, fobj, G, kappa_values(ki), default_elite, default_pool);
        end
        valid = fits < collision_threshold;
        kappa_vmean(ki) = mean(fits(valid)); kappa_std(ki) = std(fits(valid));
        fprintf('    kappa=%.2f: V.Mean=%8.2f +/- %5.2f\n', kappa_values(ki), kappa_vmean(ki), kappa_std(ki));
    end

    % --- elite_frac sweep ---
    fprintf('  -- elite_frac sweep (kappa=%.2f, pool_factor=%d) --\n', default_kappa, default_pool);
    elite_vmean = nan(1,numel(elite_frac_values)); elite_std = nan(1,numel(elite_frac_values));
    for ei = 1:numel(elite_frac_values)
        fits = zeros(num_runs,1);
        for run = 1:num_runs
            fits(run) = EGSO_Sensitivity(num_population, max_gen, X_min, Xmax, dim, fobj, G, default_kappa, elite_frac_values(ei), default_pool);
        end
        valid = fits < collision_threshold;
        elite_vmean(ei) = mean(fits(valid)); elite_std(ei) = std(fits(valid));
        fprintf('    elite_frac=%.1f: V.Mean=%8.2f +/- %5.2f\n', elite_frac_values(ei), elite_vmean(ei), elite_std(ei));
    end

    % --- pool_factor sweep ---
    fprintf('  -- pool_factor sweep (kappa=%.2f, elite_frac=%.1f) --\n', default_kappa, default_elite);
    pool_vmean = nan(1,numel(pool_factor_values)); pool_std = nan(1,numel(pool_factor_values));
    for pi = 1:numel(pool_factor_values)
        fits = zeros(num_runs,1);
        for run = 1:num_runs
            fits(run) = EGSO_Sensitivity(num_population, max_gen, X_min, Xmax, dim, fobj, G, default_kappa, default_elite, pool_factor_values(pi));
        end
        valid = fits < collision_threshold;
        pool_vmean(pi) = mean(fits(valid)); pool_std(pi) = std(fits(valid));
        fprintf('    pool_factor=%d: V.Mean=%8.2f +/- %5.2f\n', pool_factor_values(pi), pool_vmean(pi), pool_std(pi));
    end
    fprintf('\n');

    Results(mi).map_idx = map_idx;
    Results(mi).kappa_vmean = kappa_vmean; Results(mi).kappa_std = kappa_std;
    Results(mi).elite_vmean = elite_vmean; Results(mi).elite_std = elite_std;
    Results(mi).pool_vmean = pool_vmean; Results(mi).pool_std = pool_std;
end

save('results/main/SensitivityAnalysis_Data.mat', 'Results', 'kappa_values', 'elite_frac_values', 'pool_factor_values', 'map_indices', 'Map_Info');

fid = fopen('results/main/SensitivityAnalysis_Summary.txt', 'w');
fprintf(fid, 'One-at-a-time Sensitivity Analysis\n');
fprintf(fid, '=========================================================================\n\n');
for mi = 1:numel(map_indices)
    fprintf(fid, 'Map %d:\n', Results(mi).map_idx);
    fprintf(fid, '  kappa       : '); fprintf(fid, '%8.2f', kappa_values); fprintf(fid, '\n');
    fprintf(fid, '  V.Mean      : '); fprintf(fid, '%8.2f', Results(mi).kappa_vmean); fprintf(fid, '\n\n');
    fprintf(fid, '  elite_frac  : '); fprintf(fid, '%8.1f', elite_frac_values); fprintf(fid, '\n');
    fprintf(fid, '  V.Mean      : '); fprintf(fid, '%8.2f', Results(mi).elite_vmean); fprintf(fid, '\n\n');
    fprintf(fid, '  pool_factor : '); fprintf(fid, '%8d', pool_factor_values); fprintf(fid, '\n');
    fprintf(fid, '  V.Mean      : '); fprintf(fid, '%8.2f', Results(mi).pool_vmean); fprintf(fid, '\n\n');
end
fclose(fid);
fprintf('Saved: results/main/SensitivityAnalysis_Data.mat / Summary.txt\n');
fprintf('\nDone\n');
