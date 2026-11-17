% =========================================================================
% Run_SensitivityAnalysis2.m
% One-at-a-time sensitivity analysis of eta (geodesic-fit strength), the Levy beta
% and the anti-stagnation threshold (stall_limit); together with
% Run_SensitivityAnalysis.m this covers the exploration schedule, the restart
% mechanism and the geodesic-fit parameter.
% =========================================================================
clc; clear;

[Maps, Map_Info] = GenerateEGSOMaps();
close(gcf);

map_indices = [1, 5];
max_gen = 100; num_population = 50; num_runs = 10;
collision_threshold = 1e4;

eta_values        = [0.02, 0.06, 0.10, 0.20];
levy_beta_values  = [1.0, 1.3, 1.5, 2.0];
stall_values      = [3, 5, 8, 12];

default_eta = 0.06; default_beta = 1.5; default_stall = 5;

Results = struct();

fprintf('========================================\n');
fprintf('One-at-a-time sensitivity analysis 2: eta / Levy beta / stall_limit\n');
fprintf('========================================\n\n');

for mi = 1:numel(map_indices)
    map_idx = map_indices(mi);
    G = Maps{map_idx};
    M = size(G,1);
    Xmax = M; X_min = 1; dim = Xmax - 2;
    fobj = @(x) fitness(x, G);

    fprintf('===== Map %d (%s) =====\n', map_idx, Map_Info(map_idx).size);

    fprintf('  -- eta sweep (beta=%.1f, stall=%d) --\n', default_beta, default_stall);
    eta_vmean = nan(1,numel(eta_values)); eta_std = nan(1,numel(eta_values));
    for ki = 1:numel(eta_values)
        fits = zeros(num_runs,1);
        for run = 1:num_runs
            fits(run) = EGSO_Sensitivity2(num_population, max_gen, X_min, Xmax, dim, fobj, G, eta_values(ki), default_beta, default_stall);
        end
        valid = fits < collision_threshold;
        eta_vmean(ki) = mean(fits(valid)); eta_std(ki) = std(fits(valid));
        fprintf('    eta=%.2f: V.Mean=%8.2f +/- %5.2f\n', eta_values(ki), eta_vmean(ki), eta_std(ki));
    end

    fprintf('  -- Levy beta sweep (eta=%.2f, stall=%d) --\n', default_eta, default_stall);
    beta_vmean = nan(1,numel(levy_beta_values)); beta_std = nan(1,numel(levy_beta_values));
    for bi = 1:numel(levy_beta_values)
        fits = zeros(num_runs,1);
        for run = 1:num_runs
            fits(run) = EGSO_Sensitivity2(num_population, max_gen, X_min, Xmax, dim, fobj, G, default_eta, levy_beta_values(bi), default_stall);
        end
        valid = fits < collision_threshold;
        beta_vmean(bi) = mean(fits(valid)); beta_std(bi) = std(fits(valid));
        fprintf('    beta=%.1f: V.Mean=%8.2f +/- %5.2f\n', levy_beta_values(bi), beta_vmean(bi), beta_std(bi));
    end

    fprintf('  -- stall_limit sweep (eta=%.2f, beta=%.1f) --\n', default_eta, default_beta);
    stall_vmean = nan(1,numel(stall_values)); stall_std = nan(1,numel(stall_values));
    for si = 1:numel(stall_values)
        fits = zeros(num_runs,1);
        for run = 1:num_runs
            fits(run) = EGSO_Sensitivity2(num_population, max_gen, X_min, Xmax, dim, fobj, G, default_eta, default_beta, stall_values(si));
        end
        valid = fits < collision_threshold;
        stall_vmean(si) = mean(fits(valid)); stall_std(si) = std(fits(valid));
        fprintf('    stall_limit=%d: V.Mean=%8.2f +/- %5.2f\n', stall_values(si), stall_vmean(si), stall_std(si));
    end
    fprintf('\n');

    Results(mi).map_idx = map_idx;
    Results(mi).eta_vmean = eta_vmean; Results(mi).eta_std = eta_std;
    Results(mi).beta_vmean = beta_vmean; Results(mi).beta_std = beta_std;
    Results(mi).stall_vmean = stall_vmean; Results(mi).stall_std = stall_std;
end

save('results/main/SensitivityAnalysis2_Data.mat', 'Results', 'eta_values', 'levy_beta_values', 'stall_values', 'map_indices', 'Map_Info');

fid = fopen('results/main/SensitivityAnalysis2_Summary.txt', 'w');
fprintf(fid, 'One-at-a-time Sensitivity Analysis 2: eta / Levy beta / stall_limit\n');
fprintf(fid, '=========================================================================\n\n');
for mi = 1:numel(map_indices)
    fprintf(fid, 'Map %d:\n', Results(mi).map_idx);
    fprintf(fid, '  eta          : '); fprintf(fid, '%8.2f', eta_values); fprintf(fid, '\n');
    fprintf(fid, '  V.Mean       : '); fprintf(fid, '%8.2f', Results(mi).eta_vmean); fprintf(fid, '\n\n');
    fprintf(fid, '  Levy beta    : '); fprintf(fid, '%8.1f', levy_beta_values); fprintf(fid, '\n');
    fprintf(fid, '  V.Mean       : '); fprintf(fid, '%8.2f', Results(mi).beta_vmean); fprintf(fid, '\n\n');
    fprintf(fid, '  stall_limit  : '); fprintf(fid, '%8d', stall_values); fprintf(fid, '\n');
    fprintf(fid, '  V.Mean       : '); fprintf(fid, '%8.2f', Results(mi).stall_vmean); fprintf(fid, '\n\n');
end
fclose(fid);
fprintf('Saved: results/main/SensitivityAnalysis2_Data.mat / Summary.txt\n');
fprintf('\nDone\n');
