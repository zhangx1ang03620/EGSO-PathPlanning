% Run_CleanTiming_MainMaps.m
% Clean single-process re-measurement of EGSO's solve time on the five benchmark maps,
% to check whether the Time(s) column of EGSO_Comparison_Results.mat was inflated by
% the five parallel processes; RRT*/Informed-RRT* times are reused for reference.
clc; clear;

[Maps, Map_Info] = GenerateEGSOMaps();
close(gcf);

max_gen = 100; num_population = 50; num_runs = 10;
CleanTime = zeros(5,1);
CleanTime_Std = zeros(5,1);

fprintf('========================================\n');
fprintf('EGSO clean single-process timing (five benchmark maps, %d runs each)\n', num_runs);
fprintf('========================================\n\n');

for map_idx = 1:5
    G = Maps{map_idx};
    M = size(G,1);
    Xmax = M; X_min = 1; dim = Xmax - 2;
    fobj = @(x) fitness(x, G);

    times = zeros(num_runs,1);
    for run = 1:num_runs
        tic;
        EGSO(num_population, max_gen, X_min, Xmax, dim, fobj, G);
        times(run) = toc;
    end
    CleanTime(map_idx) = mean(times);
    CleanTime_Std(map_idx) = std(times);
    fprintf('Map %d (%s): clean solve time = %.3fs (mean of %d), std=%.3f  [original polluted value: see EGSO_Comparison_Results.mat]\n', ...
        map_idx, Map_Info(map_idx).size, CleanTime(map_idx), num_runs, CleanTime_Std(map_idx));
end

if exist('results/main/EGSO_Comparison_Results.mat', 'file')
    CR = load('results/main/EGSO_Comparison_Results.mat');
    fprintf('\n--- Comparison with the Time(s) column of the parallel experiment ---\n');
    for map_idx = 1:5
        orig_mean = mean(CR.Results(map_idx).EGSO.Time);
        fprintf('Map %d: original(parallel)=%.3fs   clean(serial)=%.3fs   ratio=%.2fx\n', ...
            map_idx, orig_mean, CleanTime(map_idx), orig_mean/CleanTime(map_idx));
    end
end

save('results/main/CleanTiming_MainMaps.mat', 'CleanTime', 'CleanTime_Std', 'Map_Info');

if exist('results/main/EVFPrecompute_MainMaps.mat', 'file') && exist('results/sampling/Sampling_Baselines.mat', 'file')
    EV = load('results/main/EVFPrecompute_MainMaps.mat');
    SB = load('results/sampling/Sampling_Baselines.mat');
    fid = fopen('results/main/EndToEnd_Runtime_Table_Clean.txt', 'w');
    fprintf(fid, 'End-to-End Runtime (clean, single-process): EGSO (EVF precompute + solve) vs RRT*/Informed-RRT*\n');
    fprintf(fid, '=========================================================================\n\n');
    header = sprintf('%-8s %-14s %-14s %-16s %-14s %-18s', 'Map', 'EVF-pre(s)', 'EGSO-solve(s)', 'EGSO-total(s)', 'RRT*(s)', 'Informed-RRT*(s)');
    fprintf('%s\n', header); fprintf(fid, '%s\n', header);
    for map_idx = 1:5
        evf_t = EV.EVF_Precompute_Time(map_idx);
        total_t = evf_t + CleanTime(map_idx);
        rrt_t = mean(SB.Sampling(map_idx).RRTStar.Time);
        irrt_t = mean(SB.Sampling(map_idx).InformedRRTStar.Time);
        line = sprintf('%-8d %-14.4f %-14.4f %-16.4f %-14.4f %-18.4f', map_idx, evf_t, CleanTime(map_idx), total_t, rrt_t, irrt_t);
        fprintf('%s\n', line); fprintf(fid, '%s\n', line);
    end
    fclose(fid);
    fprintf('\nSaved: results/main/EndToEnd_Runtime_Table_Clean.txt\n');
end

fprintf('\nDone\n');
