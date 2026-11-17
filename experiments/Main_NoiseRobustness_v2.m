% =========================================================================
% Main_NoiseRobustness_v2.m
% Execution-noise robustness experiment, averaged over runs.
%
% Main_NoiseRobustness.m and Main_ClearanceTunable_NoiseTest.m run the Monte Carlo
% only on the single best-fitness path of 20/10 runs, but clearance and noise tolerance
% are sensitive to the particular geometry of that path and do not represent a typical
% run: on Map 4 the sigma=0.1 survival of the "best" path changed from 27.7% to 44.7%
% between two independent sets of 10 runs, only because a different best path was found.
% A single path is therefore an unreliable estimate.
%
% Here EGSO (baseline: min_clearance=1.0) and EGSO (elevated: 2.0) are each run
% 10 times, a Monte Carlo is run on every path, and the survival curves are averaged
% (the mean represents a typical run of the setting, not the luckiest one). A* is
% deterministic, has a single path and is not affected.
% =========================================================================
clc; clear;

if ~exist('results/noise', 'dir')
    mkdir('results/noise');
end

[Maps, Map_Info] = GenerateEGSOMaps();
close(gcf);

if ~exist('results/main/Classical_Baselines.mat', 'file')
    error('Main_NoiseRobustness_v2:Missing', 'results/main/Classical_Baselines.mat (A* paths) is required');
end
CB = load('results/main/Classical_Baselines.mat');

max_gen        = 100;
num_population = 50;
num_runs       = 10;
collision_threshold = 1e4;
clearance_settings = [1.0, 2.0];
sigma_levels = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0];
num_trials = 200;   % Monte Carlo trials per path and sigma (10 runs keep the aggregate stable)

MeanClear_AStar = nan(5,1); MinClear_AStar = nan(5,1);
Survival_AStar  = nan(5, numel(sigma_levels));

MeanClear_EGSO = nan(5, numel(clearance_settings));   % mean over 10 runs
MinClear_EGSO  = nan(5, numel(clearance_settings));
Survival_EGSO_mean = nan(5, numel(clearance_settings), numel(sigma_levels));
Survival_EGSO_std  = nan(5, numel(clearance_settings), numel(sigma_levels));
VMean_EGSO = nan(5, numel(clearance_settings));

fprintf('========================================\n');
fprintf('Execution-noise robustness v2 (mean over 10 runs, not the single best path)\n');
fprintf('========================================\n\n');

t_start = tic;
for map_idx = 1:5
    G = Maps{map_idx};
    M = size(G, 1);
    Xmax = M; X_min = 1; dim = Xmax - 2;
    fobj = @(x) fitness(x, G);
    obstacleMask = (G == 1);
    T_clear = FastMarchingField(G, obstacleMask, 'clearance');

    fprintf('===== Map %d (%s) =====\n', map_idx, Map_Info(map_idx).size);

    %% ---- A* (deterministic, single path) ----
    pathA = CB.Classical(map_idx).AStar.Path;
    [MeanClear_AStar(map_idx), MinClear_AStar(map_idx)] = pathClearanceStats(pathA, T_clear, M);
    for sgi = 1:numel(sigma_levels)
        Survival_AStar(map_idx, sgi) = monteCarloSurvival(pathA, G, sigma_levels(sgi), num_trials*2); % A* has a single path; more trials to reduce variance
    end
    fprintf('  A* (deterministic): MeanClr=%.3f MinClr=%.3f\n', MeanClear_AStar(map_idx), MinClear_AStar(map_idx));

    %% ---- EGSO at each clearance setting, averaged over num_runs ----
    for si = 1:numel(clearance_settings)
        mc = clearance_settings(si);
        run_survival = nan(num_runs, numel(sigma_levels));
        run_meanclr  = nan(num_runs, 1);
        run_minclr   = nan(num_runs, 1);
        run_fit      = nan(num_runs, 1);

        t0 = tic;
        for run = 1:num_runs
            [Best_score, ~, ~, Best_path] = EGSO_ClearanceTunable(num_population, max_gen, X_min, Xmax, dim, fobj, G, mc);
            run_fit(run) = Best_score;
            if Best_score < collision_threshold && ~isempty(Best_path)
                [run_meanclr(run), run_minclr(run)] = pathClearanceStats(Best_path, T_clear, M);
                for sgi = 1:numel(sigma_levels)
                    run_survival(run, sgi) = monteCarloSurvival(Best_path, G, sigma_levels(sgi), num_trials);
                end
            end
        end
        t_elapsed = toc(t0);

        VMean_EGSO(map_idx, si) = mean(run_fit(run_fit < collision_threshold));
        MeanClear_EGSO(map_idx, si) = mean(run_meanclr, 'omitnan');
        MinClear_EGSO(map_idx, si)  = mean(run_minclr, 'omitnan');
        Survival_EGSO_mean(map_idx, si, :) = mean(run_survival, 1, 'omitnan');
        Survival_EGSO_std(map_idx, si, :)  = std(run_survival, 0, 1, 'omitnan');

        fprintf('  [min_clearance=%.1f, %d runs] V.Mean=%.2f, MeanClr=%.3f(avg), MinClr=%.3f(avg)  (%.1fs)\n', ...
            mc, num_runs, VMean_EGSO(map_idx,si), MeanClear_EGSO(map_idx,si), MinClear_EGSO(map_idx,si), t_elapsed);
        fprintf('    Survival(%%) mean over %d runs @ sigma=%s:\n      ', num_runs, mat2str(sigma_levels));
        fprintf('%7.1f', squeeze(Survival_EGSO_mean(map_idx, si, :)));
        fprintf('\n      (std: '); fprintf('%6.1f', squeeze(Survival_EGSO_std(map_idx, si, :))); fprintf(')\n');
    end
    fprintf('  A* Survival(%%) @ sigma: '); fprintf('%7.1f', Survival_AStar(map_idx,:)); fprintf('\n\n');
end
fprintf('Total time: %.1f min\n', toc(t_start)/60);

save('results/noise/NoiseRobustness_v2_Data.mat', 'clearance_settings', 'sigma_levels', ...
     'MeanClear_AStar', 'MinClear_AStar', 'Survival_AStar', ...
     'MeanClear_EGSO', 'MinClear_EGSO', 'Survival_EGSO_mean', 'Survival_EGSO_std', 'VMean_EGSO', ...
     'Map_Info', 'num_runs', 'num_trials');

%% Text summary
fid = fopen('results/noise/NoiseRobustness_v2_Summary.txt', 'w');
fprintf(fid, 'Execution-Noise Robustness v2: A* vs EGSO (baseline min_clearance=1.0, elevated=2.0)\n');
fprintf(fid, 'EGSO numbers are MEAN over %d independent runs (not a single best-fitness path)\n', num_runs);
fprintf(fid, '=========================================================================\n\n');
for map_idx = 1:5
    fprintf(fid, 'Map %d (%s):\n', map_idx, Map_Info(map_idx).size);
    fprintf(fid, '  A*:  MeanClr=%.3f MinClr=%.3f\n', MeanClear_AStar(map_idx), MinClear_AStar(map_idx));
    fprintf(fid, '       Survival(%%): '); fprintf(fid, '%7.1f', Survival_AStar(map_idx,:)); fprintf(fid, '\n');
    for si = 1:numel(clearance_settings)
        fprintf(fid, '  EGSO(min_clearance=%.1f): V.Mean=%.2f MeanClr=%.3f MinClr=%.3f\n', ...
            clearance_settings(si), VMean_EGSO(map_idx,si), MeanClear_EGSO(map_idx,si), MinClear_EGSO(map_idx,si));
        fprintf(fid, '       Survival(%%) mean: '); fprintf(fid, '%7.1f', squeeze(Survival_EGSO_mean(map_idx,si,:))); fprintf(fid, '\n');
    end
    len_cost_pct = 100*(VMean_EGSO(map_idx,2)-VMean_EGSO(map_idx,1))/VMean_EGSO(map_idx,1);
    fprintf(fid, '  Length cost of min_clearance 1.0->2.0: %.2f%%\n\n', len_cost_pct);
end
fclose(fid);
fprintf('Saved: results/noise/NoiseRobustness_v2_Data.mat / Summary.txt\n');

fprintf('\n========================================\nDone\n========================================\n');


%% ════════════════════════════════════════════════════════════════════════════
function [meanClr, minClr] = pathClearanceStats(path, T_clear, M)
    n = size(path, 1);
    vals = [];
    for i = 1:n-1
        p1 = path(i,:); p2 = path(i+1,:);
        dist = norm(p2 - p1);
        nsteps = max(2, ceil(dist / 0.5));
        for t = linspace(0, 1, nsteps)
            p = p1 + t * (p2 - p1);
            r = max(1, min(M, round(p(1))));
            c = max(1, min(M, round(p(2))));
            v = T_clear(r, c);
            if isfinite(v)
                vals(end+1) = v; %#ok<AGROW>
            end
        end
    end
    if isempty(vals)
        meanClr = NaN; minClr = NaN;
    else
        meanClr = mean(vals);
        minClr  = min(vals);
    end
end

function sr = monteCarloSurvival(path, G, sigma, num_trials)
    M = size(G, 1);
    n = size(path, 1);
    survive_count = 0;
    for trial = 1:num_trials
        noisy = path + sigma * randn(n, 2);
        noisy(1,:) = path(1,:);
        noisy(end,:) = path(end,:);
        ok = true;
        for i = 1:n-1
            if ~collisionFreeSegmentCheck(noisy(i,:), noisy(i+1,:), G, M)
                ok = false;
                break;
            end
        end
        if ok
            survive_count = survive_count + 1;
        end
    end
    sr = 100 * survive_count / num_trials;
end

function ok = collisionFreeSegmentCheck(p1, p2, G, M)
    dist = norm(p2 - p1);
    nsteps = max(2, ceil(dist / 0.5));
    for t = linspace(0, 1, nsteps)
        p = p1 + t * (p2 - p1);
        r = round(p(1)); c = round(p(2));
        if r < 1 || r > M || c < 1 || c > M
            ok = false; return;
        end
        if G(r, c) == 1
            ok = false; return;
        end
    end
    ok = true;
end
