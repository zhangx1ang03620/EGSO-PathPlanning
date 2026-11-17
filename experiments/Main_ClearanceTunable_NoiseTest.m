% =========================================================================
% Main_ClearanceTunable_NoiseTest.m
% Tests whether raising CALS's min_clearance
% threshold turns EGSO's execution-noise robustness into an explicit,
% tunable dial -- something A* has no equivalent of -- rather than an
% incidental, map-dependent side effect of the default min_clearance=1.0.
%
% For each of the 5 main maps, runs EGSO_ClearanceTunable at
% min_clearance in {1.0 (baseline, should reproduce Main_NoiseRobustness
% .m's EGSO numbers), 2.0 (elevated)}, 10 independent runs each, then
% re-runs the identical Monte Carlo noise simulation (same sigma grid,
% same collision checker) on the best-run path from each setting, and
% reports the length cost of the elevated setting alongside the
% robustness gain.
% =========================================================================
clc; clear;

if ~exist('results/noise', 'dir')
    mkdir('results/noise');
end

[Maps, Map_Info] = GenerateEGSOMaps();
close(gcf);

max_gen        = 100;
num_population = 50;
num_runs       = 10;
collision_threshold = 1e4;
clearance_settings = [1.0, 2.0];
sigma_levels = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0];
num_trials = 300;

SR_data    = zeros(5, numel(clearance_settings));
Mean_data  = nan(5, numel(clearance_settings));
BestPaths  = cell(5, numel(clearance_settings));
MeanClear  = nan(5, numel(clearance_settings));
MinClear   = nan(5, numel(clearance_settings));
Survival   = nan(5, numel(clearance_settings), numel(sigma_levels));

fprintf('========================================\n');
fprintf('Tunable clearance threshold (min_clearance) robustness experiment\n');
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

    for si = 1:numel(clearance_settings)
        mc = clearance_settings(si);
        fits = zeros(num_runs, 1);
        paths = cell(num_runs, 1);
        t0 = tic;
        for run = 1:num_runs
            [Best_score, ~, ~, Best_path] = EGSO_ClearanceTunable(num_population, max_gen, X_min, Xmax, dim, fobj, G, mc);
            fits(run) = Best_score;
            paths{run} = Best_path;
        end
        t_elapsed = toc(t0);

        valid = fits < collision_threshold;
        SR_data(map_idx, si) = sum(valid)/num_runs*100;
        if any(valid)
            valid_idx = find(valid);
            [~, bl] = min(fits(valid_idx));
            best_run = valid_idx(bl);
            Mean_data(map_idx, si) = mean(fits(valid));
            BestPaths{map_idx, si} = paths{best_run};

            [MeanClear(map_idx,si), MinClear(map_idx,si)] = pathClearanceStats(paths{best_run}, T_clear, M);
            for sgi = 1:numel(sigma_levels)
                Survival(map_idx, si, sgi) = monteCarloSurvival(paths{best_run}, G, sigma_levels(sgi), num_trials);
            end
        end

        fprintf('  [min_clearance=%.1f] SR=%6.1f%%, V.Mean=%8.2f, MeanClr=%.3f, MinClr=%.3f  (%.1fs)\n', ...
            mc, SR_data(map_idx,si), Mean_data(map_idx,si), MeanClear(map_idx,si), MinClear(map_idx,si), t_elapsed);
    end
    fprintf('  Survival(%%) @ sigma=%s:\n', mat2str(sigma_levels));
    for si = 1:numel(clearance_settings)
        fprintf('    min_clearance=%.1f: ', clearance_settings(si));
        fprintf('%7.1f', squeeze(Survival(map_idx, si, :)));
        fprintf('\n');
    end
    fprintf('\n');
end
fprintf('Total time: %.1f min\n', toc(t_start)/60);

save('results/noise/ClearanceTunable_Data.mat', 'clearance_settings', 'sigma_levels', 'SR_data', 'Mean_data', ...
     'MeanClear', 'MinClear', 'Survival', 'Map_Info', 'num_runs');

fid = fopen('results/noise/ClearanceTunable_Summary.txt', 'w');
fprintf(fid, 'Tunable min_clearance: length cost vs noise-robustness gain \n');
fprintf(fid, '=========================================================================\n\n');
for map_idx = 1:5
    fprintf(fid, 'Map %d (%s):\n', map_idx, Map_Info(map_idx).size);
    for si = 1:numel(clearance_settings)
        fprintf(fid, '  min_clearance=%.1f: SR=%.1f%%, V.Mean=%.2f, MeanClr=%.3f, MinClr=%.3f\n', ...
            clearance_settings(si), SR_data(map_idx,si), Mean_data(map_idx,si), MeanClear(map_idx,si), MinClear(map_idx,si));
        fprintf(fid, '    Survival(%%) @ sigma=%s: ', mat2str(sigma_levels));
        fprintf(fid, '%7.1f', squeeze(Survival(map_idx, si, :)));
        fprintf(fid, '\n');
    end
    len_cost_pct = 100*(Mean_data(map_idx,2)-Mean_data(map_idx,1))/Mean_data(map_idx,1);
    fprintf(fid, '  Length cost of min_clearance 1.0->2.0: %.2f%%\n\n', len_cost_pct);
end
fclose(fid);
fprintf('Saved: results/noise/ClearanceTunable_Data.mat / Summary.txt\n');
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
