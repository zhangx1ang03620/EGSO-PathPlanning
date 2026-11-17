% =========================================================================
% Main_NoiseRobustness.m
% Execution-noise robustness of EGSO and A* paths.
% Answers, with measurable quantities, when EGSO is preferable to A*.
%
% No re-planning, dynamic obstacles or perception loop: reuses the final A* and EGSO
% paths (results/main/Classical_Baselines.mat, results/main/EGSO_Comparison_
% Results.mat) and the clearance field T_clear already solved by the EVF module:
%   1. (deterministic) CALS smooths under an explicit clearance constraint while
%      A* returns a pure shortest path; is the clearance of EGSO paths larger?
%   2. (Monte Carlo) under execution noise (independent Gaussian perturbation of every
%      waypoint, modelling actuation/localisation error) paths that hug obstacles fail
%      more often; do EGSO paths survive (stay collision-free) more often than A* paths?
%
% This prediction follows directly from the CALS mechanism.
% =========================================================================
clc; clear;

if ~exist('results/noise', 'dir')
    mkdir('results/noise');
end

if ~exist('results/main/Classical_Baselines.mat', 'file') || ~exist('results/main/EGSO_Comparison_Results.mat', 'file')
    error('Main_NoiseRobustness:Missing', 'results/main/Classical_Baselines.mat and results/main/EGSO_Comparison_Results.mat are required');
end

CB = load('results/main/Classical_Baselines.mat');
CR = load('results/main/EGSO_Comparison_Results.mat');
[Maps, Map_Info] = GenerateEGSOMaps();
close(gcf);

collision_threshold = 1e4;
sigma_levels = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0];
num_trials = 300;

MeanClear_AStar = nan(5,1); MinClear_AStar = nan(5,1);
MeanClear_EGSO  = nan(5,1); MinClear_EGSO  = nan(5,1);
SurvivalRate_AStar = nan(5, numel(sigma_levels));
SurvivalRate_EGSO  = nan(5, numel(sigma_levels));

fprintf('========================================\n');
fprintf('Path robustness under execution noise\n');
fprintf('========================================\n\n');

for map_idx = 1:5
    G = Maps{map_idx};
    M = size(G, 1);
    obstacleMask = (G == 1);
    T_clear = FastMarchingField(G, obstacleMask, 'clearance');

    pathA = CB.Classical(map_idx).AStar.Path;   % Nx2 [row col], already smoothed

    valid_idx = find(CR.Results(map_idx).EGSO.Best < collision_threshold);
    [~, best_local] = min(CR.Results(map_idx).EGSO.Best(valid_idx));
    best_run = valid_idx(best_local);
    pathE = CR.Results(map_idx).EGSO.Paths{best_run};  % Nx2 [row col], already smoothed

    %% ---- 1. Deterministic clearance statistics ----
    [MeanClear_AStar(map_idx), MinClear_AStar(map_idx)] = pathClearanceStats(pathA, T_clear, M);
    [MeanClear_EGSO(map_idx),  MinClear_EGSO(map_idx)]  = pathClearanceStats(pathE, T_clear, M);

    fprintf('Map %d (%s):\n', map_idx, Map_Info(map_idx).size);
    fprintf('  Clearance  A*: mean=%.3f min=%.3f   |   EGSO: mean=%.3f min=%.3f\n', ...
        MeanClear_AStar(map_idx), MinClear_AStar(map_idx), MeanClear_EGSO(map_idx), MinClear_EGSO(map_idx));

    %% ---- 2. Monte Carlo survival under execution noise ----
    for si = 1:numel(sigma_levels)
        sigma = sigma_levels(si);
        SurvivalRate_AStar(map_idx, si) = monteCarloSurvival(pathA, G, sigma, num_trials);
        SurvivalRate_EGSO(map_idx, si)  = monteCarloSurvival(pathE, G, sigma, num_trials);
    end
    fprintf('  Survival rate vs sigma:\n');
    fprintf('    sigma      : '); fprintf('%7.1f', sigma_levels); fprintf('\n');
    fprintf('    A*   SR(%%) : '); fprintf('%7.1f', SurvivalRate_AStar(map_idx,:)); fprintf('\n');
    fprintf('    EGSO SR(%%) : '); fprintf('%7.1f', SurvivalRate_EGSO(map_idx,:)); fprintf('\n\n');
end

save('results/noise/NoiseRobustness_Data.mat', 'MeanClear_AStar', 'MinClear_AStar', 'MeanClear_EGSO', 'MinClear_EGSO', ...
     'SurvivalRate_AStar', 'SurvivalRate_EGSO', 'sigma_levels', 'num_trials', 'Map_Info');

%% Text summary
fid = fopen('results/noise/NoiseRobustness_Summary.txt', 'w');
fprintf(fid, 'Execution-Noise Robustness: A* vs EGSO\n');
fprintf(fid, 'Monte Carlo: %d trials per (map, sigma); waypoint-wise isotropic Gaussian perturbation N(0,sigma^2)\n', num_trials);
fprintf(fid, '=========================================================================\n\n');
fprintf(fid, '--- Deterministic clearance margin along nominal (unperturbed) path ---\n');
fprintf(fid, '%-8s %-14s %-14s %-14s %-14s\n', 'Map', 'A*.MeanClr', 'A*.MinClr', 'EGSO.MeanClr', 'EGSO.MinClr');
for map_idx = 1:5
    fprintf(fid, '%-8d %-14.3f %-14.3f %-14.3f %-14.3f\n', map_idx, MeanClear_AStar(map_idx), MinClear_AStar(map_idx), ...
        MeanClear_EGSO(map_idx), MinClear_EGSO(map_idx));
end
fprintf(fid, '\n--- Monte Carlo survival rate (%%) vs execution-noise sigma (grid cells) ---\n');
for map_idx = 1:5
    fprintf(fid, '\nMap %d:\n', map_idx);
    fprintf(fid, '  sigma      : '); fprintf(fid, '%7.1f', sigma_levels); fprintf(fid, '\n');
    fprintf(fid, '  A*   SR(%%) : '); fprintf(fid, '%7.1f', SurvivalRate_AStar(map_idx,:)); fprintf(fid, '\n');
    fprintf(fid, '  EGSO SR(%%) : '); fprintf(fid, '%7.1f', SurvivalRate_EGSO(map_idx,:)); fprintf(fid, '\n');
end
fclose(fid);
fprintf('Saved: results/noise/NoiseRobustness_Data.mat / Summary.txt\n');

%% Figure: survival curves (mean over the five maps + per map)
fig = figure('Position',[100,100,900,600],'Visible','off');
hold on;
avgA = mean(SurvivalRate_AStar, 1);
avgE = mean(SurvivalRate_EGSO, 1);
plot(sigma_levels, avgA, '-o', 'Color', [0.1 0.1 0.1], 'LineWidth', 2, 'MarkerFaceColor', [0.1 0.1 0.1]);
plot(sigma_levels, avgE, '-s', 'Color', [0.90 0.10 0.10], 'LineWidth', 2, 'MarkerFaceColor', [0.90 0.10 0.10]);
xlabel('Execution-noise \sigma (grid cells)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Collision-free survival rate (%)', 'FontSize', 12, 'FontWeight', 'bold');
title('Execution-Noise Robustness: A* vs EGSO (averaged over 5 maps)', 'FontSize', 13, 'FontWeight', 'bold');
legend({'A*', 'EGSO'}, 'Location', 'southwest', 'FontSize', 11);
grid on; set(gcf,'Color','w');
print(fig, 'results/noise/NoiseRobustness_Survival', '-dpng', '-r300');
close(fig);
fprintf('Saved: results/noise/NoiseRobustness_Survival.png\n');

fprintf('\n========================================\n');
fprintf('Execution-noise robustness experiment finished\n');
fprintf('========================================\n');


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
        % Start and goal are not perturbed (noise acts on the motion, not on the known endpoints)
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
