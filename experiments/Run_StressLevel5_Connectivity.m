% =========================================================================
% Run_StressLevel5_Connectivity.m
% Main_LocalMinima_StressTest.m reports SR and V.Mean of EGSO/EGSO_noEVF/GWO/WOA at
% complexity levels 1-5. This script adds (a) whether the start (1,1) and goal (M,M)
% are connected at each level and (b) the deterministic A* and FMM-descent results.
% At level 5 every method has 0% success; this separates "method failure" (feasible
% map, path not found) from "infeasible instance" (start and goal not connected).
%
% Uses exactly the map generation, seeds and feasibility-retry logic of
% Main_LocalMinima_StressTest.m, level by level:
%   1. FastMarchingField checks whether T(1,1) is finite (= connectivity)
%   2. deterministic A* and FMM-descent (4-connected greedy descent on the same T field)
% =========================================================================
clc; clear;

if ~exist('results/stress', 'dir')
    mkdir('results/stress');
end

M = 60;
complexity_levels = [1, 2, 3, 4, 5];
S = [1, 1]; E = [M, M];

fprintf('========================================\n');
fprintf('Complexity levels 1-5: connectivity + A*/FMM-descent baselines\n');
fprintf('========================================\n\n');

Connectivity = zeros(1, length(complexity_levels));
T_at_start   = nan(1, length(complexity_levels));
AStar_Length = nan(1, length(complexity_levels));
AStar_Time   = nan(1, length(complexity_levels));
FMM_Length   = nan(1, length(complexity_levels));
FMM_Time     = nan(1, length(complexity_levels));

for ci = 1:length(complexity_levels)
    nstages = complexity_levels(ci);
    seed = 90000 + nstages;
    G = GenerateDecoyStressMap(M, nstages, seed);

    % Same feasibility-retry logic as Main_LocalMinima_StressTest.m, so the maps checked
    % here are exactly those of the stress test
    goalMask = false(M, M); goalMask(M, M) = true;
    T = FastMarchingField(G, goalMask, 'navigate');
    attempt = 0;
    while ~isfinite(T(1,1)) && attempt < 5
        attempt = attempt + 1;
        seed = seed + 1000;
        G = GenerateDecoyStressMap(M, nstages, seed);
        T = FastMarchingField(G, goalMask, 'navigate');
    end

    Connectivity(ci) = isfinite(T(1,1));
    T_at_start(ci)   = T(1,1);

    fprintf('Complexity level %d: connected=%d, T_nav(1,1)=%.2f\n', nstages, Connectivity(ci), T_at_start(ci));

    if Connectivity(ci)
        tic; [pathA, ~] = AStarPath(G, S, E); AStar_Time(ci) = toc;
        if ~isempty(pathA)
            pA = GenerateSmoothPath(pathA, G); pA = GenerateSmoothPath(pA, G);
            len = 0; for i=1:size(pA,1)-1, len = len + norm(pA(i+1,:)-pA(i,:)); end
            AStar_Length(ci) = len;
        end

        tic; [pathF, ~] = FMMDescentPath(G, S, E); FMM_Time(ci) = toc;
        if ~isempty(pathF)
            pF = GenerateSmoothPath(pathF, G); pF = GenerateSmoothPath(pF, G);
            len = 0; for i=1:size(pF,1)-1, len = len + norm(pF(i+1,:)-pF(i,:)); end
            FMM_Length(ci) = len;
        end
        fprintf('  A*: length=%.2f (%.3fs)   FMM-descent: length=%.2f (%.3fs)\n', ...
            AStar_Length(ci), AStar_Time(ci), FMM_Length(ci), FMM_Time(ci));
    else
        fprintf('  [infeasible instance] start and goal are not connected at this level; A*/FMM-descent skipped\n');
    end
end

save('results/stress/StressLevel5_Connectivity.mat', 'complexity_levels', 'Connectivity', ...
     'T_at_start', 'AStar_Length', 'AStar_Time', 'FMM_Length', 'FMM_Time');

fid = fopen('results/stress/StressLevel5_Connectivity.txt', 'w');
fprintf(fid, 'Start-Goal Connectivity and Deterministic Baselines by Stress-Test Complexity Level\n');
fprintf(fid, '=========================================================================\n\n');
fprintf(fid, '%-10s %-14s %-14s %-14s %-14s %-14s\n', 'Level', 'Connected', 'A*.Length', 'A*.Time(s)', 'FMM.Length', 'FMM.Time(s)');
for ci = 1:length(complexity_levels)
    fprintf(fid, '%-10d %-14d %-14.2f %-14.3f %-14.2f %-14.3f\n', complexity_levels(ci), Connectivity(ci), ...
        AStar_Length(ci), AStar_Time(ci), FMM_Length(ci), FMM_Time(ci));
end
fclose(fid);
fprintf('\nSaved: results/stress/StressLevel5_Connectivity.mat / .txt\n');
