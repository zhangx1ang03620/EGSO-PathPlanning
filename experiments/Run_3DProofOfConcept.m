% =========================================================================
% Run_3DProofOfConcept.m
%
% Minimal proof-of-concept: transfers the EGSO paradigm (EVF + GSI + elite-scout
% search + CALS) from 2D grids to 3D voxel scenes.
% Note: this is a deliberately minimal demo with a simplified 3D path encoding
%
% (straight-segment collision checks between columns, see CollisionFreeSegment3D.m)
% rather than a full port of the cell-by-cell connectivity repair of the 2D pipeline;
% it tests transferability, not a 3D benchmark of the same rigour as the 2D study.
%
% =========================================================================
clc; clear; close all;

if ~exist('results/3d', 'dir')
    mkdir('results/3d');
end

fprintf('========================================\n');
fprintf('3D voxel proof-of-concept\n');
fprintf('========================================\n\n');

R = 20; C = 20; H = 8;
seed = 4201;
G = GenerateEGSO3DMap(R, C, H, seed);
S = [1,1,1]; E = [R,C,H];

goalMask = false(R,C,H); goalMask(E(1),E(2),E(3)) = true;
T_check = FastMarchingField3D(G, goalMask, 'navigate');
fprintf('map size: %dx%dx%d, obstacle density=%.1f%%, feasible=%d (T(S)=%.2f)\n', ...
        R, C, H, 100*sum(G(:))/numel(G), isfinite(T_check(1,1,1)), T_check(1,1,1));

dim = C - 2;
fobj = @(x) fitness3D(x, G);
num_runs = 10;
max_gen = 100;
num_population = 50;
collision_threshold = 1e4;

%% -- FMM-descent3D (deterministic reference, 1 run) -------------------------
tic;
[path_fmm, len_fmm] = FMMDescentPath3D(G, S, E);
t_fmm = toc;
fprintf('\nFMM-descent3D: length=%.2f, time=%.3fs, success=%d\n', len_fmm, t_fmm, ~isempty(path_fmm));

%% -- EGSO3D / GWO3D (stochastic, %d independent runs) -----------------------
Algo_Names = {'EGSO3D', 'GWO3D'};
Results3D = struct();
BestPaths = struct();

for a = 1:numel(Algo_Names)
    name = Algo_Names{a};
    fn = str2func(name);
    fits = zeros(num_runs,1);
    times = zeros(num_runs,1);
    best_pos_all = cell(num_runs,1);
    for run = 1:num_runs
        tic;
        [bs, bp, ~] = fn(num_population, max_gen, dim, fobj, G);
        times(run) = toc;
        fits(run) = bs;
        best_pos_all{run} = bp;
        fprintf('  [%s] run %d/%d: fitness=%.2f, time=%.2fs\n', name, run, num_runs, bs, times(run));
    end
    valid = fits < collision_threshold;
    Results3D.(name).SR = 100*sum(valid)/num_runs;
    if any(valid)
        Results3D.(name).Mean = mean(fits(valid));
        Results3D.(name).Best = min(fits(valid));
        Results3D.(name).Worst = max(fits(valid));
        Results3D.(name).Std  = std(fits(valid));
        [~, best_run_local] = min(fits + (~valid)*1e9);
        bp = best_pos_all{best_run_local};
        row = round(bp(1:dim)); hgt = round(bp(dim+1:end));
        waypoints = [S; [row(:), (2:dim+1)', hgt(:)]; E];
        BestPaths.(name) = GenerateSmoothPath3D(waypoints, G);
    else
        Results3D.(name).Mean = Inf; Results3D.(name).Best = Inf;
        Results3D.(name).Worst = Inf; Results3D.(name).Std = 0;
        BestPaths.(name) = [];
    end
    Results3D.(name).Time = mean(times);
    fprintf('  >> %s: SR=%.1f%%  V.Best=%.2f  V.Mean=%.2f  V.Worst=%.2f  Time=%.2fs\n\n', ...
            name, Results3D.(name).SR, Results3D.(name).Best, Results3D.(name).Mean, Results3D.(name).Worst, Results3D.(name).Time);
end

save('results/3d/ProofOfConcept3D_Data.mat', 'G', 'S', 'E', 'Results3D', 'BestPaths', ...
     'path_fmm', 'len_fmm', 't_fmm', 'num_runs', 'max_gen', 'num_population', 'R','C','H');

%% -- Text summary -----------------------------------------------------------
fid = fopen('results/3d/ProofOfConcept3D_Summary.txt', 'w');
fprintf(fid, '3D Voxel Proof-of-Concept: %dx%dx%d, density=%.1f%%, %d independent runs\n\n', ...
        R,C,H, 100*sum(G(:))/numel(G), num_runs);
fprintf(fid, '%-12s %10s %10s %10s %10s %10s %10s\n', 'Algorithm', 'SR(%)', 'V.Best', 'V.Mean', 'V.Worst', 'V.Std', 'Time(s)');
fprintf(fid, '%-12s %10s %10.2f %10s %10s %10s %10.3f\n', 'FMM-descent', '---', len_fmm, '---', '---', '---', t_fmm);
for a = 1:numel(Algo_Names)
    name = Algo_Names{a};
    fprintf(fid, '%-12s %10.1f %10.2f %10.2f %10.2f %10.2f %10.3f\n', name, Results3D.(name).SR, ...
        Results3D.(name).Best, Results3D.(name).Mean, Results3D.(name).Worst, Results3D.(name).Std, Results3D.(name).Time);
end
fclose(fid);
fprintf('Saved: results/3d/ProofOfConcept3D_Summary.txt\n');

%% -- LaTeX table ------------------------------------------------------------
fid = fopen('results/3d/ProofOfConcept3D_LaTeX.txt', 'w');
fprintf(fid, '%% Table: Minimal 3D voxel proof-of-concept (%d runs, %dx%dx%d grid)\n', num_runs, R,C,H);
fprintf(fid, '\\begin{table}[h]\n\\caption{Minimal 3D voxel proof-of-concept: EGSO3D vs. FMM-descent3D vs. GWO3D.}\\label{tab:3dpoc}\n');
fprintf(fid, '\\begin{tabular*}{\\textwidth}{@{\\extracolsep\\fill}lcccc@{}}\n\\toprule\n');
fprintf(fid, 'Algorithm & SR (\\%%) & V.Best & V.Mean & Time (s)\\\\\n\\midrule\n');
fprintf(fid, 'FMM-descent3D & --- & %.2f & --- & %.3f\\\\\n', len_fmm, t_fmm);
for a = 1:numel(Algo_Names)
    name = Algo_Names{a};
    fprintf(fid, '%s & %.1f & %.2f & %.2f & %.3f\\\\\n', name, Results3D.(name).SR, Results3D.(name).Best, Results3D.(name).Mean, Results3D.(name).Time);
end
fprintf(fid, '\\botrule\n\\end{tabular*}\n\\end{table}\n');
fclose(fid);
fprintf('Saved: results/3d/ProofOfConcept3D_LaTeX.txt\n');

%% -- 3D view: obstacle voxels + best path of each of the three methods -------
fig = figure('Position', [80,80,900,750], 'Visible', 'off', 'Color', 'w');
hold on;

[or_, oc_, oh_] = ind2sub(size(G), find(G==1));
scatter3(oc_, or_, oh_, 60, [0.6 0.6 0.6], 'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor', 'none');

if ~isempty(path_fmm)
    plot3(path_fmm(:,2), path_fmm(:,1), path_fmm(:,3), '--', 'Color', [0.2 0.45 0.75], 'LineWidth', 2);
end
if isfield(BestPaths,'GWO3D') && ~isempty(BestPaths.GWO3D)
    p = BestPaths.GWO3D;
    plot3(p(:,2), p(:,1), p(:,3), '-.', 'Color', [0.2 0.6 0.2], 'LineWidth', 2);
end
if isfield(BestPaths,'EGSO3D') && ~isempty(BestPaths.EGSO3D)
    p = BestPaths.EGSO3D;
    plot3(p(:,2), p(:,1), p(:,3), '-', 'Color', [0.90 0.10 0.10], 'LineWidth', 3);
end

plot3(S(2),S(1),S(3),'go','MarkerSize',10,'MarkerFaceColor','g');
plot3(E(2),E(1),E(3),'bs','MarkerSize',10,'MarkerFaceColor','b');

xlabel('Column (traversal axis)', 'FontSize', 11, 'FontWeight','bold');
ylabel('Row', 'FontSize', 11, 'FontWeight','bold');
zlabel('Height', 'FontSize', 11, 'FontWeight','bold');
legend({'Obstacle voxels', sprintf('FMM-descent3D (%.1f)', len_fmm), ...
        sprintf('GWO3D (%.1f)', Results3D.GWO3D.Best), sprintf('EGSO3D (%.1f)', Results3D.EGSO3D.Best)}, ...
        'Location', 'northeast', 'FontSize', 9);
title('3D voxel proof-of-concept: obstacle blocks and best paths found', 'FontSize', 12, 'FontWeight', 'bold');
view(-35, 22);
grid on; box on;
set(gcf, 'Color', 'w');
print(fig, 'results/3d/Fig_3DProofOfConcept', '-dpng', '-r300');
close(fig);
fprintf('Saved: results/3d/Fig_3DProofOfConcept.png\n');

fprintf('\n========================================\n');
fprintf('3D proof-of-concept finished\n');
fprintf('========================================\n');
