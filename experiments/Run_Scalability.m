% =========================================================================
% Run_Scalability.m
%
% Scalability study. Pushes the map size far beyond the 90x90 of the benchmark set
% and measures:
%   1. how the EVF (T_nav + T_clear) precomputation time actually grows with M,
%      compared with the theoretical worst case O((M*M)^2)
%   2. run time, success rate and path quality of EGSO and two representative
%      baselines (GWO, SSA) as the map grows
%   3. the run time of A* (deterministic graph-search reference) against the
%      growth of EGSO's precomputation time, i.e. whether the comparison with A*
%      still holds at large M
%
% =========================================================================
clc; clear; close all;

if ~exist('results/scalability', 'dir')
    mkdir('results/scalability');
end

fprintf('========================================\n');
fprintf('Scalability study\n');
fprintf('========================================\n\n');

sizes           = [90, 120, 150, 180, 220];
num_runs        = 5;
max_gen         = 100;
num_population  = 50;
collision_threshold = 1e4;

Algo_Names    = {'EGSO', 'GWO', 'SSA'};
num_algos     = numel(Algo_Names);

Precompute_Time = zeros(size(sizes));
AStar_Length    = zeros(size(sizes));
AStar_Time      = zeros(size(sizes));
SR_data     = zeros(num_algos, numel(sizes));
Mean_data   = nan(num_algos, numel(sizes));
Best_data   = nan(num_algos, numel(sizes));
Time_data   = zeros(num_algos, numel(sizes));

for si = 1:numel(sizes)
    M = sizes(si);
    seed = 77000 + M;
    G = GenerateScalabilityMap(M, seed);

    goalMask = false(M, M); goalMask(M, M) = true;
    T_check = FastMarchingField(G, goalMask, 'navigate');
    attempt = 0;
    while ~isfinite(T_check(1,1)) && attempt < 5
        attempt = attempt + 1;
        seed = seed + 5000;
        G = GenerateScalabilityMap(M, seed);
        T_check = FastMarchingField(G, goalMask, 'navigate');
    end
    fprintf('===== M=%d (feasible=%d, density=%.1f%%) =====\n', M, isfinite(T_check(1,1)), 100*sum(G(:))/numel(G));

    % ---- EVF precompute timing (isolated measurement) ----
    obstacleMask = (G == 1);
    tic;
    T_nav_iso = FastMarchingField(G, goalMask, 'navigate'); %#ok<NASGU>
    T_clear_iso = FastMarchingField(G, obstacleMask, 'clearance'); %#ok<NASGU>
    Precompute_Time(si) = toc;
    fprintf('  EVF precompute (T_nav+T_clear): %.3fs\n', Precompute_Time(si));

    % ---- A* deterministic reference ----
    S = [1,1]; E = [M,M];
    tic;
    [pathA, ~] = AStarPath(G, S, E);
    AStar_Time(si) = toc;
    if ~isempty(pathA)
        pathA_s = GenerateSmoothPath(pathA, G);
        pathA_s = GenerateSmoothPath(pathA_s, G);
        lenA = 0;
        for i = 1:size(pathA_s,1)-1
            lenA = lenA + norm(pathA_s(i+1,:) - pathA_s(i,:));
        end
    else
        lenA = Inf;
    end
    AStar_Length(si) = lenA;
    fprintf('  A*: length=%.2f, time=%.3fs\n', lenA, AStar_Time(si));

    % ---- Stochastic algorithms ----
    dim = M - 2;
    fobj = @(x) fitness(x, G);
    for a = 1:num_algos
        fn = str2func(Algo_Names{a});
        fits = zeros(num_runs, 1);
        times = zeros(num_runs, 1);
        for run = 1:num_runs
            tic;
            [Best_score, ~, ~] = fn(num_population, max_gen, 1, M, dim, fobj, G);
            times(run) = toc;
            fits(run) = Best_score;
        end
        valid = fits < collision_threshold;
        SR_data(a, si)   = sum(valid) / num_runs * 100;
        Time_data(a, si) = mean(times);
        if any(valid)
            Mean_data(a, si) = mean(fits(valid));
            Best_data(a, si) = min(fits(valid));
        end
        fprintf('  %-6s: SR=%6.1f%%  V.Mean=%9.2f  V.Best=%9.2f  Time=%7.2fs\n', ...
                Algo_Names{a}, SR_data(a,si), Mean_data(a,si), Best_data(a,si), Time_data(a,si));
    end
    fprintf('\n');
end

save('results/scalability/Scalability_Data.mat', 'sizes', 'num_runs', 'max_gen', 'num_population', ...
     'Precompute_Time', 'AStar_Length', 'AStar_Time', 'Algo_Names', 'SR_data', 'Mean_data', 'Best_data', 'Time_data');

%% -- Text summary -----------------------------------------------------------
fid = fopen('results/scalability/Scalability_Summary.txt', 'w');
fprintf(fid, 'Scalability experiment: map size M in %s, %d runs per stochastic algorithm\n\n', mat2str(sizes), num_runs);
fprintf(fid, '%-6s %14s %12s %12s', 'M', 'Precompute(s)', 'A*len', 'A*time(s)');
for a = 1:num_algos
    fprintf(fid, ' %10s.SR%% %10s.Mean %10s.Time', Algo_Names{a}, Algo_Names{a}, Algo_Names{a});
end
fprintf(fid, '\n');
for si = 1:numel(sizes)
    fprintf(fid, '%-6d %14.3f %12.2f %12.3f', sizes(si), Precompute_Time(si), AStar_Length(si), AStar_Time(si));
    for a = 1:num_algos
        fprintf(fid, ' %14.1f %14.2f %14.2f', SR_data(a,si), Mean_data(a,si), Time_data(a,si));
    end
    fprintf(fid, '\n');
end
fclose(fid);
fprintf('Saved: results/scalability/Scalability_Summary.txt\n');

%% -- Figure 1: precomputation time vs M (log-log, with an O(M^4) reference) --
fig1 = figure('Position',[100,100,800,550],'Visible','off','Color','w');
loglog(sizes, Precompute_Time, '-o', 'Color', [0.85 0.10 0.10], 'LineWidth', 2, 'MarkerFaceColor', [0.85 0.10 0.10]);
hold on;
ref = Precompute_Time(1) * (sizes/sizes(1)).^4;
loglog(sizes, ref, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
p = polyfit(log(sizes), log(Precompute_Time), 1);
loglog(sizes, exp(p(2)) * sizes.^p(1), ':', 'Color', [0.2 0.45 0.75], 'LineWidth', 1.5);
xlabel('Map size M (M\times M grid)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('EVF precompute time (s)', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Measured', sprintf('Worst-case O(M^4) reference'), sprintf('Empirical fit O(M^{%.2f})', p(1))}, ...
       'Location', 'northwest', 'FontSize', 10);
title('EVF precomputation cost vs. map scale', 'FontSize', 13, 'FontWeight', 'bold');
grid on;
print(fig1, 'results/scalability/Fig_PrecomputeTime_vs_M', '-dpng', '-r300');
close(fig1);
fprintf('Saved: results/scalability/Fig_PrecomputeTime_vs_M.png (empirical exponent = %.2f)\n', p(1));

%% -- Figure 2: total run time vs M (EGSO/GWO/SSA/A*) -------------------------
fig2 = figure('Position',[100,100,800,550],'Visible','off','Color','w');
colors = [0.90 0.10 0.10; 0.00 0.45 0.74; 0.20 0.60 0.20];
hold on;
for a = 1:num_algos
    plot(sizes, Time_data(a,:), '-o', 'Color', colors(a,:), 'LineWidth', 2, 'MarkerFaceColor', colors(a,:));
end
plot(sizes, AStar_Time, '--s', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
xlabel('Map size M', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Average total run time (s)', 'FontSize', 12, 'FontWeight', 'bold');
legend([Algo_Names, {'A*'}], 'Location', 'northwest', 'FontSize', 10);
title('Total algorithm run time vs. map scale', 'FontSize', 13, 'FontWeight', 'bold');
grid on;
print(fig2, 'results/scalability/Fig_RunTime_vs_M', '-dpng', '-r300');
close(fig2);
fprintf('Saved: results/scalability/Fig_RunTime_vs_M.png\n');

%% -- Figure 3: SR / gap to A* vs M --------------------------------------------
fig3 = figure('Position',[100,100,1400,550],'Visible','off','Color','w');
subplot(1,2,1);
hold on;
for a = 1:num_algos
    plot(sizes, SR_data(a,:), '-o', 'Color', colors(a,:), 'LineWidth', 2, 'MarkerFaceColor', colors(a,:));
end
xlabel('Map size M', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Success Rate (%)', 'FontSize', 12, 'FontWeight', 'bold');
legend(Algo_Names, 'Location', 'southwest', 'FontSize', 10);
title('Success rate vs. map scale', 'FontSize', 12, 'FontWeight', 'bold');
ylim([-5 105]); grid on;

subplot(1,2,2);
gap_pct = 100 * (Mean_data(1,:) - AStar_Length) ./ AStar_Length;
plot(sizes, gap_pct, '-o', 'Color', [0.90 0.10 0.10], 'LineWidth', 2, 'MarkerFaceColor', [0.90 0.10 0.10]);
xlabel('Map size M', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('EGSO V.Mean gap to A* (%)', 'FontSize', 12, 'FontWeight', 'bold');
title('EGSO--A* optimality gap vs. map scale', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
sgtitle('Path quality and success rate hold up as map scale increases', 'FontSize', 13, 'FontWeight', 'bold');
print(fig3, 'results/scalability/Fig_Quality_vs_M', '-dpng', '-r300');
close(fig3);
fprintf('Saved: results/scalability/Fig_Quality_vs_M.png\n');

fprintf('\n========================================\n');
fprintf('Scalability study finished\n');
fprintf('========================================\n');
