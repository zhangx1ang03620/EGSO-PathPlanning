% =========================================================================
% Local-minima sensitivity stress test
%
% Sweeps the number of decoy stages (1..5, the local-minimum complexity level) on a
% parameterised multi-stage decoy-corridor map and compares EGSO (Full), EGSO w/o EVF
% (APF guidance) and two representative swarm baselines (GWO/WOA) in success rate (SR)
% and path length (V.Mean) as the complexity increases.
%
%
%
% =========================================================================
clc; clear; close all;

if ~exist('results/stress', 'dir')
    mkdir('results/stress');
end

fprintf('========================================\n');
fprintf('EGSO local-minima sensitivity stress test\n');
fprintf('========================================\n\n');

%% Experiment parameters (paper setting)
M               = 60;
complexity_levels = [1, 2, 3, 4, 5];
num_runs        = 15;
max_gen         = 100;
num_population  = 50;
collision_threshold = 1e4;

Algo_Names     = {'EGSO', 'EGSO_noEVF', 'GWO', 'WOA'};
Display_Names  = {'EGSO (Full)', 'EGSO w/o EVF (APF-guided)', 'GWO', 'WOA'};
num_algos      = length(Algo_Names);

colors  = [0.90 0.10 0.10; 0.00 0.45 0.74; 0.20 0.60 0.20; 0.60 0.40 0.00];
markers = {'o', 's', '^', 'd'};

SR_data   = zeros(num_algos, length(complexity_levels));
Mean_data = nan(num_algos, length(complexity_levels));

fprintf('Configuration:\n');
fprintf('  map size: %dx%d\n', M, M);
fprintf('  complexity levels (decoy stages): %s\n', mat2str(complexity_levels));
fprintf('  algorithms: %s\n', strjoin(Algo_Names, ', '));
fprintf('  runs: %d, iterations: %d, population: %d\n\n', num_runs, max_gen, num_population);

%% Main loop over complexity levels
for ci = 1:length(complexity_levels)
    nstages = complexity_levels(ci);
    seed = 90000 + nstages;
    G = GenerateDecoyStressMap(M, nstages, seed);

    % Check feasibility with the EVF fast-marching solver; retry with another seed if infeasible
    goalMask = false(M, M); goalMask(M, M) = true;
    T = FastMarchingField(G, goalMask, 'navigate');
    attempt = 0;
    while ~isfinite(T(1,1)) && attempt < 5
        attempt = attempt + 1;
        seed = seed + 1000;
        G = GenerateDecoyStressMap(M, nstages, seed);
        T = FastMarchingField(G, goalMask, 'navigate');
    end
    fprintf('Complexity level %d (decoy stages): feasible=%d (T(1,1)=%.2f)\n', nstages, isfinite(T(1,1)), T(1,1));

    Xmax = M; X_min = 1; dim = Xmax - 2;
    fobj = @(x) fitness(x, G);

    for a = 1:num_algos
        fn = str2func(Algo_Names{a});
        fits = zeros(num_runs, 1);
        for run = 1:num_runs
            [Best_score, ~, ~] = fn(num_population, max_gen, X_min, Xmax, dim, fobj, G);
            fits(run) = Best_score;
        end
        valid = fits < collision_threshold;
        SR_data(a, ci) = sum(valid) / num_runs * 100;
        if any(valid)
            Mean_data(a, ci) = mean(fits(valid));
        end
        fprintf('  [%d decoy stages] %-14s: SR=%6.1f%%, V.Mean=%8.2f\n', ...
                nstages, Algo_Names{a}, SR_data(a, ci), Mean_data(a, ci));
    end
end

%% Map figure (increasing complexity)
fig_maps = figure('Position', [50, 50, 1600, 350], 'Visible', 'off', 'Color', 'w');
for ci = 1:length(complexity_levels)
    nstages = complexity_levels(ci);
    G = GenerateDecoyStressMap(M, nstages, 90000+nstages);
    subplot(1, length(complexity_levels), ci);
    imagesc(G);
    colormap(gca, [1 1 1; 0.20 0.45 0.70]);
    axis equal tight;
    set(gca, 'XTick', [], 'YTick', []);
    title(sprintf('%d-stage decoy', nstages), 'FontSize', 10, 'FontWeight', 'bold');
end
sgtitle('Stress-Test Maps: Increasing Local-Minima Complexity', 'FontSize', 13, 'FontWeight', 'bold');
print(fig_maps, 'results/stress/StressTest_Maps', '-dpng', '-r300');
close(fig_maps);

%% SR versus complexity
fig1 = figure('Position', [100, 100, 800, 500], 'Visible', 'off');
hold on;
for a = 1:num_algos
    plot(complexity_levels, SR_data(a, :), '-', 'Color', colors(a, :), 'LineWidth', 2, ...
         'Marker', markers{a}, 'MarkerSize', 8, 'MarkerFaceColor', colors(a, :));
end
xlabel('Decoy-trap stages (local-minima complexity level)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Success Rate (%)', 'FontSize', 12, 'FontWeight', 'bold');
title('Local-Minima Susceptibility Stress Test: Success Rate', 'FontSize', 13, 'FontWeight', 'bold');
legend(Display_Names, 'Location', 'southwest', 'FontSize', 10);
xticks(complexity_levels);
ylim([-5 105]);
grid on;
set(gcf, 'Color', 'w');
print(fig1, 'results/stress/StressTest_SR', '-dpng', '-r300');
close(fig1);

%% V.Mean versus complexity
fig2 = figure('Position', [100, 100, 800, 500], 'Visible', 'off');
hold on;
for a = 1:num_algos
    plot(complexity_levels, Mean_data(a, :), '-', 'Color', colors(a, :), 'LineWidth', 2, ...
         'Marker', markers{a}, 'MarkerSize', 8, 'MarkerFaceColor', colors(a, :));
end
xlabel('Decoy-trap stages (local-minima complexity level)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Valid Mean Path Length', 'FontSize', 12, 'FontWeight', 'bold');
title('Local-Minima Susceptibility Stress Test: Path Quality', 'FontSize', 13, 'FontWeight', 'bold');
legend(Display_Names, 'Location', 'northwest', 'FontSize', 10);
xticks(complexity_levels);
grid on;
set(gcf, 'Color', 'w');
print(fig2, 'results/stress/StressTest_VMean', '-dpng', '-r300');
close(fig2);

%% Save data and text summary
save('results/stress/StressTest_Data.mat', 'SR_data', 'Mean_data', 'complexity_levels', ...
     'Algo_Names', 'num_runs', 'max_gen', 'num_population', 'M');

fid = fopen('results/stress/StressTest_Summary.txt', 'w');
fprintf(fid, 'EGSO Local-Minima Susceptibility Stress Test\n');
fprintf(fid, 'Config: %d runs, %d iterations, pop=%d, map size=%dx%d\n\n', ...
        num_runs, max_gen, num_population, M, M);
fprintf(fid, '%-10s', 'Stages');
for a = 1:num_algos
    fprintf(fid, ' %28s', Display_Names{a});
end
fprintf(fid, '\n');
for ci = 1:length(complexity_levels)
    fprintf(fid, '%-10d', complexity_levels(ci));
    for a = 1:num_algos
        fprintf(fid, '   SR=%6.1f%%, V.Mean=%8.2f', SR_data(a, ci), Mean_data(a, ci));
    end
    fprintf(fid, '\n');
end
fclose(fid);

fprintf('\n========================================\n');
fprintf('Stress test finished. Results saved in results/stress/\n');
fprintf('========================================\n');
