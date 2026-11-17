% QUICKSTART  Run EGSO on one benchmark map and plot the result.
%
% Solves Map 1 (single decoy dead-end trap, 50x50) with the paper setting
% N = 50, T = 100, and shows the Eikonal navigation field, the local-minima
% census of that field, and the best path found. Takes about 5 s.
%
%   >> cd <repository root>
%   >> setup_egso
%   >> quickstart

setup_egso();
map_idx = 1;
[Maps, Map_Info] = GenerateEGSOMaps();
close(gcf);                               % GenerateEGSOMaps draws all maps
G = Maps{map_idx};
M = size(G, 1);
goal = [M, M];

%% 1. Eikonal viscosity field (EVF) and its local-minima census
goalMask = false(size(G)); goalMask(goal(1), goal(2)) = true;
T_nav = FastMarchingField(G, goalMask, 'navigate');
[trap_rate, num_sinks] = GreedyDescentCensus(T_nav, G, goal);
fprintf('EVF on Map %d: trap rate = %.2f%%, spurious sinks = %d\n', map_idx, trap_rate, num_sinks);

%% 2. EGSO
rng(2026);
fobj = @(x) fitness(x, G);
tic;
[best_len, best_x, curve] = EGSO(50, 100, 1, M, M - 2, fobj, G);
fprintf('EGSO on Map %d: path length = %.4f (%.1f s)\n', map_idx, best_len, toc);

%% 3. Decode the decision vector into the smoothed path
route = [1, round(best_x), M];
path = generateContinuousRoute(route, G);
path = GenerateSmoothPath(GenerateSmoothPath(path, G), G);

%% 4. Plot
figure('Color', 'w', 'Position', [100 100 1250 420]);
subplot(1, 3, 1);
Tshow = T_nav; Tshow(~isfinite(Tshow)) = NaN;
imagesc(Tshow); axis image; colorbar; colormap(gca, parula);
title('Eikonal navigation field T_{nav}');
subplot(1, 3, 2);
imagesc(1 - G); colormap(gca, gray); axis image; hold on;
drawPath(path, [0.85 0.1 0.1]);
title(sprintf('EGSO path, length %.2f', best_len));
subplot(1, 3, 3);
semilogy(curve, 'LineWidth', 1.5); grid on;
xlabel('Iteration'); ylabel('Best fitness'); title('Convergence');
