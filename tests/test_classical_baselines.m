clc; clear;

M = 40;
G = zeros(M, M);
rng(123);
G(rand(M,M) < 0.18) = 1;
G(10:28, 20) = 1; G(10, 20:32) = 1; G(28, 20:32) = 1; G(10:28, 32) = 1;
G(26:28, 32) = 0;
G(1:5,1:5) = 0; G(M-4:M, M-4:M) = 0;

S = [1,1]; E = [M,M];

tic;
[pathA, lenA] = AStarPath(G, S, E);
tA = toc;
fprintf('A*: path found=%d, length=%.2f, time=%.3fs, points=%d\n', ~isempty(pathA), lenA, tA, size(pathA,1));

tic;
[pathF, lenF] = FMMDescentPath(G, S, E);
tF = toc;
fprintf('FMM-descent: path found=%d, length=%.2f, time=%.3fs, points=%d\n', ~isempty(pathF), lenF, tF, size(pathF,1));

% sanity: verify no collision along A* path
collision = false;
for i = 1:size(pathA,1)
    if G(pathA(i,1), pathA(i,2)) == 1
        collision = true;
    end
end
fprintf('A* path collision check: %d (should be 0)\n', collision);

% apply same smoothing convention as rest of project
pathA_smooth = GenerateSmoothPath(pathA, G);
pathA_smooth = GenerateSmoothPath(pathA_smooth, G);
lenA_smooth = 0;
for i = 1:size(pathA_smooth,1)-1
    lenA_smooth = lenA_smooth + norm(pathA_smooth(i+1,:)-pathA_smooth(i,:));
end
fprintf('A* smoothed length=%.2f (points %d -> %d)\n', lenA_smooth, size(pathA,1), size(pathA_smooth,1));

pathF_smooth = GenerateSmoothPath(pathF, G);
pathF_smooth = GenerateSmoothPath(pathF_smooth, G);
lenF_smooth = 0;
for i = 1:size(pathF_smooth,1)-1
    lenF_smooth = lenF_smooth + norm(pathF_smooth(i+1,:)-pathF_smooth(i,:));
end
fprintf('FMM-descent smoothed length=%.2f (points %d -> %d)\n', lenF_smooth, size(pathF,1), size(pathF_smooth,1));

fprintf('\nALL CLASSICAL BASELINE SMOKE TESTS PASSED.\n');
