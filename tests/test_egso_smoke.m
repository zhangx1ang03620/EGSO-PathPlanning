clc; clear;
rng(123);

%% Build a simple 40x40 test map (temporary, not a paper map) with a U-trap
M = 40;
G = zeros(M, M);
G(rand(M,M) < 0.18) = 1;
% carve a concave U-trap near the middle to stress-test EVF guidance
G(10:28, 20) = 1;
G(10, 20:32) = 1;
G(28, 20:32) = 1;
G(10:28, 32) = 1;
G(26:28, 32) = 0;  % small opening near the bottom
% safety zones
G(1:5,1:5) = 0;
G(M-4:M, M-4:M) = 0;

Xmax = M; X_min = 1; dim = Xmax - 2;
fobj = @(x) fitness(x, G);

num_runs = 3;
fprintf('Running EGSO smoke test: %d runs, N=30, T=60, on 40x40 U-trap map\n', num_runs);

for run = 1:num_runs
    tic;
    [Best_score, Best_pos, Curve] = EGSO(30, 60, X_min, Xmax, dim, fobj, G);
    elapsed = toc;

    Best_pos = round(Best_pos);
    Best_pos = max(Best_pos, 1);
    Best_pos = min(Best_pos, Xmax);
    route = [1 Best_pos Xmax];
    path = generateContinuousRoute(route, G);

    is_collision = Best_score >= 1e4;
    fprintf('  Run %d: fitness=%.4f (%s), time=%.2fs, curve monotonic=%d\n', ...
        run, Best_score, ternary(is_collision,'COLLISION','collision-free'), elapsed, ...
        all(diff(Curve) <= 1e-9));

    assert(~is_collision, 'EGSO must find a collision-free path on this simple test map');
    assert(Best_score > 30 && Best_score < 500, 'fitness should be in a plausible path-length range');
    assert(~isempty(path), 'generateContinuousRoute must produce a valid path');
end

fprintf('\nALL EGSO SMOKE TESTS PASSED.\n');

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
