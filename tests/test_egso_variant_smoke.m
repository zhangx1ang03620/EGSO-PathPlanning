clc; clear;
rng(123);

%% Build the same simple 40x40 U-trap test map as test_egso_smoke.m
M = 40;
G = zeros(M, M);
G(rand(M,M) < 0.18) = 1;
G(10:28, 20) = 1;
G(10, 20:32) = 1;
G(28, 20:32) = 1;
G(10:28, 32) = 1;
G(26:28, 32) = 0;
G(1:5,1:5) = 0;
G(M-4:M, M-4:M) = 0;

Xmax = M; X_min = 1; dim = Xmax - 2;
fobj = @(x) fitness(x, G);

init_modes  = {'gsi', 'random'};
field_modes = {'eikonal', 'apf', 'dijkstra'};

fprintf('Running EGSO_Variant smoke test on all 6 (init_mode, field_mode) combinations\n');
fprintf('40x40 U-trap map, N=30, T=60, 2 runs each\n\n');

for im = 1:numel(init_modes)
    for fm = 1:numel(field_modes)
        init_mode = init_modes{im};
        field_mode = field_modes{fm};
        for run = 1:2
            tic;
            [Best_score, Best_pos, Curve] = EGSO_Variant(30, 60, X_min, Xmax, dim, fobj, G, init_mode, field_mode);
            elapsed = toc;

            is_collision = Best_score >= 1e4;
            fprintf('  [%s + %-8s] run %d: fitness=%.4f (%s), time=%.2fs, monotone=%d\n', ...
                init_mode, field_mode, run, Best_score, ternary(is_collision,'COLLISION','valid'), ...
                elapsed, all(diff(Curve) <= 1e-9));

            assert(~is_collision, 'EGSO_Variant(%s,%s) must find a collision-free path on this simple test map', init_mode, field_mode);
            assert(Best_score > 30 && Best_score < 500, 'fitness should be in a plausible path-length range');
            assert(numel(Best_pos) == dim, 'Best_pos dimension mismatch');
        end
    end
end

fprintf('\nALL EGSO_Variant SMOKE TESTS PASSED (6 combinations).\n');

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
