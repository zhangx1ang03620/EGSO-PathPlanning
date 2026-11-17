clc; clear;
rng(555);

M = 40;
G = zeros(M, M);
G(rand(M,M) < 0.18) = 1;
G(1:5,1:5) = 0; G(M-4:M, M-4:M) = 0;

Xmax = M; X_min = 1; dim = Xmax - 2;
fobj = @(x) fitness(x, G);

names = {'ARO', 'INFO', 'SBOA', 'GWO', 'PSO'};
for k = 1:numel(names)
    fn = str2func(names{k});
    tic;
    [Best_score, Best_pos, Curve] = fn(30, 60, X_min, Xmax, dim, fobj, G);
    elapsed = toc;
    fprintf('%-8s: fitness=%.4f, time=%.2fs, len(curve)=%d, monotonic=%d\n', ...
        names{k}, Best_score, elapsed, numel(Curve), all(diff(Curve) <= 1e-9));
    assert(numel(Best_pos) == dim);
    assert(numel(Curve) == 60);
end
fprintf('\nALL NEW BASELINE SMOKE TESTS PASSED.\n');
