clc; clear;
rng(321);

M = 40;
G = zeros(M, M);
G(rand(M,M) < 0.18) = 1;
G(10:28, 20) = 1; G(10, 20:32) = 1; G(28, 20:32) = 1; G(10:28, 32) = 1;
G(26:28, 32) = 0;
G(1:5,1:5) = 0; G(M-4:M, M-4:M) = 0;

Xmax = M; X_min = 1; dim = Xmax - 2;
fobj = @(x) fitness(x, G);

variants = {'EGSO', 'EGSO_noEVF', 'EGSO_noGSI', 'EGSO_noCALS', 'EGSO_noAntiStag'};

for v = 1:numel(variants)
    name = variants{v};
    fn = str2func(name);
    tic;
    [Best_score, Best_pos, Curve] = fn(30, 60, X_min, Xmax, dim, fobj, G);
    elapsed = toc;
    is_collision = Best_score >= 1e4;
    fprintf('%-16s : fitness=%.4f (%s), time=%.2fs, len(curve)=%d\n', ...
        name, Best_score, tern(is_collision,'COLLISION','ok'), elapsed, numel(Curve));
    assert(numel(Curve) == 60, 'convergence curve length mismatch');
    assert(numel(Best_pos) == dim, 'Best_pos dimension mismatch');
end

fprintf('\nALL ABLATION VARIANTS RAN WITHOUT ERROR.\n');

function out = tern(c,a,b)
    if c, out=a; else, out=b; end
end
