% =========================================================================
% Run_NFE_Comparison.m
% Number of function evaluations (NFE) per run under the fixed iteration budget, counted
% with the FobjCounter wrapper, for five representative competitors: GWO/PSO (one
% evaluation per individual per iteration), ABC (employed + onlooker phases, higher NFE),
% SSA (the strongest baseline overall), ARO (the only algorithm ahead of EGSO on Map 2),
% plus EGSO itself, on Map 1 with 5 independent runs each.
% NOTE: FobjCounter does not see the direct fitness() calls inside LocalSearch.m; the
% complete counts (Table B1 of the paper) come from experiments/equal_nfe/Run_NFE_T100.m.
% =========================================================================
clc; clear;

[Maps, Map_Info] = GenerateEGSOMaps();
close(gcf);

G = Maps{1};
M = size(G,1);
Xmax = M; X_min = 1; dim = Xmax - 2;
max_gen = 100; num_population = 50; num_runs = 5;

Algos = {'EGSO', 'GWO', 'PSO', 'ABC', 'SSA', 'ARO'};
NFE_mean = zeros(numel(Algos),1);
NFE_std  = zeros(numel(Algos),1);

fprintf('========================================\n');
fprintf('NFE comparison, Map 1, N=50, T=100\n');
fprintf('========================================\n\n');

for a = 1:numel(Algos)
    name = Algos{a};
    fn = str2func(name);
    nfes = zeros(num_runs,1);
    for run = 1:num_runs
        counter = FobjCounter(@(x) fitness(x, G));
        fn(num_population, max_gen, X_min, Xmax, dim, @(x) counter.call(x), G);
        nfes(run) = counter.count;
    end
    NFE_mean(a) = mean(nfes);
    NFE_std(a)  = std(nfes);
    fprintf('%-8s: NFE = %8.1f +/- %6.1f  (theoretical N*T=%d for 1-eval/individual/iter)\n', ...
        name, NFE_mean(a), NFE_std(a), num_population*max_gen);
end

save('results/main/NFE_Comparison.mat', 'Algos', 'NFE_mean', 'NFE_std', 'num_population', 'max_gen');

fid = fopen('results/main/NFE_Comparison.txt', 'w');
fprintf(fid, 'NFE (Number of Function Evaluations) Comparison, Map1, N=50, T=100, %d runs each\n', num_runs);
fprintf(fid, '=========================================================================\n\n');
for a = 1:numel(Algos)
    fprintf(fid, '%-8s: NFE = %.1f +/- %.1f\n', Algos{a}, NFE_mean(a), NFE_std(a));
end
fclose(fid);
fprintf('\nSaved: results/main/NFE_Comparison.mat / .txt\n');
