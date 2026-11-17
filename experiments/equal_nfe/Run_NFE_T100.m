function Run_NFE_T100(map_idx, num_runs)
% Run_NFE_T100  True number of objective evaluations per run under the
% manuscript's fixed-iteration protocol (N = 50, T = 100), counting every
% call to fitness() -- including LocalSearch.m's direct calls, which the
% FobjCounter wrapper misses (Table B1 of the paper).
%
% Usage: matlab -batch "Run_NFE_T100(1, 5)"

    global NFEB
    here = fileparts(mfilename('fullpath'));
    project_root = fullfile(here, '..', '..');
    addpath(project_root); setup_egso();
    outdir = fullfile(project_root, 'runs', 'equal_nfe');
    if ~exist(fullfile(outdir, 'logs'), 'dir'), mkdir(fullfile(outdir, 'logs')); end
    cd(here);
    assert(strcmp(fileparts(which('fitness')), here), 'counting fitness.m not active');
    if nargin < 2, num_runs = 5; end

    Algos = {'EGSO','GWO','PSO','GA','DE','SSA','WOA','MFO','ALO','HHO','SCA', ...
             'ABC','ACO','BBO','CS','FA','TLBO','BA','FPA','SMA','INFO','SBOA','ARO'};
    [Maps, ~] = GenerateEGSOMaps();
    close all;
    G = Maps{map_idx};
    M = size(G, 1);
    fobj = @(x) fitness(x, G);

    NFE = zeros(num_runs, numel(Algos));
    for run = 1:num_runs
        for a = 1:numel(Algos)
            rng(5e6 + 1e4 * map_idx + 100 * run + a, 'twister');
            NFEB = struct('active', true, 'count', 0, 'max', Inf, ...
                          'best_f', Inf, 'best_x', [], 'ckpt', Inf, 'trace', []);
            feval(Algos{a}, 50, 100, 1, M, M - 2, fobj, G);
            NFE(run, a) = NFEB.count;
        end
    end
    NFEB = struct('active', false, 'count', 0, 'max', Inf, 'best_f', Inf, 'best_x', [], 'ckpt', 1, 'trace', []);
    save(fullfile(outdir, sprintf('NFE_T100_map%d.mat', map_idx)), 'Algos', 'NFE', 'map_idx', 'num_runs');
    for a = 1:numel(Algos)
        fprintf('Map %d %-5s NFE = %9.1f +/- %7.1f\n', map_idx, Algos{a}, mean(NFE(:, a)), std(NFE(:, a)));
    end
end
