function Run_EqualNFE_Job(map_idx, maxfes, run_from, run_to)
% Run_EqualNFE_Job  Equal-NFE-budget comparison.
%
% Runs all 23 algorithms on one benchmark map under a common, hard
% maximum number of objective-function evaluations (MaxFEs), for runs
% run_from..run_to, and saves one result file. Launch many of these as
% independent `matlab -batch` processes (see run_all_jobs.sh); aggregate
% with Aggregate_EqualNFE.m.
%
% Protocol
%   * Every objective evaluation is counted by the shadowing fitness.m in
%     this folder (including LocalSearch.m's direct calls), and the
%     (MaxFEs+1)-th call aborts the run with error 'NFEB:exhausted'.
%   * The reported result of a run is the best fitness evaluated within
%     the budget (anytime best), recorded identically for all algorithms.
%   * Each algorithm's iteration count T is set from its measured
%     evaluations per iteration at T = 100 (calibration/; EGSO uses
%     the revised implementation's mean count), with a
%     5% margin, so that its own T-dependent schedules (coefficient decay,
%     restart windows, EGSO's CALS onset at 0.65T, ...) run essentially to
%     completion when the budget is reached. N = 50 for every algorithm.
%
% Usage: matlab -batch "Run_EqualNFE_Job(5, 15000, 1, 10)"

    global NFEB
    here = fileparts(mfilename('fullpath'));
    project_root = fullfile(here, '..', '..');
    addpath(project_root); setup_egso();
    outdir = fullfile(project_root, 'runs', 'equal_nfe');
    if ~exist(fullfile(outdir, 'logs'), 'dir'), mkdir(fullfile(outdir, 'logs')); end
    cd(here);
    w = which('fitness');
    assert(strcmp(fileparts(w), here), 'fitness.m must resolve to the budget wrapper in %s, got %s', here, w);

    Algos = {'EGSO','GWO','PSO','GA','DE','SSA','WOA','MFO','ALO','HHO','SCA', ...
             'ABC','ACO','BBO','CS','FA','TLBO','BA','FPA','SMA','INFO','SBOA','ARO'};
    nA = numel(Algos);
    N = 50;
    ckpt = 50;
    margin = 1.05;

    [Maps, ~] = GenerateEGSOMaps();
    close all;
    G = Maps{map_idx};
    M = size(G, 1);
    dim = M - 2;
    fobj = @(x) fitness(x, G);

    nfe100 = load_nfe100(fullfile(here, 'calibration'), map_idx, Algos);

    runs = run_from:run_to;
    nR = numel(runs);
    Best    = nan(nR, nA);
    NFEused = zeros(nR, nA);
    Tused   = zeros(nR, nA);
    Truncated = false(nR, nA);
    Time    = zeros(nR, nA);
    Trace   = nan(nR, nA, floor(maxfes / ckpt));
    BestX   = cell(nR, nA);
    ErrMsg  = cell(nR, nA);

    for ri = 1:nR
        run = runs(ri);
        for a = 1:nA
            name = Algos{a};
            init = 50;
            if strcmp(name, 'EGSO'), init = 10 * N; end
            per_iter = (nfe100(a) - init) / 100;
            T = max(2, ceil(margin * (maxfes - init) / per_iter));

            seed = 1e6 * (maxfes >= 10000) + 1e4 * map_idx + 100 * run + a;
            rng(seed, 'twister');
            nfeb_reset(maxfes, ckpt);
            t0 = tic;
            try
                feval(name, N, T, 1, M, dim, fobj, G);
            catch ME
                if strcmp(ME.identifier, 'NFEB:exhausted')
                    Truncated(ri, a) = true;
                else
                    ErrMsg{ri, a} = ME.message;
                    fprintf(2, '[Map %d FE %d run %d] %s ERROR: %s\n', map_idx, maxfes, run, name, ME.message);
                end
            end
            Time(ri, a) = toc(t0);

            tr = NFEB.trace;
            last = NaN;
            for k = 1:numel(tr)            % forward-fill after an early finish
                if isnan(tr(k)), tr(k) = last; else, last = tr(k); end
            end
            if NFEB.count < maxfes          % finished before the budget
                tr(floor(NFEB.count / ckpt) + 1:end) = NFEB.best_f;
            end
            Best(ri, a)    = NFEB.best_f;
            BestX{ri, a}   = NFEB.best_x;
            NFEused(ri, a) = NFEB.count;
            Tused(ri, a)   = T;
            Trace(ri, a, :) = tr;
            nfeb_reset();
        end
        fprintf('[Map %d | MaxFEs %d] run %d done (EGSO %.2f)\n', map_idx, maxfes, run, Best(ri, 1));
    end

    out = fullfile(outdir, sprintf('EqNFE_map%d_fe%d_r%02d-%02d.mat', map_idx, maxfes, run_from, run_to));
    save(out, 'Algos', 'map_idx', 'maxfes', 'runs', 'N', 'ckpt', 'margin', 'nfe100', ...
         'Best', 'NFEused', 'Tused', 'Truncated', 'Time', 'Trace', 'BestX', 'ErrMsg');
    fprintf('saved %s\n', out);
end

function nfe100 = load_nfe100(folder, map_idx, Algos)
% Measured evaluations per run at T = 100, N = 50 (see calibration/).
    txt = fileread(fullfile(folder, sprintf('nfe_t100_map%d.txt', map_idx)));
    lines = strsplit(strtrim(txt), newline);
    nfe100 = nan(1, numel(Algos));
    for i = 1:numel(lines)
        f = strsplit(lines{i}, sprintf('\t'));
        k = find(strcmp(Algos, f{2}));
        nfe100(k) = str2double(extractAfter(f{3}, 'all='));
    end
    assert(all(isfinite(nfe100)), 'missing calibration for some algorithm');
end
