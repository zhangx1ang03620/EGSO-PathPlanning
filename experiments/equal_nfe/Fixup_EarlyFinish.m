function Fixup_EarlyFinish(threshold)
% Fixup_EarlyFinish  Enforce full budget use in the equal-NFE experiment.
%
% Some algorithms (notably FA) perform a data-dependent number of
% objective evaluations per iteration, so the iteration count derived from
% their average evaluations per iteration at T = 100 can let a run finish
% its schedule well before MaxFEs. For every run that ended without
% reaching the budget (not truncated and NFE used < threshold*MaxFEs), this
% script re-runs the same algorithm with the same random seed and a
% proportionally larger T, up to 5 attempts, until the budget is consumed,
% and overwrites that run's entries. Every change is logged in the file's
% 'Fixups' field and in logs/fixup.log.

    global NFEB
    if nargin < 1, threshold = 0.98; end
    here = fileparts(mfilename('fullpath'));
    project_root = fullfile(here, '..', '..');
    addpath(project_root); setup_egso();
    outdir = fullfile(project_root, 'runs', 'equal_nfe');
    cd(here);
    assert(strcmp(fileparts(which('fitness')), here), 'counting fitness.m not active');
    [Maps, ~] = GenerateEGSOMaps();
    close all;
    flog = fopen(fullfile(outdir, 'logs', 'fixup.log'), 'a');
    files = dir(fullfile(outdir, 'EqNFE_map*_fe*_r*.mat'));
    for fi = 1:numel(files)
        fn = fullfile(files(fi).folder, files(fi).name);
        S = load(fn);
        if ~isfield(S, 'Fixups'), S.Fixups = {}; end
        G = Maps{S.map_idx}; M = size(G, 1); dim = M - 2;
        fobj = @(x) fitness(x, G);
        changed = false;
        for ri = 1:numel(S.runs)
            for a = 1:numel(S.Algos)
                if S.Truncated(ri, a) || S.NFEused(ri, a) >= threshold * S.maxfes || ~isempty(S.ErrMsg{ri, a})
                    continue
                end
                name = S.Algos{a};
                seed = 1e6 * (S.maxfes >= 10000) + 1e4 * S.map_idx + 100 * S.runs(ri) + a;
                T = S.Tused(ri, a); used = S.NFEused(ri, a);
                for attempt = 1:5
                    T = ceil(T * S.maxfes / max(used, 1) * 1.1);
                    rng(seed, 'twister');
                    nfeb_reset(S.maxfes, S.ckpt);
                    t0 = tic; trunc = false;
                    try
                        feval(name, S.N, T, 1, M, dim, fobj, G);
                    catch ME
                        if strcmp(ME.identifier, 'NFEB:exhausted'), trunc = true; else, rethrow(ME); end
                    end
                    used = NFEB.count;
                    if trunc || used >= threshold * S.maxfes, break; end
                end
                tr = NFEB.trace; last = NaN;
                for k = 1:numel(tr)
                    if isnan(tr(k)), tr(k) = last; else, last = tr(k); end
                end
                if used < S.maxfes, tr(floor(used / S.ckpt) + 1:end) = NFEB.best_f; end
                msg = sprintf('%s run %d %s: NFE %d -> %d (T %d -> %d, attempts %d), best %.2f -> %.2f', ...
                    files(fi).name, S.runs(ri), name, S.NFEused(ri, a), used, S.Tused(ri, a), T, attempt, S.Best(ri, a), NFEB.best_f);
                fprintf('%s\n', msg); fprintf(flog, '%s\n', msg);
                S.Fixups{end + 1} = msg;
                S.Best(ri, a) = NFEB.best_f; S.BestX{ri, a} = NFEB.best_x;
                S.NFEused(ri, a) = used; S.Tused(ri, a) = T; S.Truncated(ri, a) = trunc;
                S.Time(ri, a) = toc(t0); S.Trace(ri, a, :) = tr;
                nfeb_reset();
                changed = true;
            end
        end
        if changed
            save(fn, '-struct', 'S');
        end
    end
    fclose(flog);
    fprintf('fixup done\n');
end
