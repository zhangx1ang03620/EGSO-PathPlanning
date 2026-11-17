function fx = fitness(x, G)
% fitness  Budget-enforcing, evaluation-counting wrapper.
%
% This file deliberately shadows the project-root fitness.m: MATLAB gives
% the current folder precedence over the path, so when an experiment is
% run from this folder EVERY objective evaluation goes through here --
% including the direct fitness(x,G) calls inside LocalSearch.m that the
% FobjCounter wrapper could not see. The objective itself is the
% verbatim original (fitness_orig.m); no algorithm file is modified.
%
% State lives in the global struct NFEB (see nfeb_reset.m):
%   active  - counting/budget enforcement on/off
%   count   - evaluations performed so far in the current run
%   max     - MaxFEs budget; the (max+1)-th call throws 'NFEB:exhausted'
%   best_f  - best (lowest) fitness evaluated so far (anytime best)
%   best_x  - decision vector of best_f
%   ckpt    - checkpoint spacing for the best-so-far trace
%   trace   - best-so-far fitness at NFE = ckpt, 2*ckpt, ..., max

    global NFEB
    if ~isempty(NFEB) && NFEB.active
        if NFEB.count >= NFEB.max
            error('NFEB:exhausted', 'NFE budget of %d evaluations exhausted', NFEB.max);
        end
        NFEB.count = NFEB.count + 1;
        fx = fitness_orig(x, G);
        if fx < NFEB.best_f
            NFEB.best_f = fx;
            NFEB.best_x = x;
        end
        if mod(NFEB.count, NFEB.ckpt) == 0
            NFEB.trace(NFEB.count / NFEB.ckpt) = NFEB.best_f;
        end
    else
        fx = fitness_orig(x, G);
    end
end
