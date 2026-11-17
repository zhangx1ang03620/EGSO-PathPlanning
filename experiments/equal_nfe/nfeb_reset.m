function nfeb_reset(maxfes, ckpt)
% nfeb_reset  Arm the global evaluation counter/budget used by fitness.m
% for one algorithm run. Call with no arguments to disable counting.

    global NFEB
    if nargin == 0
        NFEB = struct('active', false, 'count', 0, 'max', Inf, ...
                      'best_f', Inf, 'best_x', [], 'ckpt', 1, 'trace', []);
        return
    end
    NFEB = struct('active', true, 'count', 0, 'max', maxfes, ...
                  'best_f', Inf, 'best_x', [], 'ckpt', ckpt, ...
                  'trace', nan(1, floor(maxfes / ckpt)));
end
