function T = FastMarchingField(G, sourceMask, domainMode)
% FastMarchingField  Solve the Eikonal equation |grad T| = 1/F on a grid map
%                     via the first-order Fast Marching Method (Sethian, 1996).
%                     T is the (unique) viscosity solution of the Eikonal PDE:
%                     monotone and consistent, hence free of spurious local
%                     minima away from the source set (Sethian's consistency
%                     theorem for the upwind Godunov scheme).
%
%   T = FastMarchingField(G, sourceMask, domainMode)
%
%   G           : MxN binary occupancy grid (1 = obstacle, 0 = free)
%   sourceMask  : MxN logical, true at cells where T = 0 (wavefront origin)
%   domainMode  : 'navigate'  - obstacles are impassable (speed F=0), used to
%                                build the goal-arrival-time guidance field
%                 'clearance' - all cells traversable (F=1 everywhere), used
%                                to build the obstacle-distance/clearance field
%
%   Returns T (MxN double), Inf where unreachable from the source set.
%
%   Precomputed once per map (O((M*N)^2) worst case via full-matrix argmin,
%   negligible cost since it never re-evaluates the black-box fitness fobj).

    [M, N] = size(G);
    switch domainMode
        case 'navigate'
            speed = double(~G);          % 0 in obstacles -> impassable
        case 'clearance'
            speed = ones(M, N);          % uniform speed; obstacles act only as sources
        otherwise
            error('FastMarchingField:BadMode', 'domainMode must be ''navigate'' or ''clearance''.');
    end

    T      = inf(M, N);
    Trial  = inf(M, N);
    frozen = false(M, N);

    T(sourceMask)      = 0;
    frozen(sourceMask) = true;

    [sr, sc] = find(sourceMask);
    for k = 1:numel(sr)
        [Trial, T, frozen] = touchNeighbors(sr(k), sc(k), T, Trial, frozen, speed, M, N);
    end

    numCells = M * N;
    for iter = 1:numCells
        maskedTrial = Trial;
        maskedTrial(frozen) = Inf;
        [minVal, idx] = min(maskedTrial(:));
        if isinf(minVal)
            break;   % no more reachable trial cells
        end
        [ri, ci] = ind2sub([M, N], idx);
        frozen(ri, ci) = true;
        T(ri, ci)      = minVal;
        [Trial, T, frozen] = touchNeighbors(ri, ci, T, Trial, frozen, speed, M, N);
    end
end

%% ========================================================================
function [Trial, T, frozen] = touchNeighbors(r, c, T, Trial, frozen, speed, M, N)
    nbrs = [r-1, c; r+1, c; r, c-1; r, c+1];
    for k = 1:4
        nr = nbrs(k,1); nc = nbrs(k,2);
        if nr < 1 || nr > M || nc < 1 || nc > N
            continue;
        end
        if frozen(nr, nc) || speed(nr, nc) <= 0
            continue;
        end
        cand = solveQuadratic(nr, nc, T, frozen, speed, M, N);
        if cand < Trial(nr, nc)
            Trial(nr, nc) = cand;
        end
    end
end

%% ========================================================================
%  First-order upwind Eikonal update (Sethian 1996, Eq. for 2D grid):
%  solves (T-Tx)^2 + (T-Ty)^2 = cost^2 for the smaller upwind neighbors
%% ========================================================================
function cand = solveQuadratic(r, c, T, frozen, speed, M, N)
    Tx = inf; Ty = inf;
    if r > 1 && frozen(r-1, c), Tx = min(Tx, T(r-1, c)); end
    if r < M && frozen(r+1, c), Tx = min(Tx, T(r+1, c)); end
    if c > 1 && frozen(r, c-1), Ty = min(Ty, T(r, c-1)); end
    if c < N && frozen(r, c+1), Ty = min(Ty, T(r, c+1)); end

    cost = 1 / speed(r, c);

    if isinf(Tx) && isinf(Ty)
        cand = inf;
    elseif isinf(Ty)
        cand = Tx + cost;
    elseif isinf(Tx)
        cand = Ty + cost;
    else
        if abs(Tx - Ty) >= cost
            cand = min(Tx, Ty) + cost;
        else
            disc = 2*cost^2 - (Tx - Ty)^2;
            disc = max(disc, 0);
            cand = (Tx + Ty + sqrt(disc)) / 2;
        end
    end
end
