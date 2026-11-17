function T = FastMarchingField3D(G, sourceMask, domainMode)
% FastMarchingField3D  3D extension (6-connectivity) of FastMarchingField.m,
% used as a minimal proof-of-concept that the EVF module's viscosity-solution
% guidance is not restricted to 2D grids and can be transferred to a
% voxelized 3D feasible region (Discussion Sec. "Application Prospects",
% previously an unverified future-work claim).
%
%   T = FastMarchingField3D(G, sourceMask, domainMode)
%
%   G          : RxCxH binary occupancy voxel grid (1 = obstacle, 0 = free)
%   sourceMask : RxCxH logical, true at cells where T = 0
%   domainMode : 'navigate'  - obstacles impassable (speed F=0)
%                'clearance' - all voxels traversable (F=1 everywhere)
%
%   Returns T (RxCxH double), Inf where unreachable.

    [R, C, H] = size(G);
    switch domainMode
        case 'navigate'
            speed = double(~G);
        case 'clearance'
            speed = ones(R, C, H);
        otherwise
            error('FastMarchingField3D:BadMode', 'domainMode must be ''navigate'' or ''clearance''.');
    end

    T      = inf(R, C, H);
    Trial  = inf(R, C, H);
    frozen = false(R, C, H);

    T(sourceMask)      = 0;
    frozen(sourceMask) = true;

    [sr, sc, sh] = ind2sub([R,C,H], find(sourceMask));
    for k = 1:numel(sr)
        [Trial, T, frozen] = touchNeighbors3D(sr(k), sc(k), sh(k), T, Trial, frozen, speed, R, C, H);
    end

    numCells = R * C * H;
    for iter = 1:numCells
        maskedTrial = Trial;
        maskedTrial(frozen) = Inf;
        [minVal, idx] = min(maskedTrial(:));
        if isinf(minVal)
            break;
        end
        [ri, ci, hi] = ind2sub([R,C,H], idx);
        frozen(ri, ci, hi) = true;
        T(ri, ci, hi)      = minVal;
        [Trial, T, frozen] = touchNeighbors3D(ri, ci, hi, T, Trial, frozen, speed, R, C, H);
    end
end

%% ========================================================================
function [Trial, T, frozen] = touchNeighbors3D(r, c, h, T, Trial, frozen, speed, R, C, H)
    nbrs = [r-1,c,h; r+1,c,h; r,c-1,h; r,c+1,h; r,c,h-1; r,c,h+1];
    for k = 1:6
        nr = nbrs(k,1); nc = nbrs(k,2); nh = nbrs(k,3);
        if nr < 1 || nr > R || nc < 1 || nc > C || nh < 1 || nh > H
            continue;
        end
        if frozen(nr, nc, nh) || speed(nr, nc, nh) <= 0
            continue;
        end
        cand = solveQuadratic3D(nr, nc, nh, T, frozen, speed, R, C, H);
        if cand < Trial(nr, nc, nh)
            Trial(nr, nc, nh) = cand;
        end
    end
end

%% ========================================================================
%  First-order upwind 3D Eikonal update: sort the up-to-3 finite frozen
%  upwind neighbor values and try the highest-dimensional consistent
%  solve first, falling back to fewer dimensions (standard multi-D FMM
%  update, e.g. Sethian 1999 Sec. 4.2).
%% ========================================================================
function cand = solveQuadratic3D(r, c, h, T, frozen, speed, R, C, H)
    vals = inf(1,3);
    if r > 1 && frozen(r-1,c,h), vals(1) = min(vals(1), T(r-1,c,h)); end
    if r < R && frozen(r+1,c,h), vals(1) = min(vals(1), T(r+1,c,h)); end
    if c > 1 && frozen(r,c-1,h), vals(2) = min(vals(2), T(r,c-1,h)); end
    if c < C && frozen(r,c+1,h), vals(2) = min(vals(2), T(r,c+1,h)); end
    if h > 1 && frozen(r,c,h-1), vals(3) = min(vals(3), T(r,c,h-1)); end
    if h < H && frozen(r,c,h+1), vals(3) = min(vals(3), T(r,c,h+1)); end

    cost = 1 / speed(r,c,h);
    finiteVals = sort(vals(isfinite(vals)));
    nAvail = numel(finiteVals);

    if nAvail == 0
        cand = inf;
        return;
    end

    for k = nAvail:-1:1
        v = finiteVals(1:k);
        sumV  = sum(v);
        sumV2 = sum(v.^2);
        a = k; b = -2*sumV; cc = sumV2 - cost^2;
        disc = b^2 - 4*a*cc;
        if disc < 0
            continue;
        end
        candTry = (-b + sqrt(disc)) / (2*a);
        if k == 1 || candTry >= v(k) - 1e-9
            cand = candTry;
            return;
        end
    end
    cand = min(finiteVals) + cost;
end
