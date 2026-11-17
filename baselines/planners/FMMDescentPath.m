function [path, path_length] = FMMDescentPath(G, S, E)
% FMMDescentPath  Extracts a path by pure steepest descent on the
% precomputed Eikonal viscosity field T_nav (4-connected, the SAME field
% EGSO uses for geodesic seeding) --
% with NO population, no swarm refinement, no CALS smoothing.
%
% This is the "EGSO minus everything except the raw guidance field"
% reference: it answers the question "what does
% the swarm/refinement layer actually add over the guidance field alone?"
% Deterministic (not a stochastic metaheuristic).
%
%   G : MxM binary occupancy grid (1 = obstacle, 0 = free)
%   S, E : [row, col] start / goal cells
%
% Returns path (Nx2 [row col]) and its Euclidean length (Inf if no path).

M = size(G, 1);
goalMask = false(size(G));
goalMask(E(1), E(2)) = true;
T = FastMarchingField(G, goalMask, 'navigate');

if isinf(T(S(1), S(2)))
    path = [];
    path_length = Inf;
    return;
end

nbrs4 = [-1 0; 1 0; 0 -1; 0 1];

path = [S(1), S(2)];
r = S(1); c = S(2);
visited = false(M, M);
visited(r, c) = true;

maxSteps = M * M * 2;
found = false;
for step = 1:maxSteps
    if r == E(1) && c == E(2)
        found = true;
        break;
    end

    best_val = Inf; best_r = -1; best_c = -1;
    for k = 1:4
        nr = r + nbrs4(k, 1);
        nc = c + nbrs4(k, 2);
        if nr < 1 || nr > M || nc < 1 || nc > M
            continue;
        end
        if visited(nr, nc)
            continue;
        end
        tval = T(nr, nc);
        if tval < best_val
            best_val = tval; best_r = nr; best_c = nc;
        end
    end

    if best_r < 0
        break;   % no unvisited downhill neighbor -- stuck, fail
    end

    r = best_r; c = best_c;
    visited(r, c) = true;
    path = [path; r, c]; %#ok<AGROW>
end

if ~found || ~(r == E(1) && c == E(2))
    path = [];
    path_length = Inf;
    return;
end

path_length = 0;
for i = 1:size(path, 1) - 1
    path_length = path_length + norm(path(i+1, :) - path(i, :));
end

end
