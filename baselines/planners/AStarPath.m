function [path, path_length] = AStarPath(G, S, E)
% AStarPath  Standard 8-connected A* search on a binary occupancy grid.
%
% Deterministic classical reference baseline (not a stochastic
% metaheuristic): reports the true (near-)optimal 2D shortest path,
% unconstrained by the "one row per column" encoding used by every
% swarm-based comparison algorithm in this project. Included specifically
% to answer the question "how close is EGSO to the
% classical optimum, and why not just use A*?".
%
%   G : MxM binary occupancy grid (1 = obstacle, 0 = free)
%   S, E : [row, col] start / goal cells
%
% Returns path (Nx2 [row col]) and its Euclidean length (Inf if no path).

[M, N] = size(G);

gScore  = inf(M, N);
fScore  = inf(M, N);
closed  = false(M, N);
inOpen  = false(M, N);
parentR = zeros(M, N);
parentC = zeros(M, N);

gScore(S(1), S(2)) = 0;
fScore(S(1), S(2)) = octileHeuristic(S, E);
inOpen(S(1), S(2)) = true;

nbrs8 = [-1 -1; -1 0; -1 1; 0 -1; 0 1; 1 -1; 1 0; 1 1];
step_cost = [sqrt(2); 1; sqrt(2); 1; 1; sqrt(2); 1; sqrt(2)];

found = false;
maxIter = M * N * 2;
for iter = 1:maxIter
    masked = fScore;
    masked(~inOpen) = Inf;
    [minVal, idx] = min(masked(:));
    if isinf(minVal)
        break;   % open set exhausted, no path
    end
    [r, c] = ind2sub([M, N], idx);
    inOpen(r, c) = false;

    if r == E(1) && c == E(2)
        found = true;
        break;
    end
    closed(r, c) = true;

    for k = 1:8
        nr = r + nbrs8(k, 1);
        nc = c + nbrs8(k, 2);
        if nr < 1 || nr > M || nc < 1 || nc > N
            continue;
        end
        if G(nr, nc) == 1 || closed(nr, nc)
            continue;
        end
        tentative_g = gScore(r, c) + step_cost(k);
        if tentative_g < gScore(nr, nc)
            gScore(nr, nc)  = tentative_g;
            fScore(nr, nc)  = tentative_g + octileHeuristic([nr, nc], E);
            parentR(nr, nc) = r;
            parentC(nr, nc) = c;
            inOpen(nr, nc)  = true;
        end
    end
end

if ~found
    path = [];
    path_length = Inf;
    return;
end

% reconstruct path by walking parent pointers back from E to S
path = [E(1), E(2)];
r = E(1); c = E(2);
while ~(r == S(1) && c == S(2))
    pr = parentR(r, c); pc = parentC(r, c);
    path = [pr, pc; path]; %#ok<AGROW>
    r = pr; c = pc;
end

path_length = 0;
for i = 1:size(path, 1) - 1
    path_length = path_length + norm(path(i+1, :) - path(i, :));
end

end

function h = octileHeuristic(a, b)
    dr = abs(a(1) - b(1));
    dc = abs(a(2) - b(2));
    h = max(dr, dc) + (sqrt(2) - 1) * min(dr, dc);
end
