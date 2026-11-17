function [path, path_length] = APF_AStarPath(G, S, E, w_pot)
% APF_AStarPath  APF-weighted A* hybrid: a guided-field
% baseline distinct from both the generic swarm/evolutionary algorithms
% and from EGSO itself, combining a hand-crafted repulsive
% potential field with A*'s graph search -- an "APF-A* hybrid",
% structurally similar to the
% APF+global-search combinations in the literature (e.g. APF hybridized
% with A* under optimized parameters).
%
% Edge cost is inflated by the destination cell's (normalized) repulsive
% APF potential, biasing the search away from obstacle-adjacent cells
% while retaining A*'s optimal-search structure (at the cost of strict
% admissibility, as is standard for potential-weighted A* variants).
%
%   G     : MxM binary occupancy grid
%   S, E  : [row, col] start/goal
%   w_pot : repulsive-cost mixing weight (default 1.0)
%
% Uses the same headline APF parameters (rho0=5, krep=50, katt=1) as the
% local-minima census, for consistency.

if nargin < 4
    w_pot = 1.0;
end

[M, N] = size(G);
U = ArtificialPotentialField(G, E, 1, 50, 5);
U_att = 0.5 * 1 * ((repmat((1:M)', 1, N) - E(1)).^2 + (repmat(1:N, M, 1) - E(2)).^2);
U_rep = U - U_att;
U_rep(~isfinite(U_rep)) = 0;
rep_max = max(U_rep(isfinite(U_rep) & G==0), [], 'all');
if isempty(rep_max) || rep_max <= 0
    rep_max = 1;
end
U_rep_norm = U_rep / rep_max;

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
        break;
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
        edge_cost = step_cost(k) * (1 + w_pot * U_rep_norm(nr, nc));
        tentative_g = gScore(r, c) + edge_cost;
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

path = [E(1), E(2)];
r = E(1); c = E(2);
while ~(r == S(1) && c == S(2))
    pr = parentR(r, c); pc = parentC(r, c);
    path = [pr, pc; path]; %#ok<AGROW>
    r = pr; c = pc;
end

% report TRUE geometric length (not the potential-inflated search cost)
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
