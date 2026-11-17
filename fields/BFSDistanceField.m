function T = BFSDistanceField(G, goalMask)
% BFSDistanceField  4-connected breadth-first-search hop-count field from
% the cell(s) marked true in goalMask, over the free cells (G==0) of a
% binary occupancy grid.
%
% This is a discrete, globally monotone "shortest-hop-to-goal" field:
% like the Eikonal viscosity field T_nav, every goal-reachable non-goal
% free cell has at least one free neighbor with a strictly smaller value
% (its BFS predecessor), so it shares T_nav's global-connectivity /
% local-minima-free property by construction. It differs from T_nav in
% that it is a purely discrete unweighted-graph distance (integer hop
% count, no continuous viscosity-solution structure, no sub-cell
% precision, no natural extension to arbitrary-angle or 3D voxel
% movement).
%
% Used as the third arm ("Dijkstra/BFS field") of the factorial ablation
% in the factorial ablation, alongside the Eikonal field and a hand-crafted APF, to
% separate "the value of using ANY globally monotone guidance field" from
% "the specific added value of the Eikonal/viscosity-solution field".
%
%   G        : MxM binary occupancy grid (1 = obstacle, 0 = free)
%   goalMask : MxM logical mask, true at the source cell(s)
%
% Returns T (MxM), Inf at obstacle cells and at free cells not reachable
% from the goal.

[M, N] = size(G);
T = inf(M, N);

[gr, gc] = find(goalMask);
num_src = numel(gr);
q = zeros(M*N, 2);
head = 1; tail = 0;

for k = 1:num_src
    if isinf(T(gr(k), gc(k)))
        T(gr(k), gc(k)) = 0;
        tail = tail + 1;
        q(tail, :) = [gr(k), gc(k)];
    end
end

nbrs4 = [-1 0; 1 0; 0 -1; 0 1];
while head <= tail
    r = q(head, 1); c = q(head, 2);
    head = head + 1;
    d = T(r, c);
    for k = 1:4
        nr = r + nbrs4(k, 1);
        nc = c + nbrs4(k, 2);
        if nr >= 1 && nr <= M && nc >= 1 && nc <= N && G(nr, nc) == 0 && isinf(T(nr, nc))
            T(nr, nc) = d + 1;
            tail = tail + 1;
            q(tail, :) = [nr, nc];
        end
    end
end

end
