function [trap_rate, num_sinks, reach_mask, sink_mask] = GreedyDescentCensus(T, G, E)
% GreedyDescentCensus  Exhaustive census of spurious local minima in a
%                      scalar guidance field T over a binary occupancy
%                      grid G, used to turn the paper's theoretical claim
%                      ("the Eikonal viscosity field is structurally free
%                      of spurious local minima away from the goal, unlike
%                      a hand-crafted APF") into a directly measured,
%                      reproducible statistic instead of a citation to
%                      Sethian's consistency theorem alone.
%
%   [trap_rate, num_sinks, reach_mask, sink_mask] = GreedyDescentCensus(T, G, E)
%
%   T : MxM scalar field (lower = more attractive), Inf/NaN at cells to
%       exclude from the traversable domain (obstacles, unreachable).
%   G : MxM binary occupancy grid (1 = obstacle, 0 = free)
%   E : [row, col] goal cell (the only point exempt from being a "sink")
%
%   For every free cell, performs 4-connected steepest descent (always
%   moving to the strictly-lowest-valued free neighbor; ties broken by
%   first match) until either the goal is reached or no strictly-lower
%   free neighbor exists (a "sink").
%
%   trap_rate : percent of free reachable-from-goal cells whose descent
%               terminates at a non-goal sink (does NOT reach E)
%   num_sinks : count of distinct non-goal cells that are themselves local
%               minima (no free neighbor has a strictly lower T)
%   reach_mask: MxM logical, true where descent from that cell reaches E
%   sink_mask : MxM logical, true at non-goal local-minima cells

    [M, N] = size(G);
    free_mask = (G == 0) & isfinite(T);
    reachable = free_mask;   % cells the field itself considers reachable

    nbrs4 = [-1 0; 1 0; 0 -1; 0 1];

    % ---- 1. exhaustive local-minima ("sink") census: no descent needed ----
    sink_mask = false(M, N);
    [rIdx, cIdx] = find(free_mask);
    for k = 1:numel(rIdx)
        r = rIdx(k); c = cIdx(k);
        if r == E(1) && c == E(2)
            continue;
        end
        here = T(r, c);
        is_sink = true;
        for nb = 1:4
            nr = r + nbrs4(nb,1); nc = c + nbrs4(nb,2);
            if nr < 1 || nr > M || nc < 1 || nc > N, continue; end
            if ~free_mask(nr, nc), continue; end
            if T(nr, nc) < here - 1e-9
                is_sink = false;
                break;
            end
        end
        if is_sink
            sink_mask(r, c) = true;
        end
    end
    num_sinks = sum(sink_mask(:));

    % ---- 2. per-cell greedy-descent reachability ----
    reach_mask = false(M, N);
    reach_mask(E(1), E(2)) = true;
    maxSteps = M * N;

    for k = 1:numel(rIdx)
        r0 = rIdx(k); c0 = cIdx(k);
        if reach_mask(r0, c0), continue; end

        path_r = zeros(maxSteps, 1); path_c = zeros(maxSteps, 1);
        r = r0; c = c0; nstep = 0;
        visited_local = false(M, N);
        success = false;
        for step = 1:maxSteps
            if reach_mask(r, c)
                success = true;
                break;
            end
            if r == E(1) && c == E(2)
                success = true;
                break;
            end
            visited_local(r, c) = true;
            nstep = nstep + 1;
            path_r(nstep) = r; path_c(nstep) = c;

            here = T(r, c);
            best_val = here; best_r = -1; best_c = -1;
            for nb = 1:4
                nr = r + nbrs4(nb,1); nc = c + nbrs4(nb,2);
                if nr < 1 || nr > M || nc < 1 || nc > N, continue; end
                if ~free_mask(nr, nc) || visited_local(nr, nc), continue; end
                if T(nr, nc) < best_val - 1e-9
                    best_val = T(nr, nc); best_r = nr; best_c = nc;
                end
            end
            if best_r < 0
                break;   % stuck: no unvisited strictly-lower neighbor
            end
            r = best_r; c = best_c;
        end

        if success
            reach_mask(E(1), E(2)) = true;
        end
        for s = 1:nstep
            reach_mask(path_r(s), path_c(s)) = success;
        end
    end

    total_free = sum(reachable(:));
    trapped = total_free - sum(reach_mask(:) & reachable(:));
    trap_rate = 100 * trapped / total_free;
end
