function T = ChamferDistanceField(G, goalMask)
% ChamferDistanceField  Classic two-pass (forward+backward raster scan)
% chamfer-3-4 distance transform, restricted to free cells (obstacle
% neighbors never contribute to propagation, so it does not "cut through
% walls" the way an unrestricted Euclidean/chamfer transform would).
%
% Answers the question "why not use the cheaper chamfer distance
% transform instead of FMM" with an actual chamfer implementation, not only the Dijkstra/BFS proxy already used in the
% factorial ablation. Unlike BFS/Dijkstra (a proper queue-based
% wavefront that always propagates in the true order of increasing
% distance, exactly like FMM's frozen-cell ordering), a chamfer transform
% only sweeps the grid in TWO fixed raster passes. Around concave,
% multi-turn obstacle geometry (this paper's engineered decoy channels),
% a single forward+backward sweep can fail to propagate the wavefront in
% the correct order, giving locally inconsistent (non-monotone) distance
% values -- this is a known, textbook limitation of chamfer transforms
% for geodesic distance approximation in non-convex domains, and is
% exactly the failure mode this function is built to exhibit and measure
% (via the same GreedyDescentCensus.m used by the local-minima census).
%
%   G        : MxM binary occupancy grid (1 = obstacle, 0 = free)
%   goalMask : MxM logical mask, true at the source (goal) cell(s)
%
% Returns T (MxM), Inf at obstacle cells and at free cells never reached
% by either pass (fully disconnected from the source under this
% propagation order).
%
% Chamfer-3-4 weights (axis-neighbor cost 3, diagonal-neighbor cost 2*3/2
% -> using integer chamfer-3-4: axis=3, diagonal=4; approximates true
% Euclidean distance scaled by ~3). Only used for descent-order/monotonicity
% testing (GreedyDescentCensus only checks relative ordering), so the
% absolute scale is immaterial.

[M, N] = size(G);
free = (G == 0);
T = inf(M, N);
T(goalMask) = 0;

AXIS = 3; DIAG = 4;

% ---- forward pass: top-left -> bottom-right ----
for r = 1:M
    for c = 1:N
        if ~free(r, c)
            continue;
        end
        best = T(r, c);
        if r > 1 && c > 1 && free(r-1, c-1), best = min(best, T(r-1, c-1) + DIAG); end
        if r > 1,                            best = min(best, T(r-1, c)   + AXIS); end
        if r > 1 && c < N && free(r-1, c+1), best = min(best, T(r-1, c+1) + DIAG); end
        if c > 1 && free(r, c-1),            best = min(best, T(r, c-1)   + AXIS); end
        T(r, c) = best;
    end
end

% ---- backward pass: bottom-right -> top-left ----
for r = M:-1:1
    for c = N:-1:1
        if ~free(r, c)
            continue;
        end
        best = T(r, c);
        if r < M && c < N && free(r+1, c+1), best = min(best, T(r+1, c+1) + DIAG); end
        if r < M,                            best = min(best, T(r+1, c)   + AXIS); end
        if r < M && c > 1 && free(r+1, c-1), best = min(best, T(r+1, c-1) + DIAG); end
        if c < N && free(r, c+1),            best = min(best, T(r, c+1)   + AXIS); end
        T(r, c) = best;
    end
end

T(~free) = inf;

end
