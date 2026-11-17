function ok = CollisionFreeSegment3D(p1, p2, G)
% CollisionFreeSegment3D  Discretized line-of-sight collision check
% between two continuous 3D points against a binary voxel grid.
% Used by the simplified 3D path encoding (straight segments between
% per-column waypoints, checked directly) adopted for the minimal 3D
% proof-of-concept in place of the 2D pipeline's cell-stepping
% generateContinuousRoute.m connectivity-repair algorithm.

    [R, C, H] = size(G);
    dist = norm(p2 - p1);
    nsteps = max(2, ceil(dist / 0.5));
    for t = linspace(0, 1, nsteps)
        p = p1 + t * (p2 - p1);
        r = max(1, min(R, round(p(1))));
        c = max(1, min(C, round(p(2))));
        h = max(1, min(H, round(p(3))));
        if G(r, c, h) == 1
            ok = false;
            return;
        end
    end
    ok = true;
end
