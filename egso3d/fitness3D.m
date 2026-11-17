function fx = fitness3D(x, G)
% fitness3D  Minimal 3D analogue of fitness.m. Decision vector x has
% length 2*dim (dim = C-2): x(1:dim) are row indices, x(dim+1:2*dim) are
% height indices, one pair per intermediate column. Start S=(1,1,1),
% goal E=(R,C,H).
%
% Simplified path model (see CollisionFreeSegment3D.m docstring): waypoints
% are connected by straight 3D segments checked directly against the
% voxel grid, rather than porting the 2D pipeline's cell-stepping
% connectivity-repair algorithm to 3D.

    [R, C, H] = size(G);
    S = [1, 1, 1];
    E = [R, C, H];
    dim = length(x) / 2;

    row = round(x(1:dim));
    hgt = round(x(dim+1:end));
    row = max(1, min(R, row));
    hgt = max(1, min(H, hgt));

    waypoints = [S; [row(:), (2:dim+1)', hgt(:)]; E];

    nB = 0;
    for j = 1:size(waypoints,1)
        p = waypoints(j,:);
        if G(p(1), p(2), p(3)) == 1
            nB = nB + 1;
        end
    end

    if nB > 0
        fx = 1e6 * nB;
        return;
    end

    for j = 1:size(waypoints,1)-1
        if ~CollisionFreeSegment3D(waypoints(j,:), waypoints(j+1,:), G)
            nB = nB + 1;
        end
    end

    if nB > 0
        fx = 1e6 * nB;
        return;
    end

    path = GenerateSmoothPath3D(waypoints, G);
    fx = 0;
    for i = 1:size(path,1)-1
        fx = fx + norm(path(i+1,:) - path(i,:));
    end
end
