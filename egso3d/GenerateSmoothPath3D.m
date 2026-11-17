function path1 = GenerateSmoothPath3D(path, G)
% GenerateSmoothPath3D  3D line-of-sight shortcutting smoothing, analogous
% to GenerateSmoothPath.m: repeatedly drop a middle waypoint whenever the
% straight segment between its neighbors is collision-free.

    if isempty(path) || size(path,1) < 3
        path1 = path;
        return;
    end

    path1 = path;
    i = 1;
    while i <= size(path1,1) - 2
        if CollisionFreeSegment3D(path1(i,:), path1(i+2,:), G)
            path1(i+1,:) = [];
        else
            i = i + 1;
        end
    end
end
