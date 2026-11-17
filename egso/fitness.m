function fx = fitness(x, G)
    S = [1 1];       % start
    E = size(G);     % goal
    route = [S(1) x E(1)];
    dim = length(route);
    nB = 0;        % collision count
    route = round(route);
    
    % Boundary check
    route = max(route, 1);
    route = min(route, size(G, 1));
    
    for j = 2 : dim-1
       if route(j) < 1 || route(j) > size(G, 1) || j > size(G, 2)
           nB = nB + 1;
       elseif G(route(j), j) == 1
           nB = nB + 1;
       end
    end
    
    if nB == 0     % collision-free: compute the path length
        path = generateContinuousRoute(route, G);
        
        if isempty(path)
            fx = 1e6;
            return;
        end
        
        path = GenerateSmoothPath(path, G);
        path = GenerateSmoothPath(path, G);
        
        if isempty(path)
            fx = 1e6;
            return;
        end
        
        fx = 0;
        for i = 1:size(path, 1)-1
            fx = fx + sqrt((path(i+1,1)-path(i,1))^2 + (path(i+1,2)-path(i,2))^2);
        end
    else
        % ===== Fix: collision penalty changed from E(1)*E(2)*nB to 1e6*nB =====
        % Original: fx = E(1) * E(2) * nB;  (only 1600 on a 40x40 map, too small)
        % Now: any colliding path has fitness >= 1e6, fully separated from valid paths (55-200)
        fx = 1e6 * nB;
    end
end
