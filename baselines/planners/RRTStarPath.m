function [path, path_length, num_nodes] = RRTStarPath(G, S, E, max_iter, informed)
% RRTStarPath  Sampling-based path planner: standard RRT* (Karaman & Frazzoli,
% 2011), with an optional Informed-RRT* (Gammell et al., 2014) sampling
% restriction once an initial solution is found.
%
% Represents the sampling-based paradigm (the RRT* family) alongside the
% graph-search (A*) and field-descent (FMM-descent) references.
%
%   G   : MxM binary occupancy grid (1 = obstacle, 0 = free)
%   S,E : [row, col] start / goal cells
%   max_iter : sample budget
%   informed : true -> Informed-RRT* (ellipsoidal sampling after first
%              solution is found); false -> plain RRT*
%
% Returns path (Nx2 [row col], empty if none found), its Euclidean length
% (Inf if none found), and the number of tree nodes actually grown.
%
% Continuous state space [1,M]x[1,M]; collision checking is done by
% discretized line-of-sight sampling against the binary grid (nearest-cell
% occupancy), consistent with the continuous-waypoint interpretation used
% throughout this project's smoothing/length post-processing.

    M = size(G, 1);
    eta          = 3.0;     % steer step size (grid units)
    goal_tol     = 2.0;     % goal-region radius
    goal_bias    = 0.08;    % probability of sampling the goal directly
    gamma_rrg    = 25.0;    % RRT* rewiring constant (tuned for this domain scale)

    max_nodes = max_iter + 1;
    pos    = zeros(max_nodes, 2);
    parent = zeros(max_nodes, 1);
    cost   = zeros(max_nodes, 1);

    pos(1,:) = S;
    parent(1) = 0;
    cost(1)   = 0;
    n = 1;

    best_goal_idx  = -1;
    best_goal_cost = inf;
    c_min = norm(E - S);
    x_center = (S + E) / 2;
    theta = atan2(E(1)-S(1), E(2)-S(2));
    Rrot = [cos(theta) -sin(theta); sin(theta) cos(theta)];

    for iter = 1:max_iter
        % ---- sampling ----
        if rand < goal_bias
            x_rand = E;
        elseif informed && best_goal_idx > 0
            x_rand = sampleEllipse(x_center, Rrot, c_min, best_goal_cost, M);
        else
            x_rand = [1 + rand*(M-1), 1 + rand*(M-1)];
        end

        % ---- nearest node ----
        d2 = sum((pos(1:n,:) - x_rand).^2, 2);
        [~, idx_near] = min(d2);
        x_near = pos(idx_near, :);

        % ---- steer ----
        dvec = x_rand - x_near;
        dist = norm(dvec);
        if dist < 1e-9
            continue;
        end
        if dist > eta
            x_new = x_near + (eta/dist) * dvec;
        else
            x_new = x_rand;
        end
        x_new = max(min(x_new, [M M]), [1 1]);

        if ~collisionFreeSegment(x_near, x_new, G)
            continue;
        end

        % ---- choose best parent among near nodes ----
        d2n = sum((pos(1:n,:) - x_new).^2, 2);
        r_n = min(gamma_rrg * sqrt(log(n+1)/(n+1)), 3*eta);
        near_idx = find(d2n <= r_n^2);
        if isempty(near_idx)
            near_idx = idx_near;
        end

        best_parent = idx_near;
        best_cost   = cost(idx_near) + norm(x_new - x_near);
        for k = 1:numel(near_idx)
            j = near_idx(k);
            c_candidate = cost(j) + norm(x_new - pos(j,:));
            if c_candidate < best_cost && collisionFreeSegment(pos(j,:), x_new, G)
                best_cost   = c_candidate;
                best_parent = j;
            end
        end

        n = n + 1;
        pos(n,:)  = x_new;
        parent(n) = best_parent;
        cost(n)   = best_cost;

        % ---- rewire ----
        for k = 1:numel(near_idx)
            j = near_idx(k);
            if j == best_parent, continue; end
            c_candidate = cost(n) + norm(pos(j,:) - x_new);
            if c_candidate < cost(j) - 1e-9 && collisionFreeSegment(x_new, pos(j,:), G)
                parent(j) = n;
                cost(j)   = c_candidate;
            end
        end

        % ---- goal check ----
        if norm(x_new - E) <= goal_tol
            if cost(n) < best_goal_cost
                best_goal_cost = cost(n);
                best_goal_idx  = n;
            end
        end
    end

    num_nodes = n;
    if best_goal_idx < 0
        path = [];
        path_length = Inf;
        return;
    end

    % ---- reconstruct path ----
    chain = best_goal_idx;
    j = best_goal_idx;
    while parent(j) ~= 0
        j = parent(j);
        chain(end+1) = j; %#ok<AGROW>
    end
    chain = fliplr(chain);
    path = pos(chain, :);
    path = [path; E];   % snap final segment exactly onto the goal cell

    path_length = 0;
    for i = 1:size(path,1)-1
        path_length = path_length + norm(path(i+1,:) - path(i,:));
    end
end

%% ========================================================================
function ok = collisionFreeSegment(p1, p2, G)
    M = size(G, 1);
    dist = norm(p2 - p1);
    nsteps = max(2, ceil(dist / 0.5));
    for t = linspace(0, 1, nsteps)
        p = p1 + t * (p2 - p1);
        r = max(1, min(M, round(p(1))));
        c = max(1, min(M, round(p(2))));
        if G(r, c) == 1
            ok = false;
            return;
        end
    end
    ok = true;
end

%% ========================================================================
function x = sampleEllipse(x_center, Rrot, c_min, c_best, M)
    if c_best <= c_min
        r1 = c_best / 2;
        r2 = r1;
    else
        r1 = c_best / 2;
        r2 = sqrt(max(c_best^2 - c_min^2, 0)) / 2;
    end
    for attempt = 1:20
        u = rand(1,2)*2 - 1;
        if sum(u.^2) <= 1
            break;
        end
    end
    scaled = [r1*u(1), r2*u(2)];
    offset = (Rrot * scaled')';
    x = x_center + offset;
    x = max(min(x, [M M]), [1 1]);
end
