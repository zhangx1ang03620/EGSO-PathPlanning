function U = ArtificialPotentialField(G, E, katt, krep, rho0)
% ArtificialPotentialField  Standard hand-crafted APF (Khatib 1986 form),
%                            built ONLY from local information (Euclidean
%                            distance to goal + Euclidean distance to the
%                            nearest obstacle cell via a distance
%                            transform), exactly the class of guidance
%                            field criticized in the paper's Related Work
%                            (Sec 2.1) as "composed entirely of local
%                            gradients ... structurally prone to spurious
%                            local minima at non-target points."
%
%   U = ArtificialPotentialField(G, E, katt, krep, rho0)
%
%   G     : MxM binary occupancy grid (1 = obstacle, 0 = free)
%   E     : [row, col] goal cell
%   katt  : attractive gain (paraboloidal well centered at E)
%   krep  : repulsive gain
%   rho0  : repulsive field radius of influence (grid units)
%
%   Returns U (MxM double), Inf at obstacle cells (excluded from the
%   traversable domain, consistent with FastMarchingField's 'navigate'
%   mode).
%
%   U(q) = 0.5*katt*d_goal(q)^2 + 0.5*krep*(1/d_obs(q) - 1/rho0)^2  [d_obs<=rho0]
%   U(q) = 0.5*katt*d_goal(q)^2                                     [d_obs> rho0]

    [M, N] = size(G);
    [cc, rr] = meshgrid(1:N, 1:M);

    d_goal = sqrt((rr - E(1)).^2 + (cc - E(2)).^2);
    U_att  = 0.5 * katt * d_goal.^2;

    d_obs = EuclideanDistanceToObstacle(G, rr, cc);

    U_rep = zeros(M, N);
    influenced = d_obs <= rho0 & d_obs > 0;
    U_rep(influenced) = 0.5 * krep * (1 ./ d_obs(influenced) - 1/rho0).^2;
    U_rep(d_obs == 0) = inf;   % obstacle cells themselves (defensive; excluded below anyway)

    U = U_att + U_rep;
    U(G == 1) = inf;
end

%% ========================================================================
%  Brute-force straight-line Euclidean distance transform to the nearest
%  obstacle cell (no Image Processing Toolbox dependency; this is exactly
%  the local, line-of-sight distance a hand-crafted APF conventionally
%  uses -- it does NOT see around corners, unlike the geodesic T_clear
%  field used internally by EGSO's CALS module).
%% ========================================================================
function d_obs = EuclideanDistanceToObstacle(G, rr, cc)
    [M, N] = size(G);
    if ~any(G(:))
        d_obs = inf(M, N);
        return;
    end
    [obs_r, obs_c] = find(G == 1);
    d2_min = inf(M, N);
    chunk = 200;   % process obstacle cells in chunks to bound memory
    for start = 1:chunk:numel(obs_r)
        idx = start:min(start+chunk-1, numel(obs_r));
        for k = idx
            d2 = (rr - obs_r(k)).^2 + (cc - obs_c(k)).^2;
            d2_min = min(d2_min, d2);
        end
    end
    d_obs = sqrt(d2_min);
end
