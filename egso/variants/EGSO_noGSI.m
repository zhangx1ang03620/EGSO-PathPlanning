function [Best_score, Best_pos, Convergence_curve] = EGSO_noGSI(SearchAgents_no, Max_iter, lb, ub, dim, fobj, G)
% =========================================================================
% EGSO ablation variant: without geodesic seeding initialization (w/o GSI)
%
% Initialisation falls back to a uniform random feasible pool (a free cell drawn
% uniformly in every column, without T_nav weighting or continuity bias).
% All other components (EVF, elite perturbation, scouts, CALS, anti-stagnation)
% are unchanged. Isolates the contribution of the GSI initialisation.
% =========================================================================

%% -- Hyperparameters ---------------------------------------------------------
M = size(G, 1);
E = [size(G,1), size(G,2)];

Ne            = max(3, floor(0.2 * SearchAgents_no));
p_pert        = 0.3;
stall_limit   = 5;
min_clearance = 1.0;

%% -- The EVF fields are still computed (EVF itself is not ablated here) -----
goalMask = false(size(G));
goalMask(E(1), E(2)) = true;
T_nav = FastMarchingField(G, goalMask, 'navigate');

obstacleMask = (G == 1);
T_clear = FastMarchingField(G, obstacleMask, 'clearance');

%% ══════════════════════════════════════════════════════════════════════════
%  Module 2 (ablated): uniform random feasible initialisation (no T_nav weighting, no continuity bias)
% ══════════════════════════════════════════════════════════════════════════
Pool_size = 10 * SearchAgents_no;
Pool      = zeros(Pool_size, dim);
Pool_fit  = inf(Pool_size, 1);

for i = 1:Pool_size
    for j = 1:dim
        col = j + 1;
        free_col = find(G(:, col) == 0);
        if isempty(free_col)
            free_col = 1;
        end
        Pool(i, j) = free_col(randi(length(free_col)));
    end
    Pool_fit(i) = fobj(Pool(i, :));
end

[Pool_fit_sorted, idx] = sort(Pool_fit);
Pool = Pool(idx, :);

Positions = Pool(1:SearchAgents_no, :);
pop_fit   = Pool_fit_sorted(1:SearchAgents_no);

Best_pos   = Positions(1, :);
Best_score = pop_fit(1);

Convergence_curve = zeros(1, Max_iter);
Stagnation_counter = 0;

%% ══════════════════════════════════════════════════════════════════════════
%  Module 3: elite-scout hierarchical search (same as the full EGSO)
% ══════════════════════════════════════════════════════════════════════════
for iter = 1:Max_iter

    [pop_fit, sort_idx] = sort(pop_fit);
    Positions = Positions(sort_idx, :);

    a = 2 * (1 - iter / Max_iter);

    for i = 1:Ne
        if rand < p_pert
            x_pert = Positions(i, :) + randn(1, dim) * a * 0.5;
            x_pert = ProjectToManifold(x_pert, lb, ub, G);
            f_pert = fobj(x_pert);
            if f_pert < pop_fit(i)
                Positions(i, :) = x_pert;
                pop_fit(i)      = f_pert;
            end
        end
    end

    if max(size(G,1), size(G,2)) > 50
        p_explore = (1 - iter/Max_iter)^3;
    else
        p_explore = 0.8 * (1 - iter/Max_iter);
    end

    for i = (Ne+1):SearchAgents_no
        if rand < p_explore
            if rand < 0.5
                levy_step = LevyFlight(dim);
                x_new = Positions(i, :) + a * 0.05 * levy_step .* (Best_pos - Positions(i, :));
            else
                candidates = setdiff(1:SearchAgents_no, i);
                r_idx = candidates(randperm(length(candidates), 2));
                r = rand;
                x_new = Positions(i, :) + r * (Best_pos - Positions(i, :)) + ...
                        (a/2) * (Positions(r_idx(1), :) - Positions(r_idx(2), :));
            end
        else
            x_new = Best_pos + randn(1, dim) * a * 0.3;
        end

        for j = 1:dim
            if x_new(j) < lb
                x_new(j) = lb + rand * (ub - lb) * 0.05;
            elseif x_new(j) > ub
                x_new(j) = ub - rand * (ub - lb) * 0.05;
            end
        end

        x_new = ProjectToManifold(x_new, lb, ub, G);
        f_new = fobj(x_new);
        if f_new < pop_fit(i)
            Positions(i, :) = x_new;
            pop_fit(i)      = f_new;
        end
    end

    [min_fit, min_idx] = min(pop_fit);
    if min_fit < Best_score
        Best_score = min_fit;
        Best_pos   = Positions(min_idx, :);
        Stagnation_counter = 0;
    else
        Stagnation_counter = Stagnation_counter + 1;
    end

    if Stagnation_counter >= stall_limit && iter > 0.2*Max_iter && iter < 0.7*Max_iter
        n_reinit = floor(0.7 * SearchAgents_no);
        for i = (SearchAgents_no - n_reinit + 1):SearchAgents_no
            for j = 1:dim
                col = min(j+1, size(G,2));
                free_col = find(G(:, col) == 0);
                if ~isempty(free_col)
                    Positions(i, j) = free_col(randi(length(free_col)));
                end
            end
            pop_fit(i) = fobj(Positions(i, :));
        end
        Stagnation_counter = 0;
    end

    if iter > 0.65 * Max_iter
        [Best_pos, Best_score] = CALS_Smooth(Best_pos, Best_score, G, fobj, T_clear, min_clearance);
        Positions(end, :) = Best_pos;
        pop_fit(end)      = Best_score;
    end

    Convergence_curve(iter) = Best_score;
end

end


function x_proj = ProjectToManifold(x, lb, ub, G)
    M      = size(G, 1);
    x_proj = max(lb, min(ub, x));
    x_proj = round(x_proj);
    dim    = length(x_proj);

    for j = 1:dim
        x_proj(j) = max(1, min(M, x_proj(j)));
        col = j + 1;
        if col <= size(G,2) && G(x_proj(j), col) == 1
            free_cells = find(G(:, col) == 0);
            if ~isempty(free_cells)
                [~, nearest] = min(abs(free_cells - x_proj(j)));
                x_proj(j)    = free_cells(nearest);
            end
        end
    end
end

function [x_out, f_out] = CALS_Smooth(x, f_current, G, fobj, T_clear, min_clearance)
    x_out = x;
    f_out = f_current;
    dim   = length(x);
    M     = size(G, 1);

    for round_iter = 1:3
        for j = 2:dim-1
            x_cand = x_out;
            x_cand(j) = round((x_out(j-1) + x_out(j+1)) / 2);
            x_cand(j) = max(min(x_cand(j), M), 1);

            col = j + 1;
            clearance_ok = (col <= size(T_clear,2)) && (T_clear(x_cand(j), col) >= min_clearance);

            if col <= size(G,2) && G(x_cand(j), col) == 0 && clearance_ok
                f_cand = fobj(x_cand);
                if f_cand < f_out
                    x_out = x_cand;
                    f_out = f_cand;
                end
            end
        end
    end
end

function L = LevyFlight(dim)
    beta  = 1.5;
    sigma = (gamma(1+beta)*sin(pi*beta/2) / ...
            (gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    u = randn(1, dim) * sigma;
    v = randn(1, dim);
    L = u ./ (abs(v).^(1/beta));
end
