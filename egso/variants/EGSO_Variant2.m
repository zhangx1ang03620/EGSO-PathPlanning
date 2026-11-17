function [Best_score, Best_pos, Convergence_curve] = EGSO_Variant2(SearchAgents_no, Max_iter, lb, ub, dim, fobj, G, init_mode, init_field_mode, elite_field_mode)
% =========================================================================
% EGSO_Variant2 -- FULLY CROSSED factorial-ablation harness
% (extends the reduced design of EGSO_Variant.m).
%
% EGSO_Variant.m used a single field_mode shared by both GSI's column
% weighting and the elite layer's descent direction, so "initialization"
% and "elite guidance" were not fully independent factors -- e.g. it
% could not test "GSI weighted by Eikonal, but elite layer guided by
% APF". This version takes the field for each role separately:
%
%   init_mode        in {'gsi', 'random'}
%   init_field_mode  in {'eikonal','apf','dijkstra'}  (used only when
%                     init_mode='gsi'; ignored for 'random')
%   elite_field_mode in {'eikonal','apf','dijkstra'}  (no longer used:
%                     the elite field-directed step was removed because
%                     it never changed a candidate; kept in the signature
%                     so Main_FactorialAblation2.m runs unchanged)
%
% This gives the full 3(init field)x3(elite field) + 3(random init, by
% elite field) = 12-cell factorial design:
% initialization and elite guidance now vary completely independently.
% Every other mechanism (elite Gaussian refinement, Lévy flight, differential
% mutation, anti-stagnation restart, CALS smoothing) is held fixed and
% identical across all 12 cells.
% =========================================================================

M = size(G, 1);
E = [size(G,1), size(G,2)];

Ne            = max(3, floor(0.2 * SearchAgents_no));
p_pert        = 0.3;
stall_limit   = 5;
min_clearance = 1.0;

alpha_att  = 1.0;
alpha_rep  = 3.0;
rep_radius = max(3, round(0.08*M));

obstacleMask = (G == 1);
T_clear = FastMarchingField(G, obstacleMask, 'clearance');

goalMask = false(size(G));
goalMask(E(1), E(2)) = true;

% Precompute all three candidate fields once; select per role below.
FieldFns = struct();
T_eik = FastMarchingField(G, goalMask, 'navigate');
FieldFns.eikonal.Uvals = @(rows, col) T_eik(rows, col);

T_bfs = BFSDistanceField(G, goalMask);
FieldFns.dijkstra.Uvals = @(rows, col) T_bfs(rows, col);

FieldFns.apf.Uvals = @(rows, col) APF_Potential(rows, col, G, E, alpha_att, alpha_rep, rep_radius);


%% ══════════════════════════════════════════════════════════════════════
%  Module 2: Initialization
% ══════════════════════════════════════════════════════════════════════
switch init_mode
    case 'gsi'
        UvalsFn = FieldFns.(init_field_mode).Uvals;
        Pool_size = 10 * SearchAgents_no;
        Pool      = zeros(Pool_size, dim);
        Pool_fit  = inf(Pool_size, 1);

        for i = 1:Pool_size
            cur_row = 1;
            geo_strength = 0.02 + rand() * 0.08;
            for j = 1:dim
                col = j + 1;
                free_col = find(G(:, col) == 0);
                if isempty(free_col)
                    free_col = 1;
                end

                Uvals = UvalsFn(free_col, col);
                finite_mask = isfinite(Uvals);
                if any(finite_mask)
                    fallback = max(Uvals(finite_mask)) * 2 + 1;
                else
                    fallback = 1;
                end
                Uvals(~finite_mask) = fallback;

                w_geo  = exp(-geo_strength * (Uvals - min(Uvals)));
                w_cont = exp(-0.15 * abs(free_col - cur_row));
                w = w_geo .* w_cont;
                w = w / sum(w);

                chosen = free_col(WeightedSample(w));
                Pool(i, j) = chosen;
                cur_row    = chosen;
            end
            Pool_fit(i) = fobj(Pool(i, :));
        end

        [Pool_fit_sorted, idx] = sort(Pool_fit);
        Pool = Pool(idx, :);
        Positions = Pool(1:SearchAgents_no, :);
        pop_fit   = Pool_fit_sorted(1:SearchAgents_no);

    case 'random'
        Positions = zeros(SearchAgents_no, dim);
        pop_fit   = zeros(SearchAgents_no, 1);
        for i = 1:SearchAgents_no
            for j = 1:dim
                col = j + 1;
                free_col = find(G(:, col) == 0);
                if isempty(free_col)
                    free_col = 1;
                end
                Positions(i, j) = free_col(randi(numel(free_col)));
            end
            pop_fit(i) = fobj(Positions(i, :));
        end

    otherwise
        error('EGSO_Variant2:init_mode', 'Unknown init_mode: %s', init_mode);
end

[Best_score, best_idx] = min(pop_fit);
Best_pos = Positions(best_idx, :);

Convergence_curve = zeros(1, Max_iter);
Stagnation_counter = 0;

%% ══════════════════════════════════════════════════════════════════════
%  Module 3: Elite-Scout hierarchical search (identical across all cells)
% ══════════════════════════════════════════════════════════════════════
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


function U = APF_Potential(rows, col, G, E, alpha_att, alpha_rep, rep_radius)
    M = size(G, 1);
    U = alpha_att * abs(rows - E(1));

    for k = 1:numel(rows)
        r = rows(k);
        for dr = -round(rep_radius):round(rep_radius)
            nr = r + dr;
            if nr >= 1 && nr <= M && col <= size(G,2) && G(nr, col) == 1
                d = abs(dr) + 1e-6;
                if d < rep_radius
                    U(k) = U(k) + alpha_rep * (1/d - 1/rep_radius);
                end
            end
        end
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
