function [Best_score, Best_pos, Convergence_curve] = EGSO3D(SearchAgents_no, Max_iter, dim, fobj, G)
% =========================================================================
% EGSO3D -- minimal 3D proof-of-concept port of EGSO.m (see EGSO.m for the
% full docstring of the four core modules). Decision vector length is
% 2*dim: row indices in x(1:dim), height indices in x(dim+1:2*dim), one
% pair per intermediate column (dim = C-2). Everything else (EVF/GSI/
% elite-scout search/CALS) mirrors the 2D algorithm, generalized from a
% 1D per-column sample to a 2D (row,height) per-column sample.
%
% This is a deliberately minimal port for a proof-of-concept demonstration
% that the EGSO design paradigm transfers to voxelized 3D domains -- NOT
% a fully re-tuned, paper-grade 3D algorithm.
% =========================================================================

[R, C, H] = size(G);
E = [R, C, H];

Ne            = max(3, floor(0.2 * SearchAgents_no));
p_pert        = 0.3;
stall_limit   = 5;
min_clearance = 1.0;

%% Module 1: EVF
goalMask = false(R,C,H); goalMask(E(1),E(2),E(3)) = true;
T_nav = FastMarchingField3D(G, goalMask, 'navigate');

obstacleMask = (G == 1);
T_clear = FastMarchingField3D(G, obstacleMask, 'clearance');

%% Module 2: GSI
Pool_size = 10 * SearchAgents_no;
Pool      = zeros(Pool_size, 2*dim);
Pool_fit  = inf(Pool_size, 1);

for i = 1:Pool_size
    cur_r = 1; cur_h = 1;
    geo_strength = 0.02 + rand() * 0.08;
    for j = 1:dim
        col = j + 1;
        slice = squeeze(G(:, col, :) == 0);
        [fr, fh] = find(slice);
        if isempty(fr)
            fr = 1; fh = 1;
        end
        Tvals = zeros(numel(fr),1);
        for kk = 1:numel(fr)
            Tvals(kk) = T_nav(fr(kk), col, fh(kk));
        end
        finite_mask = isfinite(Tvals);
        if any(finite_mask)
            fallback = max(Tvals(finite_mask)) * 2 + 1;
        else
            fallback = 1;
        end
        Tvals(~finite_mask) = fallback;

        w_geo  = exp(-geo_strength * (Tvals - min(Tvals)));
        w_cont = exp(-0.15 * (abs(fr - cur_r) + abs(fh - cur_h)));
        w = w_geo .* w_cont;
        w = w / sum(w);

        chosen = WeightedSample(w);
        Pool(i, j)     = fr(chosen);
        Pool(i, dim+j) = fh(chosen);
        cur_r = fr(chosen); cur_h = fh(chosen);
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

%% Module 3: elite-scout hierarchical search (elite Gaussian refinement)
for iter = 1:Max_iter
    [pop_fit, sort_idx] = sort(pop_fit);
    Positions = Positions(sort_idx, :);

    a = 2 * (1 - iter / Max_iter);

    for i = 1:Ne
        if rand < p_pert
            x_pert = Positions(i, :) + randn(1, 2*dim) * a * 0.5;
            x_pert = ProjectToManifold3D(x_pert, G, dim);
            f_pert = fobj(x_pert);
            if f_pert < pop_fit(i)
                Positions(i, :) = x_pert;
                pop_fit(i)      = f_pert;
            end
        end
    end

    p_explore = 0.8 * (1 - iter/Max_iter);

    for i = (Ne+1):SearchAgents_no
        if rand < p_explore
            if rand < 0.5
                levy_step = LevyFlight3D(2*dim);
                x_new = Positions(i, :) + a * 0.05 * levy_step .* (Best_pos - Positions(i, :));
            else
                candidates = setdiff(1:SearchAgents_no, i);
                r_idx = candidates(randperm(length(candidates), 2));
                r = rand;
                x_new = Positions(i, :) + r * (Best_pos - Positions(i, :)) + ...
                        (a/2) * (Positions(r_idx(1), :) - Positions(r_idx(2), :));
            end
        else
            x_new = Best_pos + randn(1, 2*dim) * a * 0.3;
        end

        x_new = ProjectToManifold3D(x_new, G, dim);
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
                col = min(j+1, C);
                slice = squeeze(G(:, col, :) == 0);
                [fr, fh] = find(slice);
                if ~isempty(fr)
                    pick = randi(numel(fr));
                    Positions(i, j)     = fr(pick);
                    Positions(i, dim+j) = fh(pick);
                end
            end
            pop_fit(i) = fobj(Positions(i, :));
        end
        Stagnation_counter = 0;
    end

    if iter > 0.65 * Max_iter
        [Best_pos, Best_score] = CALS_Smooth3D(Best_pos, Best_score, G, fobj, T_clear, min_clearance, dim);
        Positions(end, :) = Best_pos;
        pop_fit(end)      = Best_score;
    end

    Convergence_curve(iter) = Best_score;
end

end


%% ════════════════════════════════════════════════════════════════════════════
function x_proj = ProjectToManifold3D(x, G, dim)
    [R, C, H] = size(G);
    x_proj = round(x);
    for j = 1:dim
        r = max(1, min(R, x_proj(j)));
        h = max(1, min(H, x_proj(dim+j)));
        col = min(j+1, C);
        if G(r, col, h) == 1
            slice = squeeze(G(:, col, :) == 0);
            [fr, fh] = find(slice);
            if ~isempty(fr)
                d2 = (fr - r).^2 + (fh - h).^2;
                [~, nearest] = min(d2);
                r = fr(nearest); h = fh(nearest);
            end
        end
        x_proj(j) = r;
        x_proj(dim+j) = h;
    end
end

%% ════════════════════════════════════════════════════════════════════════════
function [x_out, f_out] = CALS_Smooth3D(x, f_current, G, fobj, T_clear, min_clearance, dim)
    x_out = x; f_out = f_current;
    [R, C, H] = size(G);

    for round_iter = 1:3
        for j = 2:dim-1
            x_cand = x_out;
            r_cand = round((x_out(j-1) + x_out(j+1)) / 2);
            h_cand = round((x_out(dim+j-1) + x_out(dim+j+1)) / 2);
            r_cand = max(1, min(R, r_cand));
            h_cand = max(1, min(H, h_cand));
            x_cand(j) = r_cand; x_cand(dim+j) = h_cand;

            col = min(j+1, C);
            clearance_ok = T_clear(r_cand, col, h_cand) >= min_clearance;

            if G(r_cand, col, h_cand) == 0 && clearance_ok
                f_cand = fobj(x_cand);
                if f_cand < f_out
                    x_out = x_cand; f_out = f_cand;
                end
            end
        end
    end
end

%% ════════════════════════════════════════════════════════════════════════════
function L = LevyFlight3D(n)
    beta  = 1.5;
    sigma = (gamma(1+beta)*sin(pi*beta/2) / ...
            (gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    u = randn(1, n) * sigma;
    v = randn(1, n);
    L = u ./ (abs(v).^(1/beta));
end
