% SBOA -- Secretary Bird Optimization Algorithm
% Fu, Y., Liu, D., Chen, J., & He, L. (2024).
% Adapted to the grid path-planning encoding used throughout this project
% (same interface/conventions as the other comparison algorithms).
%
% Three behavioural phases split across the iteration budget:
%   (1) searching for prey     -- random-walk exploration
%   (2) consuming prey         -- directed exploitation toward the best
%   (3) escape strategy        -- fight (local perturbation) or flight
%                                  (large random jump), chosen at random
function [Best_score, Best_pos, Convergence_curve] = SBOA(SearchAgents_no, Max_iter, lb, ub, dim, fobj, G)

Positions = zeros(SearchAgents_no, dim);
for i = 1:SearchAgents_no
    for j = 1:dim
        column = G(:, j+1);
        id = find(column == 0);
        Positions(i,j) = id(randi(length(id)));
    end
end

fitness_val = zeros(SearchAgents_no, 1);
for i = 1:SearchAgents_no
    fitness_val(i) = fobj(Positions(i,:));
end
[Best_score, idx] = min(fitness_val);
Best_pos = Positions(idx,:);

Convergence_curve = zeros(1, Max_iter);

for t = 1:Max_iter
    frac = t / Max_iter;

    for i = 1:SearchAgents_no
        r_idx = randi(SearchAgents_no);

        if frac < 1/3
            % --- phase 1: searching for prey (exploration) ---
            x_new = Positions(i,:) + (Positions(r_idx,:) - Positions(i,:)) .* randn(1, dim) * 0.5;
        elseif frac < 2/3
            % --- phase 2: consuming prey (exploitation) ---
            CF = (1 - frac)^2;
            x_new = Best_pos + CF * (rand(1, dim) - 0.5) .* (Best_pos - Positions(i,:));
        else
            % --- phase 3: escape strategy (fight or flight) ---
            if rand < 0.5
                x_new = Positions(i,:) + 0.2*(1-frac) * randn(1, dim) .* (Best_pos - Positions(i,:));
            else
                x_new = Positions(i,:) + (2*rand(1, dim)-1) .* (ub - lb) * 0.1 * (1-frac);
            end
        end

        x_new = round(x_new);
        x_new = boundaryCheck(x_new, lb, ub);

        f_new = fobj(x_new);
        if f_new < fitness_val(i)
            Positions(i,:) = x_new;
            fitness_val(i) = f_new;
            if f_new < Best_score
                Best_score = f_new;
                Best_pos   = x_new;
            end
        end
    end

    Convergence_curve(t) = Best_score;
    Best_pos   = LocalSearch(Best_pos, ub, G);
    Best_score = fobj(Best_pos);
end

end
