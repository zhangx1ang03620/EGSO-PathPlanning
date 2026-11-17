% Artificial Rabbit Optimization (ARO)
% Wang, L., Cao, Q., Zhang, Z., Mirjalili, S., & Zhao, W. (2022).
% Adapted to the grid path-planning encoding used throughout this project
% (same interface/conventions as the other comparison algorithms).
%
% Two foraging strategies balanced by an energy factor A(t) that decays
% over the run: (1) detour foraging (exploration) -- a rabbit hops toward
% a randomly chosen peer's territory; (2) random hiding (exploitation) --
% a rabbit shrinks back toward the current best with a decreasing random
% step (mimicking retreat into a nearby burrow).
function [Best_score, Best_pos, Convergence_curve] = ARO(SearchAgents_no, Max_iter, lb, ub, dim, fobj, G)

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
    A = 4 * (1 - t/Max_iter) * log(1/rand);   % energy factor, decays over the run

    for i = 1:SearchAgents_no
        r_idx = randi(SearchAgents_no);
        while r_idx == i && SearchAgents_no > 1
            r_idx = randi(SearchAgents_no);
        end

        if A > 1
            % --- detour foraging (exploration) ---
            R = zeros(1, dim);
            R(randi(dim)) = 1;
            x_new = Positions(i,:) + R .* (Positions(r_idx,:) - Positions(i,:)) ...
                    + round(0.5*(0.05+rand)) * randn(1, dim);
        else
            % --- random hiding (exploitation) ---
            g = zeros(1, dim);
            num_pick = max(1, round(dim * rand));
            pick = randperm(dim, num_pick);
            g(pick) = 1;
            shrink = (Max_iter - t + 1) / Max_iter;
            x_new = Positions(i,:) + rand * (Best_pos - Positions(i,:)) .* g ...
                    + shrink * randn(1, dim) .* (1 - g);
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
