% INFO -- weIghted meaN oF vectOrs
% Ahmadianfar, I., Heidari, A. A., Gandomi, A. H., Chu, X., & Chen, H. (2022).
% Adapted to the grid path-planning encoding used throughout this project
% (same interface/conventions as the other comparison algorithms).
%
% Candidate vectors are generated from a fitness-weighted combination of
% the current best individual and two randomly chosen peers (the
% "weighted mean of vectors" rule), refined by a vector-combining local
% search term that pulls toward the best solution with a decreasing step.
function [Best_score, Best_pos, Convergence_curve] = INFO(SearchAgents_no, Max_iter, lb, ub, dim, fobj, G)

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
[fitness_val, sidx] = sort(fitness_val);
Positions   = Positions(sidx,:);
Best_score  = fitness_val(1);
Best_pos    = Positions(1,:);
Worst_score = fitness_val(end);

Convergence_curve = zeros(1, Max_iter);

for t = 1:Max_iter
    denom = (Worst_score - Best_score) + eps;
    w = 1 - (fitness_val - Best_score) ./ denom;   % best -> weight ~1, worst -> ~0

    for i = 1:SearchAgents_no
        idxs = randperm(SearchAgents_no, min(3, SearchAgents_no));
        a = idxs(1);
        b = idxs(min(2, numel(idxs)));
        c = idxs(min(3, numel(idxs)));

        % weighted-mean rule combining best + two random peers
        WM = w(1)*Best_pos + w(a)*Positions(a,:) + w(b)*Positions(b,:);
        WM = WM / (w(1) + w(a) + w(b) + eps);

        delta = 2 * rand(1, dim) .* abs(Positions(a,:) - Positions(b,:));
        x_new = WM + (rand(1, dim) - 0.5) .* delta * (1 - t/Max_iter);

        % vector-combining local search stage
        if rand < 0.5
            x_new = x_new + randn(1, dim) .* (Best_pos - Positions(c,:)) * 0.5 * (1 - t/Max_iter);
        end

        x_new = round(x_new);
        x_new = boundaryCheck(x_new, lb, ub);

        f_new = fobj(x_new);
        if f_new < fitness_val(i)
            Positions(i,:) = x_new;
            fitness_val(i) = f_new;
        end
    end

    [fitness_val, sidx] = sort(fitness_val);
    Positions = Positions(sidx,:);
    if fitness_val(1) < Best_score
        Best_score = fitness_val(1);
        Best_pos   = Positions(1,:);
    end
    Worst_score = fitness_val(end);

    Convergence_curve(t) = Best_score;
    Best_pos   = LocalSearch(Best_pos, ub, G);
    Best_score = fobj(Best_pos);
end

end
