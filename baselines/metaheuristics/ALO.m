function [Elite_score, Elite_pos, Curve] = ALO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G)
% =========================================================================
% ALO: Ant Lion Optimizer
% Used in the path-planning comparison experiments
% =========================================================================

% Initialise parameters
Antlion_pos = rand(num_population, dimensions) .* (Xmax - X_min) + X_min;
Ant_pos = rand(num_population, dimensions) .* (Xmax - X_min) + X_min;

Antlion_fitness = zeros(1, num_population);
Ant_fitness = zeros(1, num_population);
Curve = zeros(1, max_gen);

% Fitness of the initial antlions
for i = 1:num_population
    Antlion_fitness(i) = fobj(Antlion_pos(i,:));
end

% Sort by fitness and take the initial elite antlion
[sorted_antlion_fitness, sorted_indexes] = sort(Antlion_fitness);
Antlion_pos = Antlion_pos(sorted_indexes, :);
Antlion_fitness = sorted_antlion_fitness;

Elite_pos = Antlion_pos(1, :);
Elite_score = Antlion_fitness(1);

% Main loop
for t = 1:max_gen
    
    % Shrinking ratio I, updated with the iteration ratio (the traps narrow over time)
    w = 1;
    if t > 0.1 * max_gen, w = 2; end
    if t > 0.5 * max_gen, w = 3; end
    if t > 0.75 * max_gen, w = 4; end
    if t > 0.9 * max_gen, w = 5; end
    if t > 0.95 * max_gen, w = 6; end
    
    I = 1 + 10^w * t / max_gen;
    
    % Shrink the trap
    c = X_min / I;
    d = Xmax / I;
    
    for i = 1:num_population
        
        % Select an antlion by roulette wheel
        % Minimisation: lower fitness should get a larger weight, hence the reciprocal
        weights = 1 ./ (Antlion_fitness + eps); 
        Rolette_index = RouletteWheelSelection(weights);
        if Rolette_index == 0
            Rolette_index = 1;
        end
        selected_antlion = Antlion_pos(Rolette_index, :);
        
        % Random walks around the selected antlion and around the elite
        RA = Random_walk_around_antlion(dimensions, max_gen, c, d, selected_antlion, t);
        RE = Random_walk_around_antlion(dimensions, max_gen, c, d, Elite_pos, t);
        
        % New ant position (average of the two random walks)
        Ant_pos(i,:) = (RA + RE) / 2;
        
        % Boundary check
        Ant_pos(i,:) = min(max(Ant_pos(i,:), X_min), Xmax);
        
        % Fitness of the new ants
        Ant_fitness(i) = fobj(Ant_pos(i,:));
        
        % Catching prey: if an ant is better than its antlion, the antlion moves to the ant
        if Ant_fitness(i) < Antlion_fitness(i)
            Antlion_pos(i,:) = Ant_pos(i,:);
            Antlion_fitness(i) = Ant_fitness(i);
        end
    end
    
    % Update the global elite antlion
    [current_best_score, best_idx] = min(Antlion_fitness);
    if current_best_score < Elite_score
        Elite_score = current_best_score;
        Elite_pos = Antlion_pos(best_idx, :);
    end
    
    % Keep the elite in the population (replaces the worst individual)
    Antlion_pos(num_population, :) = Elite_pos;
    Antlion_fitness(num_population) = Elite_score;
    
    % Record the convergence curve
    Curve(t) = Elite_score;
end

end

% =========================================================================
% Local function: roulette wheel selection
% =========================================================================
function index = RouletteWheelSelection(weights)
    accumulation = cumsum(weights);
    p = rand() * accumulation(end);
    index = find(accumulation >= p, 1);
end

% =========================================================================
% Local function: random walk around an antlion
% =========================================================================
function RW_t = Random_walk_around_antlion(Dim, iter, c, d, antlion, current_iter)
    RW_t = zeros(1, Dim);
    for j = 1:Dim
        % Standard ALO random walk (cumsum)
        steps = 2 * (rand(1, iter) > 0.5) - 1;
        X = [0, cumsum(steps)];
        
        % Min and max of the walk, used for normalisation
        a = min(X);
        b = max(X);
        
        % Avoid division by zero when b equals a
        if b == a
            b = a + eps;
        end
        
        % Move the trap bounds around the selected antlion
        c_j = c + antlion(j);
        d_j = d + antlion(j);
        
        % Normalisation: X_norm = (X - a) * (d - c) / (b - a) + c
        % Value at the current iteration t
        X_t = X(current_iter + 1); 
        RW_t(j) = ((X_t - a) * (d_j - c_j)) / (b - a) + c_j;
    end
end