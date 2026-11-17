function [Best_score, Best_pos, Curve] = MFO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G)
% =========================================================================
% MFO: Moth-Flame Optimization
% Used in the path-planning comparison experiments
% =========================================================================

% Convergence curve
Curve = zeros(1, max_gen);

% Initialise the moths
Moth_pos = rand(num_population, dimensions) .* (Xmax - X_min) + X_min;
Moth_fitness = zeros(1, num_population);

% Initialise flame positions and fitness
Flame_pos = zeros(num_population, dimensions);
Flame_fitness = inf(1, num_population);

% Record the global best
Best_pos = zeros(1, dimensions);
Best_score = inf;

% Fitness of the initial population
for i = 1:num_population
    Moth_fitness(i) = fobj(Moth_pos(i,:));
end

% Initial sort; the best moths become the initial flames
[sorted_fitness, sorted_indexes] = sort(Moth_fitness);
for i = 1:num_population
    Flame_pos(i,:) = Moth_pos(sorted_indexes(i), :);
end
Flame_fitness = sorted_fitness;

% Best solution of the initial generation
Best_score = Flame_fitness(1);
Best_pos = Flame_pos(1, :);

% Main loop
for t = 1:max_gen
    
    % Adaptive number of flames, decreasing over the run (exploration -> exploitation)
    Flame_no = round(num_population - t * ((num_population - 1) / max_gen));
    
    % a decreases linearly from -1 to -2 and controls the spiral shape
    a = -1 + t * ((-1) / max_gen);
    
    for i = 1:num_population
        for j = 1:dimensions
            
            % Coordinates of the current moth and its flame
            % Highly ranked moths spiral around their own flame; the others around the last remaining flame
            if i <= Flame_no
                distance_to_flame = abs(Flame_pos(i,j) - Moth_pos(i,j));
                b = 1; % logarithmic spiral shape constant
                t_param = (a - 1) * rand() + 1; % random parameter t in [a, 1]
                
                % Logarithmic spiral update
                Moth_pos(i,j) = distance_to_flame * exp(b * t_param) * cos(t_param * 2 * pi) + Flame_pos(i,j);
            else
                distance_to_flame = abs(Flame_pos(Flame_no,j) - Moth_pos(i,j));
                b = 1;
                t_param = (a - 1) * rand() + 1;
                
                % Lower-ranked moths follow the last valid flame
                Moth_pos(i,j) = distance_to_flame * exp(b * t_param) * cos(t_param * 2 * pi) + Flame_pos(Flame_no,j);
            end
        end
        
        % Boundary check (keep the moth inside the map)
        Moth_pos(i,:) = min(max(Moth_pos(i,:), X_min), Xmax);
        
        % Re-evaluate the moth
        Moth_fitness(i) = fobj(Moth_pos(i,:));
    end
    
    % Elitism: merge the current moths with the previous flames
    double_pos = [Flame_pos; Moth_pos];
    double_fitness = [Flame_fitness, Moth_fitness];
    
    % Sort the merged double-size population
    [double_fitness_sorted, I] = sort(double_fitness);
    
    % The best num_population individuals become the next flames
    Flame_fitness = double_fitness_sorted(1:num_population);
    Flame_pos = double_pos(I(1:num_population), :);
    
    % Update the global best
    if Flame_fitness(1) < Best_score
        Best_score = Flame_fitness(1);
        Best_pos = Flame_pos(1, :);
    end
    
    % Record the convergence curve
    Curve(t) = Best_score;
end

end