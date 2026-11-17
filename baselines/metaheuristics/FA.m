function [Best_score, Best_pos, Curve] = FA(num_population, max_gen, X_min, Xmax, dimensions, fobj, G)
% =========================================================================
% FA: Firefly Algorithm
% Used in the path-planning comparison experiments
% =========================================================================

% FA-specific parameters
alpha_init = 0.5;      % initial randomisation step
beta0 = 1.0;           % attractiveness at r = 0
gamma = 1.0;           % light absorption coefficient

% Initialise the population
Positions = rand(num_population, dimensions) .* (Xmax - X_min) + X_min;
Fitness = zeros(num_population, 1);
Curve = zeros(1, max_gen);

Best_pos = zeros(1, dimensions);
Best_score = inf;

% Evaluate the initial population
for i = 1:num_population
    Fitness(i) = fobj(Positions(i,:));
    if Fitness(i) < Best_score
        Best_score = Fitness(i);
        Best_pos = Positions(i,:);
    end
end

% Main loop
for t = 1:max_gen
    
    % alpha decreases linearly with the iteration count, favouring fine search late in the run
    alpha = alpha_init * (1 - t / max_gen);
    
    for i = 1:num_population
        for j = 1:num_population
            
            % If firefly j is brighter (lower fitness) than firefly i, i moves towards j
            if Fitness(j) < Fitness(i)
                
                % Euclidean distance r between the two fireflies
                r = norm(Positions(i,:) - Positions(j,:));
                
                % Attractiveness beta after distance and light absorption
                beta = beta0 * exp(-gamma * r^2);
                
                % Update the position of firefly i
                for k = 1:dimensions
                    % current position + attraction towards j + random perturbation
                    rand_val = rand() - 0.5; % random number in [-0.5, 0.5]
                    Positions(i,k) = Positions(i,k) + beta * (Positions(j,k) - Positions(i,k)) ...
                                     + alpha * (Xmax - X_min) * rand_val;
                end
                
                % Boundary check
                Positions(i,:) = min(max(Positions(i,:), X_min), Xmax);
                
                % Re-evaluate after the move
                Fitness(i) = fobj(Positions(i,:));
                
                % Update the global best
                if Fitness(i) < Best_score
                    Best_score = Fitness(i);
                    Best_pos = Positions(i,:);
                end
            end
        end
    end
    
    % Best value of this generation
    Curve(t) = Best_score;
end

end