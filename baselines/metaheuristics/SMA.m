function [Best_score, Best_pos, Curve] = SMA(num_population, max_gen, X_min, Xmax, dimensions, fobj, G)
% =========================================================================
% SMA: Slime Mould Algorithm
% Used in the path-planning comparison experiments
% =========================================================================

% Initialise the population
X = rand(num_population, dimensions) .* (Xmax - X_min) + X_min;
Fitness = zeros(num_population, 1);
weight = ones(num_population, dimensions); % vein weights of the slime mould
Curve = zeros(1, max_gen);

Best_pos = zeros(1, dimensions);
Best_score = inf;

% Initial evaluation
for i = 1:num_population
    Fitness(i) = fobj(X(i,:));
end

% Main loop
for t = 1:max_gen
    
    % Sort the population by fitness (minimisation)
    [smellOrder, SortIndex] = sort(Fitness);
    X = X(SortIndex, :);
    
    % Update the global best
    if smellOrder(1) < Best_score
        Best_score = smellOrder(1);
        Best_pos = X(1, :);
    end
    
    worstFitness = smellOrder(end);
    bestFitness = smellOrder(1);

    % Weights W (thickness of the veins)
    % Fitter individuals get thicker veins (larger weight), worse ones thinner veins
    for i = 1:num_population
        if i <= num_population / 2  % better half
            weight(i,:) = 1 + rand(1, dimensions) .* log10((bestFitness - smellOrder(i)) / (bestFitness - worstFitness + eps) + 1);
        else                        % worse half
            weight(i,:) = 1 - rand(1, dimensions) .* log10((bestFitness - smellOrder(i)) / (bestFitness - worstFitness + eps) + 1);
        end
    end
    
    % Control parameters a, vb, vc
    a = atanh(-(t / max_gen) + 1); % a decreases from inf to 0
    b = 1 - t / max_gen;           % b decreases linearly
    
    for i = 1:num_population
        % Random vb and vc
        vb = (rand(1, dimensions) * 2 * a) - a; % random numbers in [-a, a]
        vc = (rand(1, dimensions) * 2 * b) - b; % random numbers in [-b, b]
        
        % Control threshold p
        p = tanh(abs(Fitness(i) - bestFitness));
        
        % Position update
        for j = 1:dimensions
            r1 = rand(); % probability of a random re-initialisation
            r2 = rand(); % selects the update rule
            
            if r1 < 0.03 % z = 0.03: small chance of a fully random restart
                X(i,j) = rand() * (Xmax - X_min) + X_min;
            else
                if r2 < p
                    % Update from the best position, the weights and two random individuals
                    rand_idx1 = randi([1, num_population]);
                    rand_idx2 = randi([1, num_population]);
                    X(i,j) = Best_pos(j) + vb(j) * (weight(i,j) * X(rand_idx1, j) - X(rand_idx2, j));
                else
                    % Local search: contract towards the current position
                    X(i,j) = vc(j) * X(i,j);
                end
            end
        end
        
        % Boundary check
        X(i,:) = min(max(X(i,:), X_min), Xmax);
        
        % Re-evaluate the new position
        Fitness(i) = fobj(X(i,:));
    end
    
    % Best value of this generation
    Curve(t) = Best_score;
end

end