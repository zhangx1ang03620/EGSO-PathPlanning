function [Best_score, Best_pos, Curve] = ABC(num_population, max_gen, X_min, Xmax, dimensions, fobj, G)
% =========================================================================
% ABC: Artificial Bee Colony
% Used in the path-planning comparison experiments
% =========================================================================

% ABC-specific parameters
num_food_sources = round(num_population / 2); % number of food sources (usually half the population size)
limit = round(0.6 * dimensions * num_food_sources); % abandonment threshold of a food source

% Initialise variables
Positions = rand(num_food_sources, dimensions) .* (Xmax - X_min) + X_min;
Fitness = zeros(1, num_food_sources);
Trials = zeros(1, num_food_sources); % number of consecutive non-improvements of each food source
Curve = zeros(1, max_gen);

Best_pos = zeros(1, dimensions);
Best_score = inf;

% Evaluate the initial population
for i = 1:num_food_sources
    Fitness(i) = fobj(Positions(i,:));
    if Fitness(i) < Best_score
        Best_score = Fitness(i);
        Best_pos = Positions(i,:);
    end
end

% Main loop
for t = 1:max_gen
    
    % ===================== 1. Employed bees =====================
    for i = 1:num_food_sources
        % Pick a random neighbouring food source k ~= i
        k = randi([1, num_food_sources]);
        while k == i
            k = randi([1, num_food_sources]);
        end
        
        % Pick a random dimension j to update
        j = randi([1, dimensions]);
        
        % New candidate from the difference between the current source and the random neighbour
        phi = (rand() - 0.5) * 2; % random number in [-1, 1]
        new_pos = Positions(i,:);
        new_pos(j) = Positions(i,j) + phi * (Positions(i,j) - Positions(k,j));
        
        % Boundary check
        new_pos = min(max(new_pos, X_min), Xmax);
        
        % Evaluate the new candidate
        new_fitness = fobj(new_pos);
        
        % Greedy selection
        if new_fitness < Fitness(i)
            Positions(i,:) = new_pos;
            Fitness(i) = new_fitness;
            Trials(i) = 0; % improved: reset the counter
        else
            Trials(i) = Trials(i) + 1; % not improved: increment the counter
        end
    end
    
    % ===================== 2. Onlooker bees =====================
    % Selection probabilities (minimisation: lower fitness -> higher probability)
    fit_probs = 1 ./ (Fitness + eps);
    prob = fit_probs / sum(fit_probs);
    
    i = 1; 
    m = 0; % number of onlooker bees dispatched so far
    while m < num_food_sources
        % Roulette wheel: an onlooker visits the source if the random number is below its probability
        if rand() < prob(i)
            m = m + 1;
            
            % Same local search rule as the employed bees
            k = randi([1, num_food_sources]);
            while k == i
                k = randi([1, num_food_sources]);
            end
            
            j = randi([1, dimensions]);
            phi = (rand() - 0.5) * 2;
            new_pos = Positions(i,:);
            new_pos(j) = Positions(i,j) + phi * (Positions(i,j) - Positions(k,j));
            
            % Boundary check
            new_pos = min(max(new_pos, X_min), Xmax);
            
            % Evaluate the new candidate and select greedily
            new_fitness = fobj(new_pos);
            if new_fitness < Fitness(i)
                Positions(i,:) = new_pos;
                Fitness(i) = new_fitness;
                Trials(i) = 0;
            else
                Trials(i) = Trials(i) + 1;
            end
        end
        
        i = i + 1;
        if i > num_food_sources
            i = 1; % wrap around
        end
    end
    
    % ===================== 3. Scout bees =====================
    for i = 1:num_food_sources
        % A source not improved for more than `limit` trials is abandoned and re-initialised at random
        if Trials(i) >= limit
            Positions(i,:) = rand(1, dimensions) .* (Xmax - X_min) + X_min;
            Fitness(i) = fobj(Positions(i,:));
            Trials(i) = 0;
        end
    end
    
    % ===================== Record and update the best solution =====================
    for i = 1:num_food_sources
        if Fitness(i) < Best_score
            Best_score = Fitness(i);
            Best_pos = Positions(i,:);
        end
    end
    
    % Best value of this generation
    Curve(t) = Best_score;
end

end