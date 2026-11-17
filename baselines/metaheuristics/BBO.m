function [Best_score, Best_pos, Curve] = BBO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G)
% =========================================================================
% BBO: Biogeography-Based Optimization
% Used in the path-planning comparison experiments
% =========================================================================

% BBO parameters
KeepRate = 0.2; % elitism ratio (the best solutions are not modified)
KeepNum = max(1, round(KeepRate * num_population));
MutationProbability = 0.05; % mutation probability

% Initialise the habitats (population)
Habitats = rand(num_population, dimensions) .* (Xmax - X_min) + X_min;
Costs = zeros(num_population, 1);
Curve = zeros(1, max_gen);

% Initial evaluation
for i = 1:num_population
    Costs(i) = fobj(Habitats(i,:));
end

% Sort by fitness, ascending (minimisation)
[Costs, SortOrder] = sort(Costs);
Habitats = Habitats(SortOrder, :);

% Record the global best
Best_pos = Habitats(1, :);
Best_score = Costs(1);

% Main loop
for t = 1:max_gen
    
    % Immigration rate (lambda) and emigration rate (mu) of every habitat
    % Better rank -> larger emigration rate mu and smaller immigration rate lambda
    rank = 1:num_population;
    mu = (num_population - rank) ./ num_population;
    lambda = rank ./ num_population;
    
    NewHabitats = Habitats;
    
    for i = 1:num_population
        % Elitism: the first KeepNum habitats are not modified
        if i > KeepNum
            for j = 1:dimensions
                
                % Migration
                % The immigration rate lambda decides whether a feature is replaced
                if rand() < lambda(i)
                    % Source habitat chosen by roulette wheel on the emigration rate mu
                    Source_idx = RouletteWheelSelection(mu);
                    NewHabitats(i, j) = Habitats(Source_idx, j);
                end
                
                % Mutation
                if rand() < MutationProbability
                    NewHabitats(i, j) = X_min + rand() * (Xmax - X_min);
                end
            end
            
            % Boundary check
            NewHabitats(i, :) = min(max(NewHabitats(i, :), X_min), Xmax);
            
            % Re-evaluate the modified habitats
            Costs(i) = fobj(NewHabitats(i, :));
        end
    end
    
    % Re-sort the modified population
    Habitats = NewHabitats;
    [Costs, SortOrder] = sort(Costs);
    Habitats = Habitats(SortOrder, :);
    
    % Update the global best
    if Costs(1) < Best_score
        Best_score = Costs(1);
        Best_pos = Habitats(1, :);
    end
    
    % Best value of this generation
    Curve(t) = Best_score;
end

end

% =========================================================================
% Local function: roulette wheel selection
% =========================================================================
function index = RouletteWheelSelection(probs)
    accumulation = cumsum(probs);
    p = rand() * accumulation(end);
    index = find(accumulation >= p, 1);
    if isempty(index)
        index = 1;
    end
end