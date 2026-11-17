function [Best_score, Best_pos, Curve] = FPA(num_population, max_gen, X_min, Xmax, dimensions, fobj, G)
% =========================================================================
% FPA: Flower Pollination Algorithm
% Used in the path-planning comparison experiments
% =========================================================================

% FPA-specific parameters
p = 0.8; % switch probability between global and local pollination

% Initialise the population (flowers)
X = rand(num_population, dimensions) .* (Xmax - X_min) + X_min;
Fitness = zeros(num_population, 1);
Curve = zeros(1, max_gen);

% Initial evaluation
for i = 1:num_population
    Fitness(i) = fobj(X(i,:));
end

% Initial best solution
[Best_score, idx] = min(Fitness);
Best_pos = X(idx, :);

% Main loop
for t = 1:max_gen
    for i = 1:num_population
        
        % Random number in [0, 1]
        if rand() < p
            % ===================== 1. Global pollination =====================
            % Step from a Levy flight
            L = Levy(dimensions);
            
            % Move towards the global best with a Levy-distributed step
            X_new = X(i,:) + L .* (X(i,:) - Best_pos);
            
        else
            % ===================== 2. Local pollination =====================
            % Two other distinct flowers chosen at random
            j = randi(num_population);
            k = randi(num_population);
            while j == k || j == i
                j = randi(num_population);
            end
            while k == i || k == j
                k = randi(num_population);
            end
            
            % epsilon ~ U(0, 1)
            epsilon = rand();
            
            % Move along the difference between the two random flowers
            X_new = X(i,:) + epsilon * (X(j,:) - X(k,:));
        end
        
        % Boundary check
        X_new = min(max(X_new, X_min), Xmax);
        
        % Evaluate the new position
        f_new = fobj(X_new);
        
        % Greedy selection: keep the new flower if it is better
        if f_new < Fitness(i)
            X(i,:) = X_new;
            Fitness(i) = f_new;
        end
        
        % Update the global best
        if f_new < Best_score
            Best_score = f_new;
            Best_pos = X_new;
        end
    end
    
    % Best value of this generation
    Curve(t) = Best_score;
end

end

% =========================================================================
% Local function: Levy flight step
% =========================================================================
function L = Levy(d)
    beta = 1.5;
    % Standard deviation of the normal distribution
    sigma = (gamma(1 + beta) * sin(pi * beta / 2) / ...
            (gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2)))^(1 / beta);
    
    % Two normally distributed random vectors u and v
    u = randn(1, d) * sigma;
    v = randn(1, d);
    
    % Levy step
    step = u ./ abs(v).^(1 / beta);
    
    % Scaled by 0.01 to avoid overly large jumps that hinder convergence
    L = 0.01 * step;
end