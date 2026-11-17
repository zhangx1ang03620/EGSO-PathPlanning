function [Best_score, Best_pos, Curve] = SCA(num_population, max_gen, X_min, Xmax, dimensions, fobj, G)
% =========================================================================
% SCA: Sine Cosine Algorithm
% Used in the path-planning comparison experiments
% =========================================================================

% Initialise parameters
Best_pos = zeros(1, dimensions);
Best_score = inf;
Curve = zeros(1, max_gen);

% Initialise the population
X = rand(num_population, dimensions) .* (Xmax - X_min) + X_min;

% Constant a, controls the decrease of r1
a = 2;

% Main loop
for t = 1:max_gen
    
    % Boundary check and fitness evaluation
    for i = 1:num_population
        % Boundary handling
        Flag4ub = X(i,:) > Xmax;
        Flag4lb = X(i,:) < X_min;
        X(i,:) = (X(i,:) .* (~(Flag4ub + Flag4lb))) + Xmax .* Flag4ub + X_min .* Flag4lb;
        
        % Fitness
        fitness = fobj(X(i,:));
        
        % Update the global best
        if fitness < Best_score
            Best_score = fitness;
            Best_pos = X(i,:);
        end
    end
    
    % Update r1 (decreases linearly from a to 0)
    % r1 sets the range of the next move
    r1 = a - t * (a / max_gen);
    
    % Update every individual
    for i = 1:num_population
        for j = 1:dimensions
            
            % Random numbers controlling exploration and exploitation
            r2 = (2 * pi) * rand(); % r2 in [0, 2*pi], direction of the move
            r3 = 2 * rand();        % r3 in [0, 2], random weight of the destination
            r4 = rand();            % r4 in [0, 1], switches between the sine and cosine rules
            
            % Position update
            if r4 < 0.5
                % Sine rule
                X(i,j) = X(i,j) + r1 * sin(r2) * abs(r3 * Best_pos(j) - X(i,j));
            else
                % Cosine rule
                X(i,j) = X(i,j) + r1 * cos(r2) * abs(r3 * Best_pos(j) - X(i,j));
            end
        end
    end
    
    % Record the convergence curve
    Curve(t) = Best_score;
end

end