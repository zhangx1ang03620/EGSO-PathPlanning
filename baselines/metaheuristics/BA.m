function [Best_score, Best_pos, Curve] = BA(num_population, max_gen, X_min, Xmax, dimensions, fobj, G)
% =========================================================================
% BA: Bat Algorithm
% Used in the path-planning comparison experiments
% =========================================================================

% BA-specific parameters
Qmin = 0;         % minimum pulse frequency
Qmax = 2;         % maximum pulse frequency
alpha = 0.9;      % loudness decay factor
gamma = 0.9;      % pulse-rate increase factor

% Initialise the population
X = rand(num_population, dimensions) .* (Xmax - X_min) + X_min; % positions
v = zeros(num_population, dimensions);                          % velocities
A = ones(num_population, 1);                                    % loudness (initially 1)
r0 = ones(num_population, 1) * 0.5;                             % initial pulse emission rate
r = r0;                                                         % current pulse emission rate

Fitness = zeros(num_population, 1);
Curve = zeros(1, max_gen);

Best_pos = zeros(1, dimensions);
Best_score = inf;

% Evaluate the initial population
for i = 1:num_population
    Fitness(i) = fobj(X(i,:));
    if Fitness(i) < Best_score
        Best_score = Fitness(i);
        Best_pos = X(i,:);
    end
end

% Main loop
for t = 1:max_gen
    for i = 1:num_population
        
        % 1. Update frequency, velocity and position
        Q = Qmin + (Qmax - Qmin) * rand();
        v(i,:) = v(i,:) + (X(i,:) - Best_pos) * Q;
        X_new = X(i,:) + v(i,:);
        
        % 2. Local search
        % If the random number exceeds the pulse rate, take a local random walk around the best solution
        if rand() > r(i)
            % The mean loudness of all bats sets the local step size
            step_size = 0.01 * mean(A) * randn(1, dimensions);
            X_new = Best_pos + step_size;
        end
        
        % Boundary check
        X_new = min(max(X_new, X_min), Xmax);
        
        % Fitness of the new position
        f_new = fobj(X_new);
        
        % 3. Acceptance rule
        % Accept if the new solution is better and the random number is below the loudness A
        if (rand() < A(i)) && (f_new < Fitness(i))
            X(i,:) = X_new;
            Fitness(i) = f_new;
            
            % After a successful move, decrease the loudness and increase the pulse rate
            A(i) = alpha * A(i);
            r(i) = r0(i) * (1 - exp(-gamma * t));
        end
        
        % 4. Update the global best
        if Fitness(i) < Best_score
            Best_score = Fitness(i);
            Best_pos = X(i,:);
        end
    end
    
    % Best value of this generation
    Curve(t) = Best_score;
end

end