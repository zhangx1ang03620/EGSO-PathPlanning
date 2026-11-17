function [Rabbit_Energy, Rabbit_Location, Curve] = HHO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G)
% =========================================================================
% HHO: Harris Hawks Optimization
% Used in the path-planning comparison experiments
% =========================================================================

% Initialise parameters
Rabbit_Location = zeros(1, dimensions);
Rabbit_Energy = inf; 
Curve = zeros(1, max_gen);

% Initialise the population
X = rand(num_population, dimensions) .* (Xmax - X_min) + X_min;

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
        
        % Update the prey (best solution)
        if fitness < Rabbit_Energy
            Rabbit_Energy = fitness;
            Rabbit_Location = X(i,:);
        end
    end
    
    % Escaping energy E of the prey
    E1 = 2 * (1 - (t / max_gen)); % E1 decreases linearly from 2 to 0
    
    % Update every hawk
    for i = 1:num_population
        
        E0 = 2 * rand() - 1; % random number in [-1, 1]
        Escaping_Energy = E1 * E0; % escaping energy E
        
        % Exploration phase
        if abs(Escaping_Energy) >= 1
            q = rand();
            rand_Hawk_index = floor(num_population * rand() + 1);
            X_rand = X(rand_Hawk_index, :);
            
            if q >= 0.5
                % Update from a random hawk
                X(i,:) = X_rand - rand() * abs(X_rand - 2 * rand() * X(i,:));
            elseif q < 0.5
                % Update from the population mean with random scaling
                X_m = mean(X); % mean position of the population
                X(i,:) = (Rabbit_Location - X_m) - rand() * (X_min + rand() * (Xmax - X_min));
            end
            
        % Exploitation phase
        elseif abs(Escaping_Energy) < 1
            r = rand(); % chance of a successful escape
            
            % Strategy 1: soft besiege
            if r >= 0.5 && abs(Escaping_Energy) >= 0.5
                J = 2 * (1 - rand()); % random jump strength
                X(i,:) = (Rabbit_Location - X(i,:)) - Escaping_Energy * abs(J * Rabbit_Location - X(i,:));
                
            % Strategy 2: hard besiege
            elseif r >= 0.5 && abs(Escaping_Energy) < 0.5
                X(i,:) = Rabbit_Location - Escaping_Energy * abs(Rabbit_Location - X(i,:));
                
            % Strategy 3: soft besiege with progressive rapid dives
            elseif r < 0.5 && abs(Escaping_Energy) >= 0.5
                J = 2 * (1 - rand());
                Y = Rabbit_Location - Escaping_Energy * abs(J * Rabbit_Location - X(i,:));
                
                if fobj(Y) < fobj(X(i,:))
                    X(i,:) = Y;
                else
                    S = rand(1, dimensions);
                    Z = Y + S .* Levy(dimensions);
                    if fobj(Z) < fobj(X(i,:))
                        X(i,:) = Z;
                    end
                end
                
            % Strategy 4: hard besiege with progressive rapid dives
            elseif r < 0.5 && abs(Escaping_Energy) < 0.5
                J = 2 * (1 - rand());
                X_m = mean(X);
                Y = Rabbit_Location - Escaping_Energy * abs(J * Rabbit_Location - X_m);
                
                if fobj(Y) < fobj(X(i,:))
                    X(i,:) = Y;
                else
                    S = rand(1, dimensions);
                    Z = Y + S .* Levy(dimensions);
                    if fobj(Z) < fobj(X(i,:))
                        X(i,:) = Z;
                    end
                end
            end
        end
    end
    
    % Record the convergence curve
    Curve(t) = Rabbit_Energy;
end

end

% =========================================================================
% Local function: Levy flight
% =========================================================================
function o = Levy(d)
    beta = 1.5;
    sigma = (gamma(1 + beta) * sin(pi * beta / 2) / (gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2)))^(1 / beta);
    u = randn(1, d) * sigma;
    v = randn(1, d);
    step = u ./ abs(v).^(1 / beta);
    o = step;
end