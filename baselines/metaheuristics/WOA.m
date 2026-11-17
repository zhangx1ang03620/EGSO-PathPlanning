function [Leader_score, Leader_pos, Curve] = WOA(num_population, max_gen, X_min, Xmax, dimensions, fobj, G)
% =========================================================================
% WOA: Whale Optimization Algorithm
% Used in the path-planning comparison experiments
% =========================================================================

% Initialise the leader (best solution)
Leader_pos = zeros(1, dimensions);
Leader_score = inf; 

% Convergence curve
Curve = zeros(1, max_gen);

% Initialise the population
Positions = rand(num_population, dimensions) .* (Xmax - X_min) + X_min;

t = 1; % current iteration
while t <= max_gen
    
    for i = 1:size(Positions, 1)
        
        % Boundary check (keep positions within [X_min, Xmax])
        % Equivalent: Positions(i,:) = Bounds(Positions(i,:), X_min, Xmax);
        Flag4ub = Positions(i,:) > Xmax;
        Flag4lb = Positions(i,:) < X_min;
        Positions(i,:) = (Positions(i,:) .* (~(Flag4ub + Flag4lb))) + Xmax .* Flag4ub + X_min .* Flag4lb;
        
        % Fitness
        fitness = fobj(Positions(i,:));
        
        % Update the leader
        if fitness < Leader_score
            Leader_score = fitness; 
            Leader_pos = Positions(i,:);
        end
    end
    
    % a decreases linearly from 2 to 0
    a = 2 - t * ((2) / max_gen);
    
    % a2 decreases linearly from -1 to -2 (used for l)
    a2 = -1 + t * ((-1) / max_gen);
    
    % Update every whale
    for i = 1:size(Positions, 1)
        r1 = rand(); % r1 in [0, 1]
        r2 = rand(); % r2 in [0, 1]
        
        A = 2 * a * r1 - a;  % Eq. (2.3)
        C = 2 * r2;          % Eq. (2.4)
        
        b = 1;               % spiral shape constant
        l = (a2 - 1) * rand() + 1;   % l in [a2, 1]
        
        p = rand();          % p in [0, 1]
        
        for j = 1:size(Positions, 2)
            if p < 0.5   
                if abs(A) >= 1
                    % Search for prey (exploration)
                    rand_leader_index = floor(num_population * rand() + 1);
                    X_rand = Positions(rand_leader_index, :);
                    D_X_rand = abs(C * X_rand(j) - Positions(i,j)); % Eq. (2.7)
                    Positions(i,j) = X_rand(j) - A * D_X_rand;      % Eq. (2.8)
                    
                elseif abs(A) < 1
                    % Encircling prey (exploitation)
                    D_Leader = abs(C * Leader_pos(j) - Positions(i,j)); % Eq. (2.1)
                    Positions(i,j) = Leader_pos(j) - A * D_Leader;      % Eq. (2.2)
                end
                
            elseif p >= 0.5
                % Spiral bubble-net attacking
                distance2Leader = abs(Leader_pos(j) - Positions(i,j));
                % Eq. (2.5)
                Positions(i,j) = distance2Leader * exp(b .* l) .* cos(l .* 2 .* pi) + Leader_pos(j);
            end
        end
    end
    
    Curve(t) = Leader_score;
    t = t + 1;
end

end