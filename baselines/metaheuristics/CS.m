function [Best_score, Best_pos, Curve] = CS(num_population, max_gen, X_min, Xmax, dimensions, fobj, G)
% =========================================================================
% CS: Cuckoo Search
% Used in the path-planning comparison experiments
% =========================================================================

% CS-specific parameters
pa = 0.25; % discovery probability of an alien egg

% Initialise the nests (population)
nest = rand(num_population, dimensions) .* (Xmax - X_min) + X_min;
fitness = zeros(num_population, 1);
Curve = zeros(1, max_gen);

% Initial evaluation
for i = 1:num_population
    fitness(i) = fobj(nest(i,:));
end

% Initial best solution
[fmin, K] = min(fitness);
Best_pos = nest(K, :);
Best_score = fmin;

% Main loop
for t = 1:max_gen
    
    % ===================== 1. New nests by Levy flight =====================
    new_nest = get_cuckoos(nest, Best_pos, X_min, Xmax);
    
    % Evaluate the new nests and replace the old ones if better (greedy)
    for i = 1:num_population
        fnew = fobj(new_nest(i,:));
        if fnew <= fitness(i)
            fitness(i) = fnew;
            nest(i,:) = new_nest(i,:);
        end
    end
    
    % Best solution so far
    [fmin, K] = min(fitness);
    if fmin < Best_score
        Best_score = fmin;
        Best_pos = nest(K, :);
    end
    
    % ===================== 2. Discover and abandon alien nests =====================
    % With probability pa a host discovers an alien egg and the nest is rebuilt at a new position
    new_nest = empty_nests(nest, X_min, Xmax, pa);
    
    % Keep the rebuilt nest if it is better
    for i = 1:num_population
        fnew = fobj(new_nest(i,:));
        if fnew <= fitness(i)
            fitness(i) = fnew;
            nest(i,:) = new_nest(i,:);
        end
    end
    
    % Update the global best
    [fmin, K] = min(fitness);
    if fmin < Best_score
        Best_score = fmin;
        Best_pos = nest(K, :);
    end
    
    % Best value of this generation
    Curve(t) = Best_score;
end

end

% =========================================================================
% Local function: new solutions by Levy flight
% =========================================================================
function nest = get_cuckoos(nest, best, Lb, Ub)
    n = size(nest, 1);
    % Levy flight constants
    beta = 1.5;
    sigma = (gamma(1 + beta) * sin(pi * beta / 2) / (gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2)))^(1 / beta);
    
    for j = 1:n
        s = nest(j, :);
        % Levy step (Mantegna's algorithm)
        u = randn(size(s)) * sigma;
        v = randn(size(s));
        step = u ./ abs(v).^(1 / beta);
        
        % Step scale, proportional to the distance between the current position and the best
        stepsize = 0.01 * step .* (s - best);
        
        % Update the position
        s = s + stepsize .* randn(size(s));
        
        % Boundary check
        nest(j, :) = min(max(s, Lb), Ub);
    end
end

% =========================================================================
% Local function: abandon discovered nests and build new ones locally
% =========================================================================
function new_nest = empty_nests(nest, Lb, Ub, pa)
    n = size(nest, 1);
    % Logical mask of the abandoned dimensions/nests (entries above pa are reset)
    K = rand(size(nest)) > pa;
    
    % The difference of two random permutations of the population guides the new nests (diversity)
    stepsize = rand(size(nest)) .* (nest(randperm(n), :) - nest(randperm(n), :));
    new_nest = nest + stepsize .* K;
    
    % Boundary check
    for j = 1:n
        new_nest(j, :) = min(max(new_nest(j, :), Lb), Ub);
    end
end