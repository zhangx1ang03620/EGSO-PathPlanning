function [Best_score, Best_pos, Curve] = TLBO(num_population, max_gen, X_min, Xmax, dimensions, fobj, G)
% =========================================================================
% TLBO: Teaching-Learning-Based Optimization
% Used in the path-planning comparison experiments
% =========================================================================

% Initialise the population (learners)
X = rand(num_population, dimensions) .* (Xmax - X_min) + X_min;
Fitness = zeros(num_population, 1);
Curve = zeros(1, max_gen);

% Initial evaluation
for i = 1:num_population
    Fitness(i) = fobj(X(i,:));
end

% Initial best solution (teacher)
[Best_score, idx] = min(Fitness);
Best_pos = X(idx, :);

% Main loop
for t = 1:max_gen
    
    % ===================== 1. Teacher phase =====================
    % Mean of the class in every dimension
    Mean_result = mean(X);
    
    % The best individual is the teacher
    [~, teacher_idx] = min(Fitness);
    Teacher = X(teacher_idx, :);
    
    for i = 1:num_population
        % Teaching factor TF, randomly 1 or 2
        TF = round(1 + rand()); 
        
        % Difference mean: the teacher pulls the class towards itself
        r = rand(1, dimensions);
        Difference_Mean = r .* (Teacher - TF * Mean_result);
        
        % Update the learner
        X_new = X(i,:) + Difference_Mean;
        
        % Boundary check
        X_new = min(max(X_new, X_min), Xmax);
        
        % Evaluate and accept if better (greedy)
        f_new = fobj(X_new);
        if f_new < Fitness(i)
            X(i,:) = X_new;
            Fitness(i) = f_new;
        end
    end
    
    % ===================== 2. Learner phase =====================
    for i = 1:num_population
        % Pick another learner at random (not itself)
        k = randi(num_population);
        while k == i
            k = randi(num_population);
        end
        
        % Mutual learning update
        r = rand(1, dimensions);
        if Fitness(i) < Fitness(k)
            % If the learner is better than its peer, move away from the peer
            X_new = X(i,:) + r .* (X(i,:) - X(k,:));
        else
            % If the peer is better, move towards the peer
            X_new = X(i,:) + r .* (X(k,:) - X(i,:));
        end
        
        % Boundary check
        X_new = min(max(X_new, X_min), Xmax);
        
        % Evaluate and accept if better (greedy)
        f_new = fobj(X_new);
        if f_new < Fitness(i)
            X(i,:) = X_new;
            Fitness(i) = f_new;
        end
    end
    
    % ===================== Record and update the best solution =====================
    [current_best, best_idx] = min(Fitness);
    if current_best < Best_score
        Best_score = current_best;
        Best_pos = X(best_idx, :);
    end
    
    % Best value of this generation
    Curve(t) = Best_score;
end

end