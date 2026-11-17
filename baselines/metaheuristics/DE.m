function [bestf,bestx,BestCost]  = DE(nPop,MaxIt,VarMin,VarMax,nVar,CostFunction,G)
%% DE Parameters

% nVar= Number of Decision Variables
VarSize=[1 nVar];   % Decision Variables Matrix Size
% VarMin=Lower Bound of Decision Variables
% VarMax= Upper Bound of Decision Variables
% MaxIt= Maximum Number of Iterations
% nPop= Population Size
beta_min=0.2;   % Lower Bound of Scaling Factor
beta_max=0.8;   % Upper Bound of Scaling Factor
pCR=0.2;        % Crossover Probability
%% Initialization
empty_individual.Position=[];
empty_individual.Cost=[];
BestSol.Cost=inf;
pop=repmat(empty_individual,nPop,1);
for i=1:nPop
    for j = 1:nVar
       column = G(:,j+1);      % one column of the map
       id = find(column == 0); % free cells of this column
       x(1,j) =  id(randi(length(id))); % pick a random free cell
       id = [];
    end
    pop(i).Position=x;
    pop(i).Cost=CostFunction(pop(i).Position);
    if pop(i).Cost<BestSol.Cost
        BestSol=pop(i);
    end
end

BestCost=zeros(MaxIt,1);
%% DE Main Loop
for it=1:MaxIt
    
    for i=1:nPop
        
        x=pop(i).Position;
        
        A=randperm(nPop);
        
        A(A==i)=[];
        
        a=A(1);
        b=A(2);
        c=A(3);
        
        % Mutation
        % ===== Replaced unifrnd (Statistics Toolbox) with base MATLAB =====
        % Original: beta=unifrnd(beta_min,beta_max,VarSize);
        beta = beta_min + (beta_max - beta_min) * rand(VarSize);
        % =================================================
        
        y=pop(a).Position+beta.*(pop(b).Position-pop(c).Position);
        y = max(y, VarMin);
		y = min(y, VarMax);
		
        % Crossover
        z=zeros(size(x));
        j0=randi([1 numel(x)]);
        for j=1:numel(x)
            if j==j0 || rand<=pCR
                z(j)=y(j);
            else
                z(j)=x(j);
            end
        end
        
        NewSol.Position=z;
        NewSol.Cost=CostFunction(NewSol.Position);
        
        if NewSol.Cost<pop(i).Cost
            pop(i)=NewSol;
            
            if pop(i).Cost<BestSol.Cost
               BestSol=pop(i);
            end
        end
        
    end
    
    % Update Best Cost
    
    BestSol.Position = LocalSearch(BestSol.Position,VarMax(1),G); % greedy local search on the global best
    BestSol.Cost = CostFunction(BestSol.Position);    
    BestCost(it)=BestSol.Cost;
end
%% Show Results

bestx = BestSol.Position;
bestf = BestCost(end);
end
