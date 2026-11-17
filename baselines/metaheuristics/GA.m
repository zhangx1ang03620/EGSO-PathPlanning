
%_________________________________________________________________________%
function [Best_score,Best_pos,curve]=GA(pop,Max_iter,lb,ub,dim,fobj,G)
%% Parameters
popsize=pop;              % population size
lenchrom=dim;              % chromosome length
fun = fobj;  % fitness function
pc=0.8;                  % crossover probability
pm=0.05;                  % mutation probability
if(max(size(ub)) == 1)
   ub = ub.*ones(dim,1);
   lb = lb.*ones(dim,1);  
end
maxgen=Max_iter;   % number of generations
% population
bound=[lb,ub];  % variable range

%% Initial population
for i=1:popsize
    % random initial population
    GApop(i,:)=Code(lenchrom,bound,G);       % random individual
    % fitness
    [fitness(i)]=fun(GApop(i,:));            % chromosome fitness
end

% best chromosome
[bestfitness bestindex]=min(fitness);
zbest=GApop(bestindex,:);   % global best
gbest=GApop;                % personal best
fitnessgbest=fitness;       % personal best fitness
fitnesszbest=bestfitness;   % global best fitness
%% Iterations
for i=1:maxgen
%         disp(['iteration ',num2str(i)])
        % GA selection
        GApop=Select2(GApop,fitness,popsize);

        % GA crossover
        GApop=Cross(pc,lenchrom,GApop,popsize,bound);

        % GA mutation
        GApop=Mutation(pm,lenchrom,GApop,popsize,[i maxgen],bound);

        pop=GApop;
        
      for j=1:popsize
        % fitness
        [fitness(j)]=fun(pop(j,:));
        % personal best update
        if fitness(j) < fitnessgbest(j)
            gbest(j,:) = pop(j,:);
            fitnessgbest(j) = fitness(j);
        end
        
        % global best update
        if fitness(j) < fitnesszbest
            zbest = pop(j,:);
            fitnesszbest = fitness(j);
        end
        
    end

    zbest = LocalSearch(zbest,ub(1),G); % greedy local search on the global best
    fitnesszbest = fun(zbest);         
    curve(i)=fitnesszbest;  
end
Best_score = fitnesszbest;
Best_pos = zbest;




end
%% Selection
function ret=Select2(individuals,fitness,sizepop)
% Roulette-wheel selection of the chromosomes for crossover and mutation
% individuals input  : population
% fitness input  : fitness
% sizepop     input  : population size
% opts        input  : selection method
% ret         output : selected population

fitness= 1./(fitness);
sumfitness=sum(fitness);
sumf=fitness./sumfitness;
index=[];
for i=1:sizepop   % spin the wheel sizepop times
    pick=rand;
    while pick==0
        pick=rand;
    end
    for j=1:sizepop
        pick=pick-sumf(j);
        if pick<0
            index=[index j];
            break;  % interval found: chromosome i is selected (a chromosome may be selected more than once)
        end
    end
end
individualsTemp=individuals(index,:);
fitnessTemp=fitness(index);
if(size(individualsTemp,1) == 0)
    ret=individuals;
else
    ret=individualsTemp;
end
end
%% Crossover
function ret=Mutation(pmutation,lenchrom,chrom,sizepop,pop,bound)
% Mutation operator
% pcorss                input  : mutation probability
% lenchrom              input  : chromosome length
% chrom                 input  : chromosomes
% sizepop               input  : population size
% pop                   input  : current and maximum generation
% ret                   output : mutated chromosomes

for i=1:sizepop  
    % Pick a random chromosome to mutate
    pick=rand;
    while pick==0
        pick=rand;
    end
    index=ceil(pick*sizepop);
    % The mutation probability decides whether this iteration mutates
    pick=rand;
    if pick>pmutation
        continue;
    end
    flag=0;
    while flag==0
        % Mutation position
        pick=rand;
        while pick==0
            pick=rand;
        end
        pos=ceil(pick*sum(lenchrom));  % random mutation position, i.e. variable pos is mutated
        if pos<=0 
            pos = 1;
        end
        if pos>size(bound,1)
            pos = size(bound,1);
        end
        v=chrom(i,pos);
        v1=v-bound(pos,1);
        v2=bound(pos,2)-v;
        pick=rand; % start of mutation
        if pick>0.5
            delta=v2*(1-pick^((1-pop(1)/pop(2))^2));
            chrom(i,pos)=v+delta;
        else
            delta=v1*(1-pick^((1-pop(1)/pop(2))^2));
            chrom(i,pos)=v-delta;
        end   % end of mutation
        flag=test(lenchrom,bound,chrom(i,:));     % feasibility check of the chromosome
    end
end
ret=chrom;
end

%% Mutation
function ret=Cross(pcross,lenchrom,chrom,sizepop,bound)
% Crossover operator
% pcorss                input  : crossover probability
% lenchrom              input  : chromosome length
% chrom                 input  : chromosomes
% sizepop               input  : population size
% ret                   output : chromosomes after crossover

for i=1:sizepop 
    
    % Pick two random chromosomes to cross
    pick=rand(1,2);
    while prod(pick)==0
        pick=rand(1,2);
    end
    index=ceil(pick.*sizepop);
    % The crossover probability decides whether to cross
    pick=rand;
    while pick==0
        pick=rand;
    end
    if pick>pcross
        continue;
    end
    flag=0;
    while flag==0
        % Random crossover position
        pick=rand;
        while pick==0
            pick=rand;
        end
        pos=ceil(pick.*sum(lenchrom)); % random crossover position (same variable index in both chromosomes)
        pick=rand; % start of crossover
        v1=chrom(index(1),pos);
        v2=chrom(index(2),pos);
        chrom(index(1),pos)=pick*v2+(1-pick)*v1;
        chrom(index(2),pos)=pick*v1+(1-pick)*v2; % end of crossover
        flag1=test(lenchrom,bound,chrom(index(1),:));  % feasibility of chromosome 1
        flag2=test(lenchrom,bound,chrom(index(2),:));  % feasibility of chromosome 2
        if   flag1*flag2==0
            flag=0;
        else flag=1;
        end    % repeat the crossover unless both chromosomes are feasible
    end
end
ret=chrom;
end
%% Range check
function flag=test(lenchrom,bound,code)
% lenchrom   input : chromosome length
% bound      input : variable range
% code       output: encoded chromosome

flag=1;
[n,m]=size(code);

for i=1:n
    if code(i)<bound(i,1) || code(i)>bound(i,2)
        flag=0;
    end
end
end
%% Encoding
function ret=Code(lenchrom,bound,G)
% Encodes the variables as a chromosome; used to initialise a random population
% lenchrom   input : chromosome length
% bound      input : variable range
% ret        output: encoded chromosome

flag=0;
while flag==0
    % pick=rand(1,lenchrom);
    % ret=bound(:,1)'+(bound(:,2)-bound(:,1))'.*pick; % linear interpolation
    for j=1:lenchrom
       column = G(:,j+1);      % one column of the map
       id = find(column == 0); % free cells of this column
       ret(1,j) =  id(randi(length(id))); % pick a random free cell
    end

    flag=test(lenchrom,bound,ret);             % feasibility check of the chromosome
end
end
