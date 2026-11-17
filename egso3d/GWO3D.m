function [Alpha_score, Alpha_pos, Convergence_curve] = GWO3D(SearchAgents_no, Max_iter, dim, fobj, G)
% GWO3D  Minimal 3D port of GWO.m (same mechanism, decision vector length
% 2*dim: row indices then height indices), used as a classical swarm
% baseline in the 3D proof-of-concept to show the swarm-optimization
% layer still adds value over the field alone (paralleling FMMDescentPath3D).

[R, C, H] = size(G);

Alpha_pos = zeros(1,2*dim); Alpha_score = inf;
Beta_pos  = zeros(1,2*dim); Beta_score  = inf;
Delta_pos = zeros(1,2*dim); Delta_score = inf;

Positions = zeros(SearchAgents_no, 2*dim);
for i = 1:SearchAgents_no
    for j = 1:dim
        col = j+1;
        slice = squeeze(G(:,col,:) == 0);
        [fr, fh] = find(slice);
        if isempty(fr), fr=1; fh=1; end
        pick = randi(numel(fr));
        Positions(i,j)     = fr(pick);
        Positions(i,dim+j) = fh(pick);
    end
end

Convergence_curve = zeros(1, Max_iter);
lb = [ones(1,dim), ones(1,dim)];
ub = [R*ones(1,dim), H*ones(1,dim)];

l = 0;
while l < Max_iter
    for i = 1:SearchAgents_no
        Flag4ub = Positions(i,:) > ub;
        Flag4lb = Positions(i,:) < lb;
        Positions(i,:) = (Positions(i,:).*(~(Flag4ub+Flag4lb))) + ub.*Flag4ub + lb.*Flag4lb;

        fit = fobj(Positions(i,:));

        if fit < Alpha_score
            Alpha_score = fit; Alpha_pos = Positions(i,:);
        elseif fit < Beta_score
            Beta_score = fit; Beta_pos = Positions(i,:);
        elseif fit < Delta_score
            Delta_score = fit; Delta_pos = Positions(i,:);
        end
    end

    a = 2 - l*(2/Max_iter);

    for i = 1:SearchAgents_no
        for j = 1:2*dim
            r1=rand; r2=rand; A1=2*a*r1-a; C1=2*r2;
            D_alpha = abs(C1*Alpha_pos(j)-Positions(i,j)); X1 = Alpha_pos(j)-A1*D_alpha;
            r1=rand; r2=rand; A2=2*a*r1-a; C2=2*r2;
            D_beta = abs(C2*Beta_pos(j)-Positions(i,j)); X2 = Beta_pos(j)-A2*D_beta;
            r1=rand; r2=rand; A3=2*a*r1-a; C3=2*r2;
            D_delta = abs(C3*Delta_pos(j)-Positions(i,j)); X3 = Delta_pos(j)-A3*D_delta;
            Positions(i,j) = (X1+X2+X3)/3;
        end
    end
    l = l+1;
    Convergence_curve(l) = Alpha_score;
end

end
