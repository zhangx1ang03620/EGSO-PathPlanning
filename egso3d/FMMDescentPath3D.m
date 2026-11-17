function [path, path_length] = FMMDescentPath3D(G, S, E)
% FMMDescentPath3D  3D analogue of FMMDescentPath.m: pure 6-connected
% steepest descent on the precomputed 3D Eikonal field, no swarm
% refinement, no CALS. Answers the same "what does the swarm layer add
% over the field alone?" question in the 3D proof-of-concept.

    [R, C, H] = size(G);
    goalMask = false(R,C,H); goalMask(E(1),E(2),E(3)) = true;
    T = FastMarchingField3D(G, goalMask, 'navigate');

    if isinf(T(S(1),S(2),S(3)))
        path = []; path_length = Inf; return;
    end

    nbrs6 = [-1 0 0; 1 0 0; 0 -1 0; 0 1 0; 0 0 -1; 0 0 1];
    path = [S(1), S(2), S(3)];
    r=S(1); c=S(2); h=S(3);
    visited = false(R,C,H); visited(r,c,h) = true;

    maxSteps = R*C*H*2;
    found = false;
    for step = 1:maxSteps
        if r==E(1) && c==E(2) && h==E(3)
            found = true; break;
        end
        best_val = Inf; best_r=-1; best_c=-1; best_h=-1;
        for k = 1:6
            nr=r+nbrs6(k,1); nc=c+nbrs6(k,2); nh=h+nbrs6(k,3);
            if nr<1||nr>R||nc<1||nc>C||nh<1||nh>H, continue; end
            if visited(nr,nc,nh), continue; end
            tval = T(nr,nc,nh);
            if tval < best_val
                best_val=tval; best_r=nr; best_c=nc; best_h=nh;
            end
        end
        if best_r < 0
            break;
        end
        r=best_r; c=best_c; h=best_h;
        visited(r,c,h) = true;
        path = [path; r c h]; %#ok<AGROW>
    end

    if ~found
        path = []; path_length = Inf; return;
    end

    path_length = 0;
    for i = 1:size(path,1)-1
        path_length = path_length + norm(path(i+1,:)-path(i,:));
    end
end
