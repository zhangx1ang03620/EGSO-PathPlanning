function G = GenerateEGSO3DMap(R, C, H, seed)
% GenerateEGSO3DMap  Minimal 3D voxel warehouse-airspace map: a proof-of-
% concept scene for testing whether the EGSO design paradigm (auxiliary
% viscosity field threaded through the whole swarm pipeline) transfers
% beyond 2D grid path planning, as discussed (but never demonstrated) in
% the Discussion "Application Prospects" section.
%
%   R,C,H : grid size along row / column(=traversal axis) / height
%   seed  : rng seed
%
%   Floating obstacle blocks force the path to detour both horizontally
%   (row) and vertically (height) -- a genuine 3D navigation problem, not
%   a flat 2D slice embedded in 3D.

    rng(seed);
    G = zeros(R, C, H);

    num_blocks = max(2, round(C / 4));
    col_positions = round(linspace(5, C-4, num_blocks));
    for k = 1:num_blocks
        c0 = col_positions(k);
        c1 = min(C, col_positions(k)+1);
        if mod(k,2) == 1
            r0 = randi([1, max(1,R-8)]); r1 = min(R, r0+randi([4,7]));
            h0 = 1; h1 = min(H, randi([2, max(2,H-2)]));
        else
            r0 = randi([1, max(1,R-8)]); r1 = min(R, r0+randi([4,7]));
            h0 = max(1, randi([2, max(2,H-2)])); h1 = H;
        end
        G(r0:r1, c0:c1, h0:h1) = 1;
    end

    % light random clutter, kept off the start/goal safety zones
    for t = 1:round(R*C*H*0.02)
        r = randi(R); c = randi([3, C-2]); h = randi(H);
        if norm([r,c,h]-[1,1,1]) > 4 && norm([r,c,h]-[R,C,H]) > 4
            G(r,c,h) = 1;
        end
    end

    G(1:min(3,R), 1:min(3,C), 1:min(3,H)) = 0;
    G(max(1,R-2):R, max(1,C-2):C, max(1,H-2):H) = 0;
end
