function G = GenerateScalabilityMap(M, seed)
% GenerateScalabilityMap  Parameterized dense warehouse-shelf map (same
% construction family as Map5/MapE in GenerateEGSOMaps.m), scaled to an
% arbitrary size M, used to test how EVF precomputation cost and overall
% EGSO runtime actually grow with problem scale (Sec. Discussion,
% "limitation 3": the O((M*M)^2) worst-case precompute cost).

    rng(seed);
    G = zeros(M, M);

    shelf_cols = 10:8:M-10;
    for sc = shelf_cols
        G(6:M-5, sc:sc+2) = 1;
    end
    aisle_rows = round(linspace(20, M-20, 3));
    for ar = aisle_rows
        G(ar:ar+3, :) = 0;
    end
    for si = 1:numel(shelf_cols)
        sc = shelf_cols(si);
        extra_gap = 10 + mod(si*17, max(1, M-20));
        extra_gap = min(extra_gap, M-4);
        G(extra_gap:extra_gap+3, sc:sc+2) = 0;
    end

    for k = 1:round(M * 1.67)
        r = randi([3, M-2]); c = randi([3, M-2]);
        if norm([r,c]-[1,1]) > 12 && norm([r,c]-[M,M]) > 12
            G(r, c) = 1;
        end
    end

    G(1:8,   1:8)   = 0;
    G(M-7:M, M-7:M) = 0;
end
