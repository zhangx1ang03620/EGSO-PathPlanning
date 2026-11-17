function G = GenerateDecoyStressMap(M, num_stages, seed)
% GenerateDecoyStressMap  Parameterized decoy dead-end stress-test map.
%
% Places `num_stages` sequential decoy-channel dead-end traps (the same
% construction used for Map A/B in GenerateEGSOMaps.m) along the diagonal
% from S=(1,1) to E=(M,M): each stage's wall has a "decoy" gap near the
% straight-line row (locally attractive to a hand-crafted potential field)
% and a "real" gap far from it; a closure wall a few columns later seals
% the decoy channel into a dead-end pocket. Increasing num_stages raises
% the local-minima complexity level used by Main_LocalMinima_StressTest.m.

rng(seed);
G = zeros(M, M);

usable_lo = 8;
usable_hi = M - 8;
span      = usable_hi - usable_lo;
spacing   = span / (num_stages + 1);

engineered_cols = [];
for k = 1:num_stages
    col_start = round(usable_lo + k * spacing);
    col_start = max(usable_lo, min(usable_hi - 6, col_start));
    wall_cols    = col_start : col_start+1;
    closure_cols = col_start+4 : col_start+5;
    engineered_cols = [engineered_cols, wall_cols, closure_cols]; %#ok<AGROW>

    diag_row = round(1 + (M-1) * (col_start / M));
    decoy_lo = max(1, diag_row - 3);
    decoy_hi = min(M, diag_row + 3);

    if mod(k, 2) == 1
        real_lo = 2; real_hi = min(M, 8);         % real gap near top
    else
        real_lo = max(1, M-7); real_hi = M-1;     % real gap near bottom
    end

    G(:, wall_cols) = 1;
    G(decoy_lo:decoy_hi, wall_cols) = 0;   % decoy gap (near diagonal)
    G(real_lo:real_hi,   wall_cols) = 0;   % real gap (far from diagonal)

    G(:, closure_cols) = 1;
    G(real_lo:real_hi, closure_cols) = 0;  % only the real channel continues
end

% light background clutter, kept off the engineered columns
for t = 1:round(M*1.2)
    r = randi([3, M-2]); c = randi([3, M-2]);
    if ~ismember(c, engineered_cols)
        G(r, c) = 1;
    end
end

G(1:5,   1:5)   = 0;
G(M-4:M, M-4:M) = 0;

end
