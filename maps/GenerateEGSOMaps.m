function [Maps, Map_Info] = GenerateEGSOMaps()
% =========================================================================
% Generate the 5 benchmark maps of the EGSO paper. Two of them (Map A, Map B) are purpose-built "decoy dead-end"
% traps: a wall gap that looks attractive (near the straight start-goal
% line) leads into a pocket that is sealed a few columns later, while the
% true through-gap lies far from the straight line. A purely local
% attractive/repulsive potential field (evaluated independently per
% column) cannot see that the near gap is a dead end, because that
% requires GLOBAL connectivity information -- exactly what the
% Eikonal viscosity field provides and a hand-crafted APF does not.
%
%   Map A -- 50x50  Single decoy-channel dead-end trap
%   Map B -- 50x50  Multi-stage (3x) decoy corridor, increasing difficulty
%   Map C -- 60x60  Procedural Voronoi room/corridor layout
%   Map D -- 80x80  Large-scale sparse long-detour network
%   Map E -- 90x90  Dense warehouse shelf grid (high density, large scale)
%
% Encoding note: all comparison algorithms in this project represent a
% path as one row-index per column (dim = M-2), i.e. the path is a
% single-valued function of the column index. Obstacle layouts are
% therefore designed as column-wise wall/gap structures (never requiring
% the path to revisit a column), which is the representational
% family of grid path-planning benchmarks but with a fresh visual/topological
% design here.
% =========================================================================

Maps     = cell(5, 1);
Map_Info = struct();

%% ========================================================================
%  Map A: 50x50  Single decoy-channel dead-end trap
%% ========================================================================
fprintf('Generating Map A: 50x50  Decoy-channel dead-end trap...\n');
rng(31101);
M = 50;
MapA = zeros(M, M);

% light background clutter (kept away from the engineered decoy geometry)
for k = 1:70
    r = randi([3, M-2]); c = randi([3, 16]);
    MapA(r, c) = 1;
end
for k = 1:70
    r = randi([3, M-2]); c = randi([26, M-2]);
    MapA(r, c) = 1;
end

% --- decoy wall (cols 18-19): near-diagonal gap (attractive) + far gap (real) ---
MapA(:, 18:19)  = 1;
MapA(15:21, 18:19) = 0;   % decoy gap: sits right on the start-goal diagonal (~row 18)
MapA(2:7,   18:19) = 0;   % real gap: far from the diagonal, near the top border

% --- closure wall (cols 22-23): seals the decoy channel into a dead-end pocket ---
MapA(:, 22:23) = 1;
MapA(2:7, 22:23) = 0;     % only the real channel continues past the closure

% safety zones
MapA(1:5,   1:5)   = 0;
MapA(M-4:M, M-4:M) = 0;

Maps{1}             = MapA;
Map_Info(1).size    = '50x50';
Map_Info(1).density = sprintf('%.1f%%', 100*sum(MapA(:))/numel(MapA));
Map_Info(1).pattern = 'Single decoy-channel dead-end trap';

%% ========================================================================
%  Map B: 50x50  Multi-stage (3x) decoy corridor
%% ========================================================================
fprintf('Generating Map B: 50x50  Multi-stage decoy corridor...\n');
rng(31102);
M = 50;
MapB = zeros(M, M);

for k = 1:60
    r = randi([3, M-2]); c = randi([3, M-2]);
    if c < 10 || (c > 16 && c < 24) || (c > 30 && c < 38) || c > 44
        MapB(r, c) = 1;
    end
end

% --- stage 1 (cols 10-11): decoy near diagonal (~row10), real gap far (near bottom) ---
MapB(:, 10:11) = 1;
MapB(7:13,  10:11) = 0;   % decoy
MapB(40:46, 10:11) = 0;   % real
MapB(:, 14:15) = 1;
MapB(40:46, 14:15) = 0;   % closure seals decoy stage-1 pocket

% --- stage 2 (cols 24-25): decoy near diagonal (~row24), real gap far (near top) ---
MapB(:, 24:25) = 1;
MapB(21:27, 24:25) = 0;   % decoy
MapB(2:8,   24:25) = 0;   % real
MapB(:, 28:29) = 1;
MapB(2:8, 28:29) = 0;     % closure seals decoy stage-2 pocket

% --- stage 3 (cols 38-39): decoy near diagonal (~row38), real gap far (near bottom) ---
MapB(:, 38:39) = 1;
MapB(35:41, 38:39) = 0;   % decoy
MapB(44:49, 38:39) = 0;   % real
MapB(:, 42:43) = 1;
MapB(44:49, 42:43) = 0;   % closure seals decoy stage-3 pocket

MapB(1:5,   1:5)   = 0;
MapB(M-4:M, M-4:M) = 0;

Maps{2}             = MapB;
Map_Info(2).size    = '50x50';
Map_Info(2).density = sprintf('%.1f%%', 100*sum(MapB(:))/numel(MapB));
Map_Info(2).pattern = 'Multi-stage (3x) decoy corridor';

%% ========================================================================
%  Map C: 60x60  Procedural Voronoi room/corridor layout
%% ========================================================================
fprintf('Generating Map C: 60x60  Voronoi room/corridor layout...\n');
rng(31103);
M = 60;
MapC = zeros(M, M);

num_seeds = 14;
seeds = [randi([4,M-3], num_seeds, 1), randi([4,M-3], num_seeds, 1)];
[cc, rr] = meshgrid(1:M, 1:M);
regionId = zeros(M, M);
for r = 1:M
    for c = 1:M
        d2 = (seeds(:,1)-r).^2 + (seeds(:,2)-c).^2;
        [~, idx] = min(d2);
        regionId(r,c) = idx;
    end
end

% walls along region boundaries
for r = 1:M-1
    for c = 1:M-1
        if regionId(r,c) ~= regionId(r+1,c) || regionId(r,c) ~= regionId(r,c+1)
            MapC(r,c) = 1;
        end
    end
end

% punch a door gap (width 3) through every boundary wall segment that
% still separates two adjacent regions, so rooms remain traversable
for r = 3:M-2
    for c = 3:M-2
        if MapC(r,c) == 1 && mod(r + 2*c, 9) == 0
            MapC(max(1,r-1):min(M,r+1), c) = 0;
        end
    end
end

MapC(1:5,   1:5)   = 0;
MapC(M-4:M, M-4:M) = 0;

Maps{3}             = MapC;
Map_Info(3).size    = '60x60';
Map_Info(3).density = sprintf('%.1f%%', 100*sum(MapC(:))/numel(MapC));
Map_Info(3).pattern = 'Procedural Voronoi room/corridor layout';

%% ========================================================================
%  Map D: 80x80  Large-scale sparse long-detour network
%% ========================================================================
fprintf('Generating Map D: 80x80  Sparse long-detour network...\n');
rng(31104);
M = 80;
MapD = zeros(M, M);

% a few moderate sparse obstacles that force long detours (kept well
% under half the map area so the map stays genuinely "sparse")
big_blocks = [10 10 20 26; 15 45 26 61; 40 5 51 21; 45 40 59 56; 20 55 31 69];
for k = 1:size(big_blocks,1)
    r0=big_blocks(k,1); c0=big_blocks(k,2); r1=big_blocks(k,3); c1=big_blocks(k,4);
    MapD(r0:r1, c0:c1) = 1;
end
% carve a single-corridor gap through each block so it's a detour, not a
% full barrier
gaps = [15 18; 20 53; 45 13; 52 48; 25 62];
for k = 1:size(gaps,1)
    MapD(max(1,gaps(k,1)-2):min(M,gaps(k,1)+2), gaps(k,2)) = 0;
end

for k = 1:60
    r = randi([3, M-2]); c = randi([3, M-2]);
    if norm([r,c]-[1,1]) > 10 && norm([r,c]-[M,M]) > 10
        MapD(r, c) = 1;
    end
end

MapD(1:8,   1:8)   = 0;
MapD(M-7:M, M-7:M) = 0;

Maps{4}             = MapD;
Map_Info(4).size    = '80x80';
Map_Info(4).density = sprintf('%.1f%%', 100*sum(MapD(:))/numel(MapD));
Map_Info(4).pattern = 'Sparse large-scale long-detour network';

%% ========================================================================
%  Map E: 90x90  Dense warehouse shelf grid
%% ========================================================================
fprintf('Generating Map E: 90x90  Dense warehouse shelf grid...\n');
rng(31105);
M = 90;
MapE = zeros(M, M);

shelf_cols = 10:8:M-10;
for sc = shelf_cols
    MapE(6:M-5, sc:sc+2) = 1;
end
% cross-aisles: horizontal gaps through all shelves at a few rows
aisle_rows = [20 45 70];
for ar = aisle_rows
    MapE(ar:ar+3, :) = 0;
end
% stagger extra gaps so no single aisle row spans the full map identically
for si = 1:numel(shelf_cols)
    sc = shelf_cols(si);
    extra_gap = 10 + mod(si*17, 60);
    MapE(extra_gap:extra_gap+3, sc:sc+2) = 0;
end

for k = 1:150
    r = randi([3, M-2]); c = randi([3, M-2]);
    if norm([r,c]-[1,1]) > 12 && norm([r,c]-[M,M]) > 12
        MapE(r, c) = 1;
    end
end

MapE(1:8,   1:8)   = 0;
MapE(M-7:M, M-7:M) = 0;

Maps{5}             = MapE;
Map_Info(5).size    = '90x90';
Map_Info(5).density = sprintf('%.1f%%', 100*sum(MapE(:))/numel(MapE));
Map_Info(5).pattern = 'Dense warehouse shelf grid';

%% ========================================================================
%  Connectivity verification (uses the same FMM solver as EGSO itself):
%  every map must have a finite T_nav at the start cell, i.e. a feasible
%  collision-free path from S=(1,1) to E=(M,M) must exist.
%% ========================================================================
fprintf('\nVerifying feasibility (S -> E) via Eikonal field...\n');
for i = 1:5
    G = Maps{i};
    Mi = size(G,1);
    goalMask = false(Mi,Mi); goalMask(Mi,Mi) = true;
    T = FastMarchingField(G, goalMask, 'navigate');
    feasible = isfinite(T(1,1));
    fprintf('  Map %d (%s): S->E feasible = %d, T(1,1) = %.2f, density = %s\n', ...
            i, Map_Info(i).size, feasible, T(1,1), Map_Info(i).density);
    if ~feasible
        error('GenerateEGSOMaps:Infeasible', 'Map %d has no feasible path from S to E.', i);
    end
end

%% ========================================================================
%  Visualisation (publication quality)
%% ========================================================================
figure('Name', 'EGSO Test Maps', 'Position', [80, 80, 1700, 950], 'Color', 'w');

panel_titles = { ...
    'Map A: 50\times50  Decoy-channel trap', ...
    'Map B: 50\times50  Multi-stage decoy corridor', ...
    'Map C: 60\times60  Voronoi rooms', ...
    'Map D: 80\times80  Sparse long-detour', ...
    'Map E: 90\times90  Warehouse shelf grid'};

pos_list = { ...
    [0.04, 0.54, 0.28, 0.36], ...
    [0.36, 0.54, 0.28, 0.36], ...
    [0.68, 0.54, 0.28, 0.36], ...
    [0.12, 0.05, 0.34, 0.36], ...
    [0.54, 0.05, 0.34, 0.36]};

obstacle_color = [0.20 0.45 0.70];
free_color     = [0.97 0.97 0.97];

for i = 1:5
    ax = axes('Position', pos_list{i}); %#ok<LAXES>
    Mi = Maps{i};
    ms = size(Mi, 1);

    imagesc(ax, Mi);
    colormap(ax, [free_color; obstacle_color]);
    clim(ax, [0 1]);
    axis(ax, 'equal', 'tight');
    set(ax, 'XDir','normal', 'YDir','normal', 'Box','on', 'FontSize',8);
    grid(ax, 'on');
    hold(ax, 'on');

    plot(ax, 1,  1,  'o', 'MarkerSize',9, 'MarkerFaceColor',[0.18 0.70 0.22], ...
         'MarkerEdgeColor','k', 'LineWidth',1.2);
    plot(ax, ms, ms, 's', 'MarkerSize',9, 'MarkerFaceColor',[0.85 0.12 0.12], ...
         'MarkerEdgeColor','k', 'LineWidth',1.2);

    title(ax, panel_titles{i}, 'FontSize',10, 'FontWeight','bold', 'Interpreter','tex');
end

sgtitle('EGSO: Benchmark Test Maps', 'FontSize',14, 'FontWeight','bold');

%% Summary report
fprintf('\n%-6s  %-10s  %-10s  %-40s\n', 'Map','Grid','Density','Pattern');
labels = {'A','B','C','D','E'};
for i = 1:5
    fprintf('Map %-3s %-10s  %-10s  %-40s\n', labels{i}, Map_Info(i).size, Map_Info(i).density, Map_Info(i).pattern);
end
fprintf('\nAll maps generated and verified feasible.\n');

end
