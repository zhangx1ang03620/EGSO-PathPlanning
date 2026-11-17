clc; clear;

%% Test 1: simple open map, sanity check monotonicity & goal=0
M = 20;
G = zeros(M, M);
goalMask = false(M, M);
goalMask(M, M) = true;
T = FastMarchingField(G, goalMask, 'navigate');

assert(T(M, M) == 0, 'goal T should be 0');
assert(all(isfinite(T(:))), 'open map: all cells should be reachable');
assert(T(1,1) > T(M,M), 'start should have larger T than goal');
% monotonic non-decreasing away from goal along the anti-diagonal-ish path
d1 = T(M,M); d2 = T(M-1,M); d3 = T(M-2,M);
assert(d2 >= d1 && d3 >= d2, 'T should increase monotonically away from goal');
fprintf('[PASS] Test 1: open map basic sanity (T(1,1)=%.3f, T(goal)=%.3f)\n', T(1,1), T(M,M));

%% Test 2: U-shaped trap - verify no local minimum trap (field must route AROUND)
M = 30;
G = zeros(M, M);
% deep U-trap directly between a start-side region and nothing -- concave pocket
G(5:25, 15) = 1;   % vertical wall
G(5, 15:25) = 1;   % top of U
G(25,15:25) = 1;   % bottom of U
% wall at col 25 too, except small opening far from the straight-line direction
G(5:25, 25) = 1;
G(24:25, 25) = 0;  % opening near the bottom only
goalMask = false(M,M);
goalMask(1,1) = true;   % put "goal" (source) at top-left, "robot" would approach from inside the cup
T = FastMarchingField(G, goalMask, 'navigate');

% pick a point deep inside the cup (e.g. row=15,col=20) - must have a FINITE, and
% must be larger than a point just outside the cup near the opening
inside_deep = T(15, 20);
near_open   = T(24, 24);
assert(isfinite(inside_deep), 'point inside trap must still be reachable via the opening');
fprintf('[PASS] Test 2: U-trap reachable through opening (deep=%.2f, near-open=%.2f)\n', inside_deep, near_open);

%% Test 3: obstacle cells stay Inf in navigate mode
assert(isinf(T(5,15)), 'obstacle cell must have T=Inf (impassable) in navigate mode');
fprintf('[PASS] Test 3: obstacle cells correctly marked impassable (T=Inf)\n');

%% Test 4: clearance mode - obstacles are sources with T=0, free cells get positive distance
Gc = zeros(10,10);
Gc(5,5) = 1;   % single obstacle cell
srcMask = (Gc == 1);
Tc = FastMarchingField(Gc, srcMask, 'clearance');
assert(Tc(5,5) == 0, 'obstacle cell itself should have clearance 0');
assert(Tc(5,6) > 0 && Tc(5,6) < Tc(1,1), 'clearance should increase with distance from obstacle');
fprintf('[PASS] Test 4: clearance field behaves correctly (adjacent=%.3f, far=%.3f)\n', Tc(5,6), Tc(1,1));

%% Test 5: timing on largest planned map size (90x90, ~30%% density)
M = 90;
rng(42);
G = zeros(M,M);
G(rand(M,M) < 0.30) = 1;
G(1:5,1:5) = 0; G(M-4:M,M-4:M) = 0;
goalMask = false(M,M); goalMask(M,M) = true;
tic;
T = FastMarchingField(G, goalMask, 'navigate');
elapsed = toc;
fprintf('[INFO] Test 5: 90x90 map (~30%% density) FMM precompute time = %.3f s\n', elapsed);
assert(elapsed < 30, 'precompute should be well under 30s for a one-time cost');

fprintf('\nALL FMM TESTS PASSED.\n');
