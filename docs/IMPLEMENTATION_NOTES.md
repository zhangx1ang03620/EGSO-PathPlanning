# Implementation notes

## Problem encoding (shared by all 23 algorithms)

- Grid `G` is M×M, `1` = obstacle, `0` = free. Start (1,1), goal (M,M).
- A candidate `x` holds one row index per interior column, `dim = M - 2`; the route is `[1, x, M]`. A path is therefore a single-valued function of the column index, and the benchmark maps are built from column-wise wall/gap structures.
- `fitness.m` rounds `x`, counts waypoints on obstacles (`1e6` per collision), and otherwise returns the length of the path obtained by `generateContinuousRoute` (cell-by-cell connectivity repair) followed by two passes of `GenerateSmoothPath` (line-of-sight shortcutting).
- A run is valid if its best fitness is below the collision threshold `1e4`. V.Mean, V.Best, V.Worst and V.Std are computed over valid runs; SR is the share of valid runs.

## EGSO modules (`egso/EGSO.m`)

| Module | Code | Default parameters |
|---|---|---|
| EVF | `FastMarchingField(G, goalMask, 'navigate')` and `(G, obstacleMask, 'clearance')`, first-order upwind fast marching, 4-neighbour | computed once per map |
| GSI | pool of `10N` candidates in total; in every column a free row is drawn with weight `exp(-eta (T - min T)) * exp(-0.15 |r - r_prev|)`; the best `N` are kept | `eta ~ U(0.02, 0.10)` per candidate |
| Elite layer | each of the best `Ne = max(3, floor(0.2N))` individuals is perturbed with probability `p`, `x + 0.5 a randn`, and replaced if strictly better | `p = 0.3`, `a = 2(1 - t/T)` |
| Scout layer | with probability `p_explore`: Lévy flight towards the best (Mantegna, `beta = 1.5`, scale `0.05a`) or differential mutation; otherwise Gaussian search `best + 0.3 a randn`; greedy acceptance | `p_explore = (1 - t/T)^3` if M > 50, else `0.8 (1 - t/T)` |
| Anti-stagnation | after `stall_limit` iterations without improvement and only for `0.2T < t < 0.7T`, the worst 70% of the population is re-drawn uniformly from free cells | `stall_limit = 5` |
| CALS | applied to the global best in every iteration after `0.65T`: three Laplacian passes, each interior waypoint is moved to the rounded midpoint of its neighbours and the move is kept only if the cell is free, its clearance `T_clear >= c_min`, and the fitness strictly improves; the result replaces the last individual | `c_min = 1.0` |

Every candidate is projected onto the encoding (rounded, clipped to [1, M], moved to the nearest free cell of its column), so decision variables stay integer row indices throughout.

## Removed elite step

Earlier versions of `EGSO.m` and its variants also moved each elite individual by at most 0.1 cells along the normalised field gradient. Because the projection rounds every candidate back to an integer row, that move never changed a candidate, and removing it leaves every fixed-budget result bit-identical (checked on the five maps for `EGSO` and all variants). Some variant signatures still accept the corresponding arguments (`kappa` in `EGSO_Sensitivity`, `elite_field_mode` in `EGSO_Variant2`) so that the experiment scripts run unchanged; these arguments have no effect.

## Baselines

- `baselines/metaheuristics/`: the 22 competing algorithms (SSA is the sparrow search algorithm of Xue and Shen, 2020), each with the interface `[best, x, curve] = ALGO(N, T, lb, ub, dim, fobj, G)` and the hyperparameters listed in Appendix A of the paper. None was retuned for these maps.
- GWO, PSO, GA, DE, SSA, INFO, SBOA and ARO call `LocalSearch.m` once per iteration: each coordinate of the current best is replaced by a random free row and kept if the fitness improves (`dim + 1` evaluations). This is not part of the canonical algorithms and favours these eight baselines.
- `baselines/planners/`: `AStarPath` (8-connected, octile heuristic), `APF_AStarPath` (A* with repulsive-potential edge cost), `FMMDescentPath` (greedy 4-neighbour descent on `T_nav`), `RRTStarPath` (RRT* and Informed-RRT* in continuous coordinates). All their paths go through the same `GenerateSmoothPath` post-processing as the swarm algorithms.

## Counting objective evaluations

- `egso/FobjCounter.m` wraps the function handle `fobj` and therefore misses the direct `fitness(x, G)` calls made inside `LocalSearch.m`.
- `experiments/equal_nfe/fitness.m` shadows `egso/fitness.m` when MATLAB's current folder is `experiments/equal_nfe` (the current folder takes precedence over the path). It counts every call, stops a run with error `NFEB:exhausted` at the budget, and records the best-so-far trace. The scripts in that folder `cd` there and assert that the counting wrapper is active.
- In the equal-budget protocol each algorithm's iteration count is set from its measured evaluations per run at T = 100 (`experiments/equal_nfe/calibration/`, with a 5% margin), so that its T-dependent schedules run to completion. Runs that finish before 98% of the budget are repeated with the same seed and a larger T (`Fixup_EarlyFinish`). The reported value is the best fitness evaluated within the budget.

## Statistics

- `stats/ManualRankSum.m`: two-sided Wilcoxon rank-sum test with tie and continuity correction, normal approximation.
- `stats/ManualFriedman.m`: Friedman test with average ranks for ties; the p-value uses `gammainc`.
- Neither needs the Statistics and Machine Learning Toolbox.

## Known limitations

- The one-row-per-column encoding cannot represent paths that revisit a column; spiral or maze layouts that require it are outside the benchmark design.
- The 3D proof-of-concept (`egso3d/`) uses a simpler encoding (straight segments between per-column waypoints, checked by `CollisionFreeSegment3D`) instead of a port of the 2D connectivity repair. It tests transferability, not a full 3D benchmark.
- `FastMarchingField` selects the next trial cell by a full-matrix argmin, i.e. O((M·M)²) in the worst case; the measured growth up to M = 220 is far below this (see `results/scalability/`).
