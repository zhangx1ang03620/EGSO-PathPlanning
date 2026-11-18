# Maps

## Benchmark maps (`benchmark_maps/`)

The five benchmark maps of the paper (Table 2), exported from `GenerateEGSOMaps.m`:

| File | Size | Obstacle density | Layout |
|---|---|---:|---|
| `Map1.csv` | 50×50 | 11.8% | single decoy-channel dead-end trap |
| `Map2.csv` | 50×50 | 20.6% | multi-stage (3×) decoy corridor |
| `Map3.csv` | 60×60 | 8.7% | procedural Voronoi room/corridor layout |
| `Map4.csv` | 80×80 | 15.6% | sparse large-scale long-detour network |
| `Map5.csv` | 90×90 | 23.1% | dense warehouse shelf grid |

Each CSV is an M×M occupancy grid: `1` = obstacle, `0` = free. Row 1, column 1 is the start S = (1,1) and row M, column M is the goal E = (M,M), in array indexing (the figures in the paper draw the row index upwards). `benchmark_maps.mat` holds the same grids (`BenchmarkMaps.Map1` … `Map5`) and their specifications (`MapInfo`).

```matlab
G = readmatrix('maps/benchmark_maps/Map1.csv');
```

## Generators

| Script | Maps |
|---|---|
| `GenerateEGSOMaps.m` | the five benchmark maps above (fixed seeds 31101–31105) |
| `GenerateDecoyStressMap.m` | the stress-test maps of Sec. 4.9 (decoy stages 1–5, 60×60); the seeds and feasibility retries are in `experiments/Main_LocalMinima_StressTest.m` |
| `GenerateScalabilityMap.m` | the warehouse-style maps of Sec. 4.10 (M = 90 to 220) |
| `egso3d/GenerateEGSO3DMap.m` | the 20×20×8 voxel scene of Sec. 4.11 |

All generators are deterministic for a given seed, so every map used in the paper can be regenerated exactly.
