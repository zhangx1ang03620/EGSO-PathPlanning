# Results

Raw outputs of the runs reported in the paper. Each folder is written by the script(s) listed; rerunning a script overwrites the corresponding files here (the equal-budget scripts write to `runs/equal_nfe/` instead).

| Folder | Paper | Written by | Key files |
|---|---|---|---|
| `main/` | §4.3, §4.6–4.8, §4.10, §4.13 | `Run_Comparison_SingleMap` + `Aggregate_Comparison_Results`, `Run_ClassicalBaselines`, `Run_APFAStar_Baseline`, significance, timing, sensitivity and iteration-budget scripts | `EGSO_Comparison_Results.mat` (all 23 algorithms × 5 maps × 20 runs), `Statistics_Summary.txt`, `Classical_vs_EGSO_Summary.txt`, `Significance_v2_*.txt` |
| `equal_nfe/` | §4.4 | `experiments/equal_nfe/*` | `EqNFE_map*_fe*_r*.mat` (40 files), `EqualNFE_Summary.txt`, `nfe_t100/` (evaluations per run, Appendix B) |
| `ablation/` | §4.5 | `Main_Ablation_Experiment` | `Ablation_Results.txt`, `Ablation_Data.mat` |
| `factorial/`, `factorial2/`, `factorial_pooled/` | §4.5 | `Main_FactorialAblation`, `Main_FactorialAblation2`, `Run_PooledSeedingFields` | seeding-field comparison: `factorial_pooled/PooledSeedingFields_Summary.txt` |
| `localmin/` | §4.2 | `Run_LocalMinimaCensus`, `Run_ChamferCensus` | `LocalMinimaCensus_Summary.txt`, `ChamferCensus_Summary.txt` |
| `sampling/` | §4.7 | `Run_SamplingBaselines`, `Run_SamplingSignificance` | `Sampling_vs_EGSO_Summary.txt`, `Sampling_Significance.txt` |
| `stress/` | §4.9 | `Main_LocalMinima_StressTest`, `Run_StressLevel5_Connectivity` | `StressTest_Summary.txt`, `StressLevel5_Connectivity.txt` |
| `scalability/` | §4.10 | `Run_Scalability` | `Scalability_Summary.txt` |
| `3d/` | §4.11 | `Run_3DProofOfConcept` | see below |
| `noise/` | §4.12 | `Main_NoiseRobustness_v2`, `Main_ClearanceTunable_NoiseTest` | `NoiseRobustness_v2_Summary.txt` |

## Notes on specific files

- **Which implementation produced which file.** Results under the fixed iteration budget were produced with an earlier `EGSO.m` that also contained a field-directed elite step. That step never changed a candidate, and the current code gives bit-identical results (verified on all five maps for the main algorithm and all variants). The files that depend on run time or on evaluation counts were regenerated with the current code: `equal_nfe/` (EGSO columns), `equal_nfe/nfe_t100/NFE_T100_EGSO_map*.mat`, `main/CleanTiming_revised.mat`, `main/EndToEnd_Runtime_Table_Clean_revised.txt` and `3d/`.
- **3D proof-of-concept (§4.11).** The paper reports 20 runs per method: the 10 runs in `3d/` (current code) pooled with the 10 runs in `3d/previous_implementation/` (earlier code, different random streams).
- **Friedman statistic.** `main/Significance_v2_Friedman.txt` gives χ² = 74.30 with average ranks for tied failures; the paper's 79.29 breaks these ties by listing order (`Run_FriedmanPaperConvention`). `main/Significance_Friedman.txt` (75.28) is an earlier variant that replaced failures by a penalty value and is kept for completeness.
- **Sensitivity analysis set 1.** `main/SensitivityAnalysis_Summary.txt` also contains a `kappa` row: that parameter belonged to the removed elite step and has no effect, so the paper reports only the elite fraction and pool factor rows.
- **Timing.** The `Time(s)` column of `EGSO_Comparison_Results.mat` comes from five parallel processes and includes contention; the clean single-process times are in `main/EndToEnd_Runtime_Table_Clean_revised.txt`.
- `main/NFE_Comparison.txt` uses `FobjCounter`, which misses the evaluations made inside `LocalSearch.m`; use `equal_nfe/nfe_t100/` for complete counts.
