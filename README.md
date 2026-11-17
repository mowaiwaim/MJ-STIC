# MJ-STIC: MATLAB code and data

Code and simulation data for **Mean-Square Passivity and Synchronization of Markovian Jump Multi-Weighted Coupled Neural Networks via Hybrid Control**, by Zi'an Wang and Yu Zhang.

## Run the examples

Use MATLAB R2024a. Download the repository and set its root as the MATLAB Current Folder. The simulations, plots and theorem-condition checks run in base MATLAB.

To redraw Figures 2–7 from the saved data, run:

```matlab
redraw_all
```

To rerun all three examples and the controller comparisons, run:

```matlab
reproduce_all
```

New plots go into `figures/`. Rerunning the simulations updates the CSV results. Figure 1 is the network diagram; its editable TikZ source is `figures/figure1.tex`.

## Find the code and results

| Paper item | Files |
|---|---|
| Table 1: coupling matrices | `examples/`, `check_example1_unified.m` |
| Example 1, Figure 2 | `examples/passivity/` |
| Example 2, Figures 3–4 | `examples/synchronization/` |
| Figure 5: network without control | `examples/uncontrolled/` |
| Saved data for Figures 2–5 | `early_figure_data.mat` |
| Table 2: triggering times | `outputs/table2_schedule.csv` |
| Example 3, Figure 6 | `run_bounded_wait.m`, `bw_*.csv` |
| Table 3: waiting times | `outputs/table3_schedule.csv` |
| Figure 7, Table 4 | `run_comparison.m`, `outputs/comparison_*.csv` |
| Shorter maximum waiting time | `check_growth_parameters.m`, `bw_growth_parameters.csv` |

`export_timing_tables.m` generates Tables 2 and 3 from the simulation records. In `outputs/`, `comparison_summary.csv` contains the Table 4 averages, `comparison_trajectories.csv` contains the Figure 7 curves, and `comparison_paths.csv` contains the results for each path. The `zhan_*`, `shang_*` and `liu_*` files contain the three recent comparison methods' results. Their formulas and implementation details are in the three `*_method.txt` files.

## Random paths

Examples 1 and 2 use `rng(2,'twister')` with a time step of 0.001.

Example 3 uses 256 independent pilot paths to choose controller parameters, followed by 384 shared test paths to compare the methods. The test interval is 0–8 and the time step is 0.005. The Brownian and Markov seeds are 512091 and 512092 for the pilot paths, and 73129 and 73130 for the test paths. `bw_paths.m` generates these paths. The bootstrap seed is 39147.

All methods use the same network and test paths. `IntegratedMSE` is the time integral of the mean-square synchronization error, called J in the paper. The CSV files keep full numerical precision.

To check the theorem conditions using the reported matrices, run:

```matlab
check_example1_unified
bw_certificate(bw_parameters())
```
