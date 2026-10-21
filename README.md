[README.md](https://github.com/user-attachments/files/33788092/README.md)
# epidemic-compression-distortion

Code for **"Choosing the spatial scale of aggregation in epidemic models: compression versus distortion"**
(Y. R. Bhora, H. K. Das, I. Z. Kiss, K. V. Parag, E. Moro, L. M. Stolerman, V. Colizza, N. Kutz, M. Santillana).

This repository currently contains the county-network SIR simulations (Fig. 2 of the main text and the
network-density sensitivity analysis in the Supplementary Materials).

## What the simulation does

For one state (Texas, Colorado, or Massachusetts), `Fig_2_S2f_epidemic_compression_distortion_simulation.m`:

1. loads daily 2020 county-to-county flux matrices, φ_ij(t) = F_ij(t) / N_i;
2. thresholds every day at φ^th = 0.02, resets the diagonal so every row sums to one, and selects the
   least connected day of 2020 as a static network (Texas and Colorado: 12 April 2020; Massachusetts:
   17 December 2020);
3. simulates the metapopulation SIR-network model on that network (β₀ = 0.20, γ = 0.10 day⁻¹, one
   infectious individual in the most populous county) with `ode45` (relative tolerance 1e-7, absolute 1e-9);
4. aggregates prevalence to the county, MSA, and state levels;
5. repeats the simulation retaining the strongest 5%, 10%, 15%, and 100% of the links of the same day, and
   for a range of fixed thresholds, reporting the mean absolute optimal lag between each aggregate and its
   counties.

## Requirements

- MATLAB R2018b or later (R2020a or later for vector PDF export); tested with R2024b.
- No toolboxes required.

## Repository layout

```
Fig_2_S2f_epidemic_compression_distortion_simulation.m   main script
TDFluxMatrix.m                                            loads the daily flux matrices
SirGenNetworkModel.m                                      SIR-network model equations
InitialValueMobilitySIRNetwork.m                          initial conditions
data/
  Texas/          county_name_list.m, pop_vector_origin_only.csv, <County>_Mobility_Flux.csv
  Colorado/       (same)
  Massachusetts/  (same)
```

## How to run

1. Add the repository root to the MATLAB path.
2. Change the MATLAB working folder to one state's data folder, e.g. `data/Colorado`.
3. In Section 0 of the script, set `STATE = 'CO'` (`'TX'`, `'CO'`, or `'MA'`) and run it.
4. Figures and CSV tables are written to `result_figures/` inside that folder.

## Outputs

| File | Content |
|---|---|
| `Fig_MSA_state_county_minmax_*` | Fig. 2 panel: county, MSA, and state curves |
| `density_sweep_XX_rowstoch.*` | Density-sensitivity figure and table (Supplementary Materials) |
| `daily_mobility_snapshot_XX_rowstoch.*` | Daily inter-county mobility in 2020 with the selected day |
| `snapshot_date_XX_rowstoch.csv` | Selected day and share of links retained |
| `distortion_vs_threshold_XX_rowstoch.*` | Distortion as a function of the threshold |
| `lcc_vs_threshold_XX_rowstoch.*` | Largest connected component vs. threshold |

## Data

Mobility flows are derived from the GeoDS COVID19USFlows dataset (Kang et al., *Scientific Data* 7, 390,
2020; https://github.com/GeoDS/COVID19USFlows). Please cite that dataset and keep its license notice when
reusing the derived files.

## Citation

If you use this code, please cite the paper (citation to be added upon publication).

## License

See `LICENSE`.
