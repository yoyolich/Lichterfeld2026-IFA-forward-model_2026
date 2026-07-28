# IFA Forward Modelling Toolbox

A Python toolbox for simulating foraminiferal δ¹⁸O variability from ocean reanalysis data and generating synthetic Individual Foraminifera Analysis (IFA) datasets.

Developed as part of research on past climate variability in the tropical Indian Ocean.


## Overview

This repository contains **three forward-modelling scripts** designed to investigate how different modes of climate variability impact IFA-derived δ¹⁸O distributions.

**`Forward-modelled_IFA_updated.py` is the recommended entry point** — it unifies and extends the two earlier standalone scripts (kept below for reference and reproducibility of earlier results).

### Included Scripts

| Script | Description |
|---|---|
| `Forward-modelled_IFA_updated.py` | **Recommended.** Combines seasonal and interannual testing, depth-selection options, equation uncertainty propagation, and FORCIS flux weighting in one configurable tool. |
| `Season_ForwardModeled_IFA.py` | Modifies the **seasonal cycle amplitude** of temperature and salinity |
| `Interannual_variability_ForwardModeled_IFA.py` | Modifies **interannual variability** associated with ENSO/IOD-like climate events |


# Master Forward Model (recommended)

## Description

`Forward-modelled_IFA_updated.py` supersedes the two standalone scripts below by combining seasonal and interannual variability testing, depth-selection options, equation uncertainty propagation, and FORCIS flux weighting into a single configurable tool.

### What it adds beyond the original two scripts

| Feature | `Season_...` / `Interannual_...` | `Forward-modelled_IFA_updated.py` |
|---|---|---|
| Seasonal cycle test | Yes (seasonal script only) | Yes |
| Interannual (ENSO/IOD) test | Yes (interannual script only) | Yes |
| Depth selection | Fixed depth only | Fixed **or** layer-averaged depth |
| Equation uncertainty propagation (Farmer 2007, δ¹⁸Ow-S) | No | Yes (Monte Carlo draw of a, b, slope, intercept per simulation) |
| FORCIS seasonal flux weighting | No | Yes (optional; falls back to uniform sampling) |
| Sites covered per run | Single site | All four core sites in one run |
| Statistical test | — | Wilcoxon signed-rank test (original vs. modified SD) |
| Output | Excel | PDF figures (barplot, time series, single IFA draw) + console summary |

### Configuration

All options are set at the top of the script:

```python
DEPTH_MODE = "layer"          # "fixed" or "layer"
USE_FLUX_WEIGHTING = True     # True = FORCIS, False = uniform sampling
ANALYSIS_TYPE = "seasonal"    # "seasonal" or "interannual"
MODE = "augment"              # "augment" or "reduce" (interannual only)
SEASON_COEF = 1.75            # seasonal amplitude multiplier
N_IFA = 60                    # individuals per synthetic IFA sample
N_SIMS = 500                  # Monte Carlo simulations
PROPAGATE_EQ_UNCERTAINTY = True
```

Set `TEMP_FILE`, `SAL_FILE`, `EXCEL_PATH` (if using flux weighting), and `OUTPUT_DIR` to your local paths, then run:

```bash
python Forward-modelled_IFA_updated.py
```

### What it does

1. Extracts ocean reanalysis temperature/salinity (e.g. ORAS5, or any compatible product) at each of the four sites (MD96-2060, NIOP 905, U1467, SO189-39KL), at fixed depth or averaged over a depth layer.
2. Applies the chosen climate modification (seasonal amplification/damping, or interannual ENSO/IOD augmentation/reduction).
3. Converts T/S to δ¹⁸Oc via the seawater δ¹⁸O-salinity relationship and the species-specific paleotemperature equation (Farmer et al. 2007), optionally drawing equation parameters from their reported uncertainty on every Monte Carlo simulation.
4. Draws synthetic IFA populations (N=60 individuals, uniform or FORCIS flux-weighted sampling) with added analytical noise, 500 times, for both the original and the modified climate.
5. Reports the mean ± SD of the IFA standard deviation for original vs. modified conditions at each site/species, with a Wilcoxon signed-rank test for significance.
6. Saves three sets of figures per configuration: a summary barplot across all sites/species, per-site time series (original vs. modified), and a single example IFA draw (N=60) per site/species for a quick distribution check.

### Notes

- The two standalone scripts below (`Season_ForwardModeled_IFA.py`, `Interannual_variability_ForwardModeled_IFA.py`) are kept in the repository for reference/reproducibility of earlier results, but `Forward-modelled_IFA_updated.py` is the version used for all current analyses and is the one to build on.
- FORCIS flux weighting requires a sediment-trap Excel file with one sheet per site (see `SITE_CONFIG_FLUX` in the script for the expected sheet/column names). If unavailable, set `USE_FLUX_WEIGHTING = False` to fall back to uniform monthly sampling.


# 1. Seasonal Variability Forward Model (legacy)

## Description

This script:

1. Extracts temperature and salinity time series at a selected site/depth from ocean reanalysis NetCDF files
2. Computes seawater and foraminiferal calcite δ¹⁸O
3. Amplifies or reduces the seasonal cycle by a user-defined coefficient
4. Generates synthetic IFA populations by random sampling
5. Quantifies sampling uncertainty via bootstrap resampling
6. Exports modified and original time series to Excel

# 2. Interannual Variability Forward Model (legacy)

## Description

This script simulates the impact of **enhanced or reduced interannual variability** on IFA distributions by modifying temperature and salinity anomalies associated with major ENSO/IOD-type climate events.

It:

1. Extracts T/S time series from ocean reanalysis NetCDF files
2. Modifies interannual variability during predefined ENSO/IOD event years
3. Recalculates δ¹⁸Ow and δ¹⁸Oc
4. Generates Monte Carlo synthetic IFA populations
5. Computes distribution statistics and KS tests
6. Exports all results to Excel


## Variability Modification Modes

| Mode | Description |
|------|-------------|
| `augment` | Enhances ENSO/IOD-related anomalies |
| `reduce` | Dampens ENSO/IOD-related anomalies |


# Input Data Requirements

All three scripts require NetCDF temperature and salinity files from an ocean reanalysis product (e.g. **ORAS5**, or any equivalent product with compatible variable names):

| Variable | Description | NetCDF Variable Name |
|---------|-------------|----------------------|
| Temperature | Potential temperature | `votemper` |
| Salinity | Practical salinity | `vosaline` |

ORAS5 data is available from the **Copernicus Marine Service**:
https://marine.copernicus.eu


# Scientific Background

## Seawater δ¹⁸O

Surface example:

δ¹⁸Ow = 0.26 × S − 8

Thermocline example:

δ¹⁸Ow = 0.45 × S − 15.43

Foraminiferal Calcite δ¹⁸O

δ¹⁸Oc = (δ¹⁸Ow − 0.27) + (3.10 − T / 4.8)

These equations can be modified directly within the scripts.


# Citation

If you use these tools in your research, please cite:
> Lichterfeld, Y. et al., 2026.
> *Assessing the potential of Individual Foraminifera Analyses (δ18O) to reconstruct variability, seasonality and extremes in the tropical Indian Ocean.*


# License

MIT License — see `LICENSE` file for details.

# Contact

If you have any questions or remarks, please feel free to contact me on my mail, i will be happy to help you.

**Yohan Lichterfeld**
lichterfeld@cerege.fr
