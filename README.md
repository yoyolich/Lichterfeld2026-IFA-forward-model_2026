# IFA Forward Modelling Toolbox

A Python toolbox for simulating foraminiferal δ¹⁸O variability from ocean reanalysis data and generating synthetic Individual Foraminifera Analysis (IFA) datasets.

Developed as part of research on past climate variability in the tropical Indian Ocean.

---

## Overview

This repository contains **two complementary forward-modelling tools** designed to investigate how different modes of climate variability impact IFA-derived δ¹⁸O distributions.

### Included Scripts

| `Season_ForwardModeled_IFA.py` | Modifies the **seasonal cycle amplitude** of temperature and salinity |
| `Interannual_variability_ForwardModeled_IFA.py` | Modifies **interannual variability** associated with ENSO/IOD-like climate events |

---

# 1. Seasonal Variability Forward Model

## Description

This script:

1. Extracts temperature and salinity time series at a selected site/depth from ORAS5 NetCDF files  
2. Computes seawater and foraminiferal calcite δ¹⁸O  
3. Amplifies or reduces the seasonal cycle by a user-defined coefficient  
4. Generates synthetic IFA populations by random sampling  
5. Quantifies sampling uncertainty via bootstrap resampling  
6. Exports modified and original time series to Excel  



# 2. Interannual Variability Forward Model

## Description

This script simulates the impact of **enhanced or reduced interannual variability** on IFA distributions by modifying temperature and salinity anomalies associated with major ENSO/IOD-type climate events.

It:

1. Extracts T/S time series from ORAS5 NetCDF files  
2. Modifies interannual variability during predefined ENSO/IOD event years  
3. Recalculates δ¹⁸Ow and δ¹⁸Oc  
4. Generates Monte Carlo synthetic IFA populations  
5. Computes distribution statistics and KS tests  
6. Exports all results to Excel  

---

## Variability Modification Modes

| Mode | Description |
|------|-------------|
| `augment` | Enhances ENSO/IOD-related anomalies |
| `reduce` | Dampens ENSO/IOD-related anomalies |

---

# Input Data Requirements

Both scripts require NetCDF temperature and salinity files from **ORAS5 ocean reanalysis** (or equivalent products with compatible variable names):

| Variable | Description | NetCDF Variable Name |
|---------|-------------|----------------------|
| Temperature | Potential temperature | `votemper` |
| Salinity | Practical salinity | `vosaline` |

ORAS5 data is available from the **Copernicus Marine Service**:  
https://marine.copernicus.eu

---

# Scientific Background

## Seawater δ¹⁸O

Surface example:


δ¹⁸Ow = 0.26 × S − 8


Thermocline example:


δ¹⁸Ow = 0.45 × S − 15.43


Foraminiferal Calcite δ¹⁸O

δ¹⁸Oc = (δ¹⁸Ow − 0.27) + (3.10 − T / 4.8)


These equations can be modified directly within the scripts.

---

# Citation

If you use these tools in your research, please cite:

> Lichterfeld, Y. et al., 2026.  
> *Assessing the potential of Individual Foraminifera Analyses (δ18O) to reconstruct variability, seasonality and extremes in the tropical Indian Ocean.*

---

# License

MIT License — see `LICENSE` file for details.


# Contact

If you have any questions or remarks, please feel free to contact me on my mail, i will be happy to help you.

**Yohan Lichterfeld**  
lichterfeld@cerege.fr
