#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Forward-modelled_IFA_updated.py
================================
Master forward-model script for IFA-delta18Oc, combining in a single tool
everything previously split across the standalone seasonal and interannual
scripts. It supports:

  1. Depth selection: fixed depth or layer-averaged depth
  2. Species-specific equations after Farmer et al. (2007), with full
     uncertainty propagation
  3. Seawater delta18O (delta18Ow) equations, with uncertainty propagation
  4. Optional FORCIS flux weighting (seasonal calcification flux instead of
     a uniform sampling assumption)
  5. Monte Carlo IFA resampling (N=60 individuals, 500 simulations by
     default) with propagated equation uncertainty and analytical noise

This script supersedes and unifies:
  - Season_ForwardModeled_IFA.py      (seasonal cycle amplitude test)
  - Interannual_variability_ForwardModeled_IFA.py (ENSO/IOD interannual test)

Both modes ("seasonal" / "interannual") are now available from one
configuration block (see ANALYSIS_TYPE below), applied consistently across
all four core sites and both species/depths, with the same statistical
framework (Monte Carlo IFA + Wilcoxon test) used for both.

Equation references:
  delta18Oc : Farmer et al. (2007), Table 3
              G. ruber white   (a=15.4 +/- 0.675, b=4.78 +/- 0.270)
              N. dutertrei      (a=14.6 +/- 1.100, b=5.09 +/- 0.185)
  delta18Ow : Singh et al. (2010) via Thirumalai & Clemens (2020) -- surface
              (slope=0.26 +/- 0.02, intercept=-8.0 +/- 0.5)
              LeGrande & Schmidt (2006) -- thermocline
              (slope=0.447 +/- 0.002, intercept=-15.43 +/- 0.054)

Usage
-----
1. Set the paths in the CONFIGURATION block below (ocean reanalysis temperature and
   salinity NetCDF files, FORCIS flux Excel file, output directory).
   Any ocean reanalysis or model product can be used, as long as the NetCDF
   variable names (see NETCDF VARIABLE NAMES below) and dimension names match,
   or are adapted in the extraction functions. In this study, the following
   input files were used (given here only as a concrete example):
     - TEMP_FILE  : ORAS5 reanalysis, potential temperature   (TI_ORAS5.nc)
     - SAL_FILE   : ORAS5 reanalysis, salinity                (TIs_ORAS5.nc)
     - EXCEL_PATH : FORCIS sediment-trap flux data            (FORCIS_trap_Claude.xlsx)
2. Choose DEPTH_MODE ("fixed" or "layer"), ANALYSIS_TYPE ("seasonal" or
   "interannual"), MODE ("augment" or "reduce", interannual only), and
   whether to use FORCIS flux weighting (USE_FLUX_WEIGHTING).
3. Run the script: `python ForwardModel_Master.py`
4. Outputs (PDF figures + console summary) are written to OUTPUT_DIR:
   - A barplot comparing IFA standard deviation (original vs. modified
     climate) for all four sites and both species/depths, with Wilcoxon
     significance test.
   - Time series figures (original vs. modified) per site.
   - A single random IFA draw (N=60) per site/species, for a quick sanity
     check of the sampling distribution.

Requirements: xarray, numpy, pandas, matplotlib, scipy.
"""

import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import wilcoxon
import warnings
warnings.filterwarnings("ignore")

# =============================================================================
# ===                          USER CONFIGURATION                          ===
# =============================================================================

TEMP_FILE  = "<path to ocean reanalysis NetCDF file (potential temperature)>"
SAL_FILE   = "<path to ocean reanalysis NetCDF file (salinity)>"
EXCEL_PATH = "<path to FORCIS sediment-trap flux Excel file (optional)>"
OUTPUT_DIR = "<path to output directory>"

SITES = {
    "MD96-2060":  {"lat": -8.30, "lon": 40.46},
    "NIOP 905":   {"lat": 10.76, "lon": 51.95},
    "U1467":      {"lat":  4.50, "lon": 73.17},
    "SO189-39KL": {"lat": -0.78, "lon": 99.00},
}

# ── Depth selection ──────────────────────────────────────────────────────────
DEPTH_MODE = "layer"   # "fixed" or "layer"

DEPTH_SURF_FIXED  = 0.5
DEPTH_THERM_FIXED = 97.0

DEPTH_SURF_MIN,  DEPTH_SURF_MAX  = 0,   40
DEPTH_THERM_MIN, DEPTH_THERM_MAX = 70, 120

# ── Flux weighting ───────────────────────────────────────────────────────────
USE_FLUX_WEIGHTING = True   # True = FORCIS seasonal flux, False = uniform

# ── Climate modification mode ────────────────────────────────────────────────
ANALYSIS_TYPE = "seasonal"   # "seasonal" or "interannual"
MODE          = "augment"    # "augment" or "reduce" (interannual only)
SEASON_COEF   = 1.75

# ── Monte Carlo parameters ───────────────────────────────────────────────────
N_IFA          = 60
N_SIMS         = 500
NOISE_STD      = 0.08    # analytical noise (per mil)
RANDOM_SEED    = 42
PROPAGATE_EQ_UNCERTAINTY = True  # propagate uncertainty on a, b, slope, intercept

# ── Interannual event years ──────────────────────────────────────────────────
ELNINO_YEARS  = [1958, 1966, 1973, 1983, 1987, 1992, 1998, 2002, 2010, 2015, 2016]
LANINA_YEARS  = [1971, 1974, 1976, 1989, 1999, 2000, 2008, 2011, 2012]
IOD_POS_YEARS = [1961, 1963, 1967, 1972, 1982, 1994, 1997, 2006, 2007, 2012, 2015, 2017, 2018]
IOD_NEG_YEARS = [1958, 1975, 1984, 1985, 1992, 1995, 1996, 1998, 2005, 2010, 2016]
ALL_EVENT_YEARS = sorted(set(ELNINO_YEARS + LANINA_YEARS + IOD_POS_YEARS + IOD_NEG_YEARS))

ELNINO_TEMP_DELTA   =  2.5
LANINA_TEMP_DELTA   = -2.0
SALINITY_AMP_FACTOR =  0.75
ROLLING_WINDOW_MONTHS = 60
REDUCTION_PERCENT   = 75.0

np.random.seed(RANDOM_SEED)

# =============================================================================
# ===                EQUATION PARAMETERS WITH UNCERTAINTIES                 ===
# =============================================================================

# Farmer et al. (2007) Table 3 -- 95% CI -> SE = CI_width / 4
EQ_PARAMS = {
    "surface": {        # G. ruber white
        "a_mean": 15.4,  "a_se": 0.675,   # (16.8-14.1)/4
        "b_mean":  4.78, "b_se": 0.270,   # (5.26-4.18)/4
    },
    "thermocline": {    # N. dutertrei
        "a_mean": 14.6,  "a_se": 1.100,   # (16.9-12.5)/4
        "b_mean":  5.09, "b_se": 0.185,   # (5.29-4.55)/4
    },
}

# delta18Ow-S with uncertainties
D18OW_PARAMS = {
    "surface": {        # Singh et al. (2010) via Thirumalai & Clemens (2020)
        "slope_mean":     0.26,  "slope_se":     0.02,
        "intercept_mean": -8.0,  "intercept_se":  0.5,
    },
    "thermocline": {    # LeGrande & Schmidt (2006)
        "slope_mean":     0.447, "slope_se":     0.002,
        "intercept_mean": -15.43,"intercept_se":  0.054,
    },
}

# =============================================================================
# ===                     OCEAN REANALYSIS EXTRACTION                       ===
# =============================================================================

def normalize_time(times):
    return pd.DatetimeIndex(times).to_period("M").to_timestamp()

def extract_fixed(lat, lon, depth_surf, depth_therm):
    """
    Extract T/S time series at a single fixed depth per species.
    NOTE: 'votemper'/'vosaline' and the 'deptht'/'time_counter' dimension
    names follow the ORAS5/NEMO convention used in this study. If using a
    different reanalysis or model product, update these variable/dimension
    names accordingly.
    """
    ds_T = xr.open_dataset(TEMP_FILE)
    ds_S = xr.open_dataset(SAL_FILE)
    results = {}
    for species, depth in [("surface", depth_surf), ("thermocline", depth_therm)]:
        T = ds_T['votemper'].sel(lat=lat, lon=lon, method='nearest') \
                            .sel(deptht=depth, method='nearest')
        S = ds_S['vosaline'].sel(lat=lat, lon=lon, method='nearest') \
                            .sel(deptht=depth, method='nearest')
        s_T = pd.Series(T.values.flatten(), index=normalize_time(T.time_counter.values))
        s_S = pd.Series(S.values.flatten(), index=normalize_time(S.time_counter.values))
        df  = pd.concat([s_T.rename('temp'), s_S.rename('sal')], axis=1).dropna()
        df['month'] = df.index.month
        results[species] = df
    return results

def extract_layer(lat, lon, surf_min, surf_max, therm_min, therm_max):
    """Extract T/S time series averaged over a depth layer per species."""
    ds_T = xr.open_dataset(TEMP_FILE)
    ds_S = xr.open_dataset(SAL_FILE)
    results = {}
    for species, dmin, dmax in [("surface", surf_min, surf_max),
                                  ("thermocline", therm_min, therm_max)]:
        T = ds_T['votemper'].sel(lat=lat, lon=lon, method='nearest') \
                            .sel(deptht=slice(dmin, dmax)).mean(dim='deptht')
        S = ds_S['vosaline'].sel(lat=lat, lon=lon, method='nearest') \
                            .sel(deptht=slice(dmin, dmax)).mean(dim='deptht')
        s_T = pd.Series(T.values.flatten(), index=normalize_time(T.time_counter.values))
        s_S = pd.Series(S.values.flatten(), index=normalize_time(S.time_counter.values))
        df  = pd.concat([s_T.rename('temp'), s_S.rename('sal')], axis=1).dropna()
        df['month'] = df.index.month
        results[species] = df
    return results

# =============================================================================
# ===                     CLIMATE MODIFICATIONS                            ===
# =============================================================================

def apply_seasonal(df):
    """Amplify or dampen the seasonal cycle by SEASON_COEF."""
    df = df.copy()
    for var in ['temp', 'sal']:
        year  = df.index.year
        month = df.index.month
        ann_mean = df.groupby(year)[var].transform('mean')
        mon_mean = df.groupby([year, month])[var].transform('mean')
        anomaly  = mon_mean - ann_mean
        df[var]  = df[var] + anomaly * (SEASON_COEF - 1)
    return df

def apply_interannual(df):
    """Augment or reduce interannual (ENSO/IOD-related) anomalies."""
    df    = df.copy()
    years = df.index.year
    temp  = df['temp'].values.copy()
    sal   = df['sal'].values.copy()

    if MODE == 'augment':
        for y in ELNINO_YEARS:
            mask = (years == y) | (years == y + 1)
            temp[mask] += ELNINO_TEMP_DELTA
        for y in LANINA_YEARS:
            mask = (years == y) | (years == y + 1)
            temp[mask] += LANINA_TEMP_DELTA
        s_sal = pd.Series(sal, index=df.index)
        rm    = s_sal.rolling(ROLLING_WINDOW_MONTHS, center=True, min_periods=1).mean()
        anom  = s_sal - rm
        enso  = set(ELNINO_YEARS) | set(LANINA_YEARS)
        mask_s = np.array([(y in enso) or ((y-1) in enso) for y in years])
        s_sal.iloc[mask_s] += SALINITY_AMP_FACTOR * anom.iloc[mask_s]
        sal = s_sal.values

    elif MODE == 'reduce':
        mean_T = np.nanmean(temp)
        for i, y in enumerate(years):
            if y in ALL_EVENT_YEARS:
                temp[i] = mean_T + (temp[i] - mean_T) * (1 - REDUCTION_PERCENT / 100)
        for y in np.unique(years):
            idx = np.where(years == y)[0]
            if y in ALL_EVENT_YEARS:
                vals   = sal[idx]
                center = vals.mean()
                r_amp  = (vals.max() - vals.min()) * (1 - REDUCTION_PERCENT / 100)
                for i in idx:
                    sal[i] = center + np.clip(sal[i] - center, -r_amp/2, r_amp/2)

    df['temp'] = temp
    df['sal']  = sal
    return df

def apply_modification(df):
    if ANALYSIS_TYPE == 'seasonal':
        return apply_seasonal(df)
    elif ANALYSIS_TYPE == 'interannual':
        return apply_interannual(df)
    return df

# =============================================================================
# ===                     FORCIS FLUX WEIGHTING                            ===
# =============================================================================

SITE_CONFIG_FLUX = {
    "MD96-2060":  {"sheet": "MOZ_MD96-2060",  "ruber_col": "g_ruber_any_LT",
                   "dut_col": "n_dutertrei_LT",  "date_col": "date"},
    "NIOP 905":   {"sheet": "AS8_NIOP905",     "ruber_col": "g_ruber_raw",
                   "dut_col": "g_dutertrei_raw", "date_col": "date"},
    "U1467":      {"sheet": "EAS_U1467",       "ruber_col": "g_ruber_raw",
                   "dut_col": "g_dutertrei_raw", "date_col": "date"},
    "SO189-39KL": {"sheet": "JAM_SO189-39KL",  "ruber_col": "g_ruber_any_LT",
                   "dut_col": "n_dutertrei_VT",  "date_col": "date"},
}

def compute_flux_weights(site_name):
    """Compute monthly sampling weights from FORCIS sediment-trap flux data."""
    cfg = SITE_CONFIG_FLUX[site_name]
    df  = pd.read_excel(EXCEL_PATH, sheet_name=cfg["sheet"])
    df[cfg["date_col"]] = pd.to_datetime(df[cfg["date_col"]], errors='coerce')
    df["month"] = df[cfg["date_col"]].dt.month
    df = df.dropna(subset=["month"])
    def to_weights(series):
        arr = np.array([series.get(m, np.nan) for m in range(1, 13)])
        arr = np.where(np.isnan(arr), np.nanmean(arr), arr)
        return arr / arr.sum() * 12
    ruber = df.groupby("month")[cfg["ruber_col"]].mean()
    dut   = df.groupby("month")[cfg["dut_col"]].mean()
    return {"surface": to_weights(ruber), "thermocline": to_weights(dut)}

# =============================================================================
# ===                     IFA FORWARD MODEL                                ===
# =============================================================================

def compute_d18Oc_single(temp, sal, species, rng=None):
    """
    Compute delta18Oc for a T/S vector.
    If PROPAGATE_EQ_UNCERTAINTY, randomly draw a, b, slope, intercept from
    their reported uncertainty.
    """
    p_eq  = EQ_PARAMS[species]
    p_d18 = D18OW_PARAMS[species]

    if PROPAGATE_EQ_UNCERTAINTY and rng is not None:
        a         = rng.normal(p_eq["a_mean"],         p_eq["a_se"])
        b         = rng.normal(p_eq["b_mean"],         p_eq["b_se"])
        slope_w   = rng.normal(p_d18["slope_mean"],    p_d18["slope_se"])
        interc_w  = rng.normal(p_d18["intercept_mean"],p_d18["intercept_se"])
    else:
        a        = p_eq["a_mean"]
        b        = p_eq["b_mean"]
        slope_w  = p_d18["slope_mean"]
        interc_w = p_d18["intercept_mean"]

    d18Ow = slope_w * np.array(sal) + interc_w
    d18Oc = d18Ow + (a - np.array(temp)) / b
    return d18Oc


def forward_ifa(df, species, site_name, flux_weights=None):
    """
    Monte Carlo IFA resampling.
    Returns (stds, means, d18Oc_series):
      - stds, means: arrays of shape (N_SIMS,)
      - d18Oc_series: raw delta18Oc time series (mean-parameter equations,
        used for plotting only)
    """
    rng   = np.random.default_rng(RANDOM_SEED)
    stds, means = np.zeros(N_SIMS), np.zeros(N_SIMS)

    # Time series computed with mean equation parameters (no propagation),
    # used only for the time-series plot.
    d18Oc_series = compute_d18Oc_single(df['temp'].values, df['sal'].values,
                                         species, rng=None)

    # Flux weights
    if flux_weights is not None:
        fw = flux_weights[species]
        w  = np.array([fw[m - 1] for m in df['month'].values])
        w  = w / w.sum()
    else:
        w = None

    for i in range(N_SIMS):
        d18Oc = compute_d18Oc_single(df['temp'].values, df['sal'].values,
                                      species, rng=rng if PROPAGATE_EQ_UNCERTAINTY else None)
        clean = d18Oc[np.isfinite(d18Oc)]
        if len(clean) < N_IFA:
            stds[i] = np.nan; means[i] = np.nan; continue

        if w is not None:
            w_clean = w[np.isfinite(d18Oc)]
            w_clean = w_clean / w_clean.sum()
            sample  = rng.choice(clean, size=N_IFA, replace=True, p=w_clean)
        else:
            sample  = rng.choice(clean, size=N_IFA, replace=True)

        sample  += rng.normal(0, NOISE_STD, N_IFA)
        stds[i]  = np.std(sample, ddof=1)
        means[i] = np.mean(sample)

    return stds, means, d18Oc_series

# =============================================================================
# ===                              MAIN                                    ===
# =============================================================================

print("="*70)
print("  FORWARD MODEL MASTER -- IFA delta18Oc")
print(f"  Depth mode        : {DEPTH_MODE}")
print(f"  Analysis type     : {ANALYSIS_TYPE} ({MODE})")
print(f"  Flux weighting    : {'FORCIS' if USE_FLUX_WEIGHTING else 'Uniform'}")
print(f"  Eq. uncertainty   : {PROPAGATE_EQ_UNCERTAINTY}")
print(f"  N_IFA={N_IFA}, N_SIMS={N_SIMS}, noise={NOISE_STD} per mil")
print("="*70)

# Load flux weights if needed
flux_weights_all = {}
if USE_FLUX_WEIGHTING:
    print("\nLoading FORCIS flux weights...")
    for site in SITES:
        try:
            flux_weights_all[site] = compute_flux_weights(site)
            print(f"  {site} OK")
        except Exception as e:
            print(f"  {site} ERROR: {e} -- falling back to uniform flux")
            flux_weights_all[site] = None

results = {}

for site_name, coords in SITES.items():
    lat, lon = coords['lat'], coords['lon']
    print(f"\n{'-'*60}")
    print(f"  {site_name}")
    print(f"{'-'*60}")
    results[site_name] = {}

    # Extraction
    if DEPTH_MODE == "fixed":
        data = extract_fixed(lat, lon, DEPTH_SURF_FIXED, DEPTH_THERM_FIXED)
    else:
        data = extract_layer(lat, lon, DEPTH_SURF_MIN, DEPTH_SURF_MAX,
                             DEPTH_THERM_MIN, DEPTH_THERM_MAX)

    fw = flux_weights_all.get(site_name) if USE_FLUX_WEIGHTING else None

    for species in ["surface", "thermocline"]:
        df_orig = data[species].copy()
        df_mod  = apply_modification(df_orig.copy())

        # Original climate
        stds_orig, means_orig, ts_orig = forward_ifa(df_orig, species, site_name, fw)
        # Modified climate
        stds_mod,  means_mod,  ts_mod  = forward_ifa(df_mod,  species, site_name, fw)

        results[site_name][species] = {
            "orig":     {"stds": stds_orig, "means": means_orig},
            "mod":      {"stds": stds_mod,  "means": means_mod},
            "ts_orig":  pd.Series(ts_orig,  index=df_orig.index),
            "ts_mod":   pd.Series(ts_mod,   index=df_mod.index),
            "df_orig":  df_orig,
        }

        sd_orig = np.nanmean(stds_orig)
        sd_mod  = np.nanmean(stds_mod)
        se_orig = np.nanstd(stds_orig)
        se_mod  = np.nanstd(stds_mod)
        delta   = sd_mod - sd_orig

        sp_label = "G. ruber" if species == "surface" else "N. dutertrei"
        print(f"\n  {sp_label}")
        print(f"    Original  : SD = {sd_orig:.4f} +/- {se_orig:.4f} per mil")
        print(f"    Modified  : SD = {sd_mod:.4f}  +/- {se_mod:.4f} per mil")
        print(f"    DeltaSD   : {delta:+.4f} per mil  ({delta/sd_orig*100:+.1f}%)")

        # Wilcoxon test: original vs. modified
        _, p_w = wilcoxon(stds_orig[np.isfinite(stds_orig)],
                          stds_mod[np.isfinite(stds_mod)])
        sig = "***" if p_w < 0.001 else "**" if p_w < 0.01 else "*" if p_w < 0.05 else "ns"
        print(f"    Wilcoxon p: {p_w:.4f} {sig}")

# =============================================================================
# ===                     SUMMARY FIGURE (BARPLOT)                         ===
# =============================================================================

fig, axes = plt.subplots(2, 4, figsize=(18, 9), sharey=False)
fig.suptitle(f"Forward Model Master -- IFA delta18Oc SD\n"
             f"Depth: {DEPTH_MODE} | Analysis: {ANALYSIS_TYPE} ({MODE}) | "
             f"Flux: {'FORCIS' if USE_FLUX_WEIGHTING else 'Uniform'} | "
             f"Eq. uncertainty: {PROPAGATE_EQ_UNCERTAINTY}",
             fontsize=10, fontweight='bold')

sp_labels = {"surface": "G. ruber", "thermocline": "N. dutertrei"}
site_list = list(SITES.keys())
col_orig, col_mod = '#4C72B0', '#DD8452'

for row, species in enumerate(["surface", "thermocline"]):
    for col, site_name in enumerate(site_list):
        ax = axes[row, col]
        res = results[site_name][species]

        sd_o  = np.nanmean(res["orig"]["stds"])
        sd_m  = np.nanmean(res["mod"]["stds"])
        se_o  = np.nanstd(res["orig"]["stds"])
        se_m  = np.nanstd(res["mod"]["stds"])

        ax.bar([0, 1], [sd_o, sd_m], yerr=[se_o, se_m],
               color=[col_orig, col_mod], alpha=0.82, capsize=5,
               edgecolor='black', linewidth=0.6, width=0.5)

        ax.text(0, sd_o + se_o + 0.005, f'{sd_o:.3f}', ha='center', fontsize=8)
        ax.text(1, sd_m + se_m + 0.005, f'{sd_m:.3f}', ha='center', fontsize=8)

        delta = sd_m - sd_o
        ax.text(0.5, max(sd_o, sd_m) + max(se_o, se_m) + 0.025,
                f'Delta={delta:+.3f} permil', ha='center', fontsize=7.5,
                bbox=dict(boxstyle='round,pad=0.2', facecolor='lightyellow',
                          edgecolor='gray', alpha=0.8))

        ax.set_title(f"{site_name}\n{sp_labels[species]}", fontsize=8, fontweight='bold')
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['Original', f'Modified\n(+/-75%)'], fontsize=7)
        ax.set_ylabel('IFA SD (permil)' if col == 0 else '', fontsize=8)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(True, alpha=0.2, axis='y', linestyle='--')

plt.tight_layout()
Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
out_fig = f"{OUTPUT_DIR}/ForwardModel_Master_{DEPTH_MODE}_{ANALYSIS_TYPE}_{MODE}.pdf"
plt.savefig(out_fig, dpi=300, bbox_inches='tight')
print(f"\nBarplot figure saved -> {out_fig}")
plt.show()

# =============================================================================
# ===              FIGURE 2: ORIGINAL vs. MODIFIED TIME SERIES             ===
# =============================================================================

print("\nGenerating time series figures...")

for site_name in SITES:
    fig, axes = plt.subplots(2, 1, figsize=(16, 8), sharex=True)
    fig.suptitle(f"delta18Oc time series -- {site_name}\n"
                 f"Original vs {ANALYSIS_TYPE} ({MODE}) +/-{int((SEASON_COEF-1)*100)}%",
                 fontsize=12, fontweight='bold')

    sp_labels_ts = {"surface": "G. ruber (surface)", "thermocline": "N. dutertrei (thermocline)"}

    for row, species in enumerate(["surface", "thermocline"]):
        ax  = axes[row]
        res = results[site_name][species]
        ts_o = res["ts_orig"]
        ts_m = res["ts_mod"]

        ax.plot(ts_o.index, ts_o.values, color='steelblue', lw=0.9, alpha=0.9,
                label=f'Original  sigma={ts_o.std():.3f} permil')
        ax.plot(ts_m.index, ts_m.values, color='firebrick', lw=0.9, alpha=0.75,
                label=f'Modified  sigma={ts_m.std():.3f} permil')
        ax.invert_yaxis()
        ax.set_ylabel('delta18Oc (permil VPDB)', fontsize=10)
        ax.set_title(sp_labels_ts[species], fontsize=10, fontweight='bold')
        ax.legend(fontsize=9, loc='best', framealpha=0.9)
        ax.grid(True, alpha=0.2, linestyle='--')

    axes[-1].set_xlabel('Year', fontsize=10)
    plt.tight_layout()
    out_ts = f"{OUTPUT_DIR}/Timeseries_{site_name.replace(' ','_')}_{ANALYSIS_TYPE}_{MODE}.pdf"
    plt.savefig(out_ts, dpi=300, bbox_inches='tight')
    print(f"  {site_name} -> {out_ts}")
    plt.show()

# =============================================================================
# ===           FIGURE 3 + PRINT: SINGLE DRAW OF 60 IFA VALUES             ===
# =============================================================================

print("\n" + "="*70)
print("  SINGLE DRAW OF 60 IFA VALUES -- console summary")
print("="*70)

rng_single = np.random.default_rng(99)  # separate seed for this single draw

fig, axes = plt.subplots(2, 4, figsize=(18, 8))
fig.suptitle(f"Single IFA draw (N=60) -- Original conditions\n"
             f"Depth: {DEPTH_MODE} | Flux: {'FORCIS' if USE_FLUX_WEIGHTING else 'Uniform'}",
             fontsize=11, fontweight='bold')

print(f"\n{'Site':15s} {'Species':14s} {'Mean (permil)':>13s} {'SD (permil)':>13s} "
      f"{'Min (permil)':>13s} {'Max (permil)':>13s} {'N':>5s}")
print("-"*80)

for col, site_name in enumerate(SITES):
    for row, species in enumerate(["surface", "thermocline"]):
        ax  = axes[row, col]
        res = results[site_name][species]
        df_o = res["df_orig"]

        # Compute d18Oc series
        d18Oc = compute_d18Oc_single(df_o['temp'].values, df_o['sal'].values,
                                      species, rng=None)
        clean = d18Oc[np.isfinite(d18Oc)]

        # Flux weights
        fw = flux_weights_all.get(site_name) if USE_FLUX_WEIGHTING else None
        if fw is not None:
            ww = np.array([fw[species][m-1] for m in df_o['month'].values])
            ww_clean = ww[np.isfinite(d18Oc)]
            ww_clean = ww_clean / ww_clean.sum()
            sample = rng_single.choice(clean, size=N_IFA, replace=True, p=ww_clean)
        else:
            sample = rng_single.choice(clean, size=N_IFA, replace=True)
        sample += rng_single.normal(0, NOISE_STD, N_IFA)

        mean_ifa = np.mean(sample)
        sd_ifa   = np.std(sample, ddof=1)

        sp_short = "G. ruber" if species == "surface" else "N. dutertrei"
        print(f"{site_name:15s} {sp_short:14s} {mean_ifa:>13.4f} {sd_ifa:>13.4f} "
              f"{sample.min():>13.4f} {sample.max():>13.4f} {N_IFA:>5d}")

        # Scatter plot of the 60 IFA values
        color = '#4C72B0' if species == "surface" else '#DD8452'
        ax.scatter(range(N_IFA), np.sort(sample), s=20, color=color,
                   alpha=0.7, edgecolor='white', linewidth=0.3)
        ax.axhline(mean_ifa, color='black', lw=1.5, linestyle='--',
                   label=f'Mean={mean_ifa:.3f}')
        ax.axhline(mean_ifa + sd_ifa, color='gray', lw=1, linestyle=':',
                   label=f'SD={sd_ifa:.3f}')
        ax.axhline(mean_ifa - sd_ifa, color='gray', lw=1, linestyle=':')
        ax.invert_yaxis()
        ax.set_title(f"{site_name}\n{sp_short}", fontsize=8, fontweight='bold')
        ax.set_xlabel('Individual #', fontsize=7)
        ax.set_ylabel('delta18Oc (permil)' if col == 0 else '', fontsize=8)
        ax.legend(fontsize=7, loc='best', framealpha=0.8)
        ax.grid(True, alpha=0.2, linestyle='--')

print("-"*80)

plt.tight_layout()
out_ifa = f"{OUTPUT_DIR}/IFA_single_draw_{DEPTH_MODE}_{ANALYSIS_TYPE}.pdf"
plt.savefig(out_ifa, dpi=300, bbox_inches='tight')
print(f"\nIFA draw figure saved -> {out_ifa}")
plt.show()
print("\nDone.")
