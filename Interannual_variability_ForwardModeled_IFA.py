#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: yohanlichterfeld
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

======================================
Forward modelling of Individual Foraminifera Analysis (IFA) with 
interannual variability modification:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
from pathlib import Path
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import ks_2samp, gaussian_kde, skew, kurtosis
import warnings
warnings.filterwarnings("ignore")

try:
    from netCDF4 import Dataset
except ImportError:
    raise ImportError("Please install netCDF4: pip install netCDF4")

# =============================================================================
# ===                         USER PARAMETERS                               ===
# =============================================================================

LAT_SITE   = "Enter lat"        # Site latitude (degrees N)
LON_SITE   = "Enter lon"       # Site longitude (degrees E)
DEPTH_SITE = "Enter depht"        # Extraction depth (m)

TEMP_FILE  = "Your path file"    # Temperature NetCDF file
SAL_FILE   = "Your path file"    # Salinity NetCDF file

MODE       = "augment"     # 'augment': increases variability | 'reduce': reduces it


SPECIES    = "thermocline" # 'surface'     -> G. ruber     -> d18Ow = 0.26*S - 8  #You can put others equations
                           # 'thermocline' -> N. dutertrei -> d18Ow = 0.45*S - 15.43

N_IFA      = 60            # Number of IFA specimens (randomly picked)
N_SIMS     = 1000           # Number of Monte Carlo simulations

# === d18O Equations (can be modified directly here) ===
# Available variables: temp (C), sal (PSU), d18Ow (permil)
EQUATION_D18Ow_SURFACE     = "0.26 * sal - 8.0"          # G. ruber
EQUATION_D18Ow_THERMOCLINE = "0.45 * sal - 15.43"        # N. dutertrei
EQUATION_D18Oc             = "d18Ow + (16.5 - temp) / 4.8"

MEAS_ERROR_STD = 0.1       # Analytical noise (permil, 1-sigma)

OUTPUT_DIR  = "Your path file"
OUTPUT_XLSX = "Your path file.xlsx"

RANDOM_SEED = 42
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)

# =============================================================================
# ===              ENSO/IOD PARAMETERS — do not modify                      ===
# =============================================================================

ELNINO_YEARS    = [1958, 1966, 1973, 1983, 1987, 1992, 1998, 2002, 2010, 2015, 2016]
LANINA_YEARS    = [1971, 1974, 1976, 1989, 1999, 2000, 2008, 2011, 2012]
IOD_POS_YEARS   = [1961, 1963, 1967, 1972, 1982, 1994, 1997, 2006, 2007, 2012, 2015, 2017, 2018]
IOD_NEG_YEARS   = [1958, 1975, 1984, 1985, 1992, 1995, 1996, 1998, 2005, 2010, 2016]
ALL_EVENT_YEARS = sorted(set(ELNINO_YEARS + LANINA_YEARS + IOD_POS_YEARS + IOD_NEG_YEARS))

ELNINO_TEMP_DELTA     =  2.5
LANINA_TEMP_DELTA     = -2.0
SALINITY_AMP_FACTOR   =  0.75
ROLLING_WINDOW_MONTHS =  60
REDUCTION_PERCENT     =  75.0

# =============================================================================
# ===                          NETCDF FUNCTIONS                             ===
# =============================================================================

def load_and_extract(nc_file, var_name, lat_site, lon_site, depth_m):
    """
    Loads an ORAS5 NetCDF file, transposes to (time, [depth,] lat, lon) 
    and extracts the time series at (lat_site, lon_site, depth_m) 
    using bilinear interpolation.
    """
    ds = Dataset(nc_file, 'r')

    var_obj = ds.variables[var_name]
    dims    = list(var_obj.dimensions)

    def find_dim(keywords):
        for d in dims:
            if any(k in d.lower() for k in keywords):
                return d
        return None

    time_dim  = find_dim(['time'])
    depth_dim = find_dim(['depth'])
    lat_dim   = find_dim(['lat'])
    lon_dim   = find_dim(['lon'])

    if None in (time_dim, lat_dim, lon_dim):
        raise ValueError(f"Unidentified dimensions in {dims}")

    has_depth = depth_dim is not None
    depths    = np.array(ds.variables[depth_dim][:]) if has_depth else np.array([0.0])

    lats = np.array(ds.variables[lat_dim][:])
    lons = np.array(ds.variables[lon_dim][:])
    if lats.ndim == 2:
        lats = lats[:, 0]
        lons = lons[0, :]

    time_var = ds.variables[time_dim]
    try:
        from netCDF4 import num2date
        times_raw = num2date(time_var[:], units=time_var.units,
                             calendar=getattr(time_var, 'calendar', 'standard'))
        dates = pd.to_datetime([str(t) for t in times_raw])
    except Exception:
        dates = pd.date_range(start='1958-01', periods=len(time_var), freq='MS')

    data = np.array(var_obj[:])
    target = [time_dim, depth_dim, lat_dim, lon_dim] if has_depth else [time_dim, lat_dim, lon_dim]
    perm   = [dims.index(d) for d in target]
    data   = np.transpose(data, perm)
    ds.close()

    if hasattr(data, 'mask'):
        data = data.filled(np.nan)

    # Clean ORAS5 fill values (typically ~9.97e+36 or ~1e20)
    fill_threshold = 1e10
    data[np.abs(data) > fill_threshold] = np.nan

    n_times = data.shape[0]
    dates   = dates[:n_times]

    lat_sorted = np.sort(lats)
    lon_sorted = np.sort(lons)
    lat_flip   = lats[0] > lats[-1]

    if has_depth:
        depth_idx = int(np.argmin(np.abs(depths - depth_m)))
        print(f"  Requested depth: {depth_m} m  ->  ORAS5 level: {depths[depth_idx]:.1f} m")
    else:
        depth_idx = None

    series = np.full(n_times, np.nan)
    for t in range(n_times):
        sl = data[t, depth_idx, :, :] if depth_idx is not None else data[t, :, :]
        if lat_flip:
            sl = np.flipud(sl)
        try:
            interp    = RegularGridInterpolator(
                (lat_sorted, lon_sorted), sl,
                method='linear', bounds_error=False, fill_value=np.nan)
            series[t] = interp([[lat_site, lon_site]])[0]
        except Exception:
            series[t] = np.nan

    # Diagnostic: print extracted series statistics
    valid = series[~np.isnan(series)]
    print(f"  Valid values: {len(valid)}/{n_times}")
    if len(valid) == 0:
        print(f"  !!! ALL VALUES ARE NaN.")
        print(f"      -> Grid extent: lat=[{lats.min():.2f}, {lats.max():.2f}]  "
              f"lon=[{lons.min():.2f}, {lons.max():.2f}]")
        print(f"      -> Requested site: lat={lat_site}  lon={lon_site}")
        in_lat = lats.min() <= lat_site <= lats.max()
        in_lon = lons.min() <= lon_site <= lons.max()
        if not in_lat or not in_lon:
            print(f"      -> SITE OUTSIDE GRID! (lat_ok={in_lat}, lon_ok={in_lon})")
        else:
            print(f"      -> Site within grid but cell is masked (land or sea floor).")
        raise ValueError("Extraction impossible: no valid values. See diagnostic above.")
    print(f"  min={valid.min():.3f}  max={valid.max():.3f}  "
          f"std={valid.std():.4f}  mean={valid.mean():.3f}")
    if valid.std() < 1e-6:
        print("  !!! WARNING: std ~ 0 -> constant series (masked cell?).")

    return pd.Series(series.astype(float), index=dates, name=var_name)

# =============================================================================
# ===                 INTERANNUAL VARIABILITY MODIFICATION                  ===
# =============================================================================

def augment_temperature(temp, dates):
    mod   = temp.copy()
    years = pd.to_datetime(dates).year
    for y in ELNINO_YEARS:
        mod[(years == y) | (years == y + 1)] += ELNINO_TEMP_DELTA
    for y in LANINA_YEARS:
        mod[(years == y) | (years == y + 1)] += LANINA_TEMP_DELTA
    return mod

def augment_salinity(sal, dates):
    s    = pd.Series(sal, index=pd.to_datetime(dates))
    rm   = s.rolling(window=ROLLING_WINDOW_MONTHS, center=True, min_periods=1).mean()
    anom = s - rm
    mod  = s.copy()
    enso = set(ELNINO_YEARS) | set(LANINA_YEARS)
    mask = np.array([(y in enso) or ((y - 1) in enso) for y in s.index.year])
    mod.iloc[mask] += SALINITY_AMP_FACTOR * anom.iloc[mask]
    return mod.values

def reduce_temperature(temp, dates):
    mean  = np.nanmean(temp)
    mod   = temp.copy()
    years = pd.to_datetime(dates).year
    for i, y in enumerate(years):
        if y in ALL_EVENT_YEARS:
            mod[i] = mean + (temp[i] - mean) * (1 - REDUCTION_PERCENT / 100)
    return mod

def reduce_salinity(sal, dates):
    mod   = sal.copy()
    years = pd.to_datetime(dates).year
    for y in np.unique(years):
        idx = np.where(years == y)[0]
        if y in ALL_EVENT_YEARS:
            vals   = sal[idx]
            center = vals.mean()
            r_amp  = (vals.max() - vals.min()) * (1 - REDUCTION_PERCENT / 100)
            for i in idx:
                if sal[i] > center:
                    mod[i] = min(center + r_amp / 2, sal[i])
                elif sal[i] < center:
                    mod[i] = max(center - r_amp / 2, sal[i])
    return mod

# =============================================================================
# ===                         d18O AND IFA CALCULATION                      ===
# =============================================================================

def compute_d18Ow(sal, species):
    eq = EQUATION_D18Ow_SURFACE if species == 'surface' else EQUATION_D18Ow_THERMOCLINE
    return eval(eq, {"sal": sal, "np": np})

def compute_d18Oc(temp, d18Ow):
    return eval(EQUATION_D18Oc, {"temp": temp, "d18Ow": d18Ow, "np": np})

def forward_ifa(d18Oc_monthly, n_ifa, n_sims, noise_std):
    clean    = d18Oc_monthly[~np.isnan(d18Oc_monthly)]
    all_sims = np.full((n_sims, n_ifa), np.nan)
    for i in range(n_sims):
        s           = np.random.choice(clean, size=n_ifa, replace=True)
        all_sims[i] = s + np.random.normal(0, noise_std, n_ifa)
    pooled = all_sims.flatten()
    stats  = {
        'mean'    : np.nanmean(pooled),
        'std'     : np.nanstd(pooled),
        'skewness': skew(pooled[~np.isnan(pooled)]),
        'kurtosis': kurtosis(pooled[~np.isnan(pooled)]),
        'p5'      : np.nanpercentile(pooled,  5),
        'p25'     : np.nanpercentile(pooled, 25),
        'p50'     : np.nanpercentile(pooled, 50),
        'p75'     : np.nanpercentile(pooled, 75),
        'p95'     : np.nanpercentile(pooled, 95),
    }
    return all_sims, stats

# =============================================================================
# ===                                FIGURES                                ===
# =============================================================================

def plot_timeseries(dates, temp, sal, temp_mod, sal_mod, mode, site_label):
    fig, axes = plt.subplots(2, 1, figsize=(14, 6), sharex=True)

    axes[0].plot(dates, temp,     color='steelblue', lw=0.9, label='Original',     alpha=0.85)
    axes[0].plot(dates, temp_mod, color='firebrick', lw=0.9, label=f'Modified ({mode})', alpha=0.85)
    axes[0].set_ylabel("Temperature (C)", fontsize=10)
    axes[0].legend(fontsize=9)
    axes[0].grid(alpha=0.25)
    axes[0].set_title("ORAS5 Monthly Temperature", fontsize=10, loc='left')

    axes[1].plot(dates, sal,     color='steelblue',  lw=0.9, label='Original',         alpha=0.85)
    axes[1].plot(dates, sal_mod, color='darkorange',  lw=0.9, label=f'Modified ({mode})', alpha=0.85)
    axes[1].set_ylabel("Salinity (PSU)", fontsize=10)
    axes[1].legend(fontsize=9)
    axes[1].grid(alpha=0.25)
    axes[1].set_title("ORAS5 Monthly Salinity", fontsize=10, loc='left')

    plt.suptitle(f"{site_label}  |  MODE: {mode.upper()}", fontsize=11, fontweight='bold')
    plt.tight_layout()
    plt.show()


def plot_ifa_results(dates, d18Oc_orig, d18Oc_mod,
                     ifa_orig, ifa_mod, stats_orig, stats_mod,
                     mode, species, site_label):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # (a) Monthly d18Oc
    axes[0].plot(dates, d18Oc_orig, color='steelblue', lw=0.8, label='Original', alpha=0.85)
    axes[0].plot(dates, d18Oc_mod,  color='firebrick',  lw=0.8, label=f'Modified', alpha=0.85)
    axes[0].set_ylabel("d18Oc (permil VPDB)", fontsize=10)
    axes[0].set_xlabel("Date", fontsize=10)
    axes[0].legend(fontsize=9)
    axes[0].grid(alpha=0.25)
    axes[0].set_title("(a) Monthly d18Oc", fontsize=10, loc='left')
    axes[0].tick_params(axis='x', rotation=30)

    # (b) IFA Distributions
    po   = ifa_orig.flatten()
    pm   = ifa_mod.flatten()
    bins = np.linspace(min(po.min(), pm.min()) - 0.2,
                       max(po.max(), pm.max()) + 0.2, 40)
    axes[1].hist(po, bins=bins, density=True, alpha=0.4, color='steelblue',
                 label=f"Original  sigma={stats_orig['std']:.3f}")
    axes[1].hist(pm, bins=bins, density=True, alpha=0.4, color='firebrick',
                 label=f"Modified  sigma={stats_mod['std']:.3f}")
    for vals, col in [(po, 'steelblue'), (pm, 'firebrick')]:
        v    = vals[~np.isnan(vals)]
        kde  = gaussian_kde(v, bw_method='silverman')
        xkde = np.linspace(bins[0], bins[-1], 300)
        axes[1].plot(xkde, kde(xkde), color=col, lw=2)
    axes[1].set_xlabel("d18Oc (permil VPDB)", fontsize=10)
    axes[1].set_ylabel("Density", fontsize=10)
    axes[1].legend(fontsize=9)
    axes[1].grid(alpha=0.25)
    axes[1].set_title("(b) Forward modelled IFA distribution", fontsize=10, loc='left')

    # (c) Q-Q plot
    pct = np.arange(1, 100)
    axes[2].scatter(np.nanpercentile(po, pct), np.nanpercentile(pm, pct),
                    s=12, color='purple', alpha=0.6)
    lims = [min(po.min(), pm.min()), max(po.max(), pm.max())]
    axes[2].plot(lims, lims, 'k--', lw=1)
    axes[2].set_xlabel("Original Quantiles (permil)", fontsize=10)
    axes[2].set_ylabel("Modified Quantiles (permil)",  fontsize=10)
    axes[2].grid(alpha=0.25)
    axes[2].set_title("(c) Q-Q plot", fontsize=10, loc='left')

    sp_label = "G. ruber (surface)" if species == 'surface' else "N. dutertrei (thermocline)"
    plt.suptitle(f"{site_label}  |  {sp_label}  |  {mode.upper()}",
                 fontsize=11, fontweight='bold')
    plt.tight_layout()
    plt.show()

# =============================================================================
# ===                             MAIN PROGRAM                              ===
# =============================================================================

def main():

    site_label = f"lat={LAT_SITE}  lon={LON_SITE}  depth={DEPTH_SITE} m"
    print("=" * 60)
    print(f"  Site    : {site_label}")
    print(f"  Mode    : {MODE}   |   Species: {SPECIES}")
    print("=" * 60)

    # ------------------------------------------------------------------
    # 1. T/S Extraction
    # ------------------------------------------------------------------
    print("\n[1/4] Extracting T/S time series from NetCDF...")
    T_series = load_and_extract(TEMP_FILE, "votemper", LAT_SITE, LON_SITE, DEPTH_SITE)
    S_series = load_and_extract(SAL_FILE,  "vosaline", LAT_SITE, LON_SITE, DEPTH_SITE)

    df    = pd.DataFrame({'Temperature': T_series, 'vosaline': S_series}).dropna()
    df.index = pd.to_datetime(df.index)
    dates = df.index
    temp  = df['Temperature'].values
    sal   = df['vosaline'].values
    print(f"  -> {len(df)} months extracted ({dates[0].year}-{dates[-1].year})")

    # ------------------------------------------------------------------
    # 2. Interannual modification
    # ------------------------------------------------------------------
    print(f"\n[2/4] Modifying interannual variability (mode='{MODE}')...")
    if MODE == 'augment':
        temp_mod = augment_temperature(temp.copy(), dates)
        sal_mod  = augment_salinity(sal.copy(), dates)
    elif MODE == 'reduce':
        temp_mod = reduce_temperature(temp.copy(), dates)
        sal_mod  = reduce_salinity(sal.copy(), dates)
    else:
        raise ValueError(f"Invalid MODE: '{MODE}'. Use 'augment' or 'reduce'.")

    # Plot T/S time series
    plot_timeseries(dates, temp, sal, temp_mod, sal_mod, MODE, site_label)

    # ------------------------------------------------------------------
    # 3. d18O Calculation + IFA
    # ------------------------------------------------------------------
    print("\n[3/4] Calculating d18Oc and IFA forward modelling...")
    d18Ow_orig = compute_d18Ow(sal,      SPECIES)
    d18Ow_mod  = compute_d18Ow(sal_mod, SPECIES)
    d18Oc_orig = compute_d18Oc(temp,     d18Ow_orig)
    d18Oc_mod  = compute_d18Oc(temp_mod, d18Ow_mod)

    ifa_orig, stats_orig = forward_ifa(d18Oc_orig, N_IFA, N_SIMS, MEAS_ERROR_STD)
    ifa_mod,  stats_mod  = forward_ifa(d18Oc_mod,  N_IFA, N_SIMS, MEAS_ERROR_STD)

    ks_stat, ks_pval = ks_2samp(ifa_orig.flatten(), ifa_mod.flatten())
    print(f"  KS test orig vs modified: stat={ks_stat:.4f}  p={ks_pval:.4e}")

    # Plot d18Oc + IFA distributions
    plot_ifa_results(dates, d18Oc_orig, d18Oc_mod,
                     ifa_orig, ifa_mod, stats_orig, stats_mod,
                     MODE, SPECIES, site_label)

    # ------------------------------------------------------------------
    # 4. Excel Export
    # ------------------------------------------------------------------
    print("\n[4/4] Exporting to Excel...")
    out_dir  = Path(OUTPUT_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / OUTPUT_XLSX

    with pd.ExcelWriter(str(out_path), engine='openpyxl') as writer:

        pd.DataFrame({
            'Date'                : dates,
            'Original_Temperature': temp,
            'Modified_Temperature': temp_mod,
            'Original_Salinity'   : sal,
            'Modified_Salinity'   : sal_mod,
            'Original_d18Ow'      : d18Ow_orig,
            'Modified_d18Ow'      : d18Ow_mod,
            'Original_d18Oc'      : d18Oc_orig,
            'Modified_d18Oc'      : d18Oc_mod,
        }).to_excel(writer, sheet_name='TimeSeries', index=False)

        n_exp = min(N_SIMS, 100)
        pd.DataFrame(
            {**{f"Sim_{i+1:03d}_orig": ifa_orig[i] for i in range(n_exp)},
             **{f"Sim_{i+1:03d}_mod" : ifa_mod[i]  for i in range(n_exp)}}
        ).to_excel(writer, sheet_name='IFA_distributions', index=False)

        pd.DataFrame({
            'Statistic'  : list(stats_orig.keys()) + ['KS_stat', 'KS_pval'],
            'Original'   : list(stats_orig.values()) + [ks_stat, ks_pval],
            'Modified'   : list(stats_mod.values())  + [ks_stat, ks_pval],
        }).to_excel(writer, sheet_name='Stats', index=False)

    print(f"\mo  Workbook saved: {out_path}")
    print("=" * 60)
    print("Done.")


if __name__ == "__main__":
    main()
