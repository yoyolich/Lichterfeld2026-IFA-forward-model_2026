#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: yohanlichterfeld
"""

import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import shapiro

plt.rcParams.update({'font.size': 12})

# =========================
# ===  USER PARAMETERS  ===
# =========================
LAT_SITE = "enter lat"       # Site latitude
LON_SITE = "enter lon"       # Site longitude
DEPTH_SITE = 0        # Depth in meters
TEMP_FILE = "Your path file"  #ORAS5 in this study
SAL_FILE = "Your path file"   #ORAS5 in this study
SEASON_COEF = 1.75     # Seasonality amplification factor (here 75%)
N_IFA = 60             # Number of individual foraminifera (IFA) to draw randomly

# === δ18O EQUATIONS (modifiable directly here) === 
EQUATION_D18Ow = "0.26 * sal - 8"                       # Subsurface: δ18Ow = 0.45*S - 15.43
EQUATION_D18Oc = "(d18Ow - 0.27) + (3.10 - temp / 4.8)"    

# =========================


def main():
    df = extract_point_data(TEMP_FILE, SAL_FILE, LAT_SITE, LON_SITE, DEPTH_SITE)

    df['δ18Ow'] = compute_d18Ow(df['vosaline'])
    df['δ18Oc estimated'] = compute_d18Oc(df['votemper'], df['δ18Ow'])

    df = alter_seasonality(df, 'votemper', coef=SEASON_COEF)
    df = alter_seasonality(df, 'vosaline', coef=SEASON_COEF)

    df['δ18Ow altered'] = compute_d18Ow(df['Altered vosaline'])
    df['Altered δ18Oc estimated'] = compute_d18Oc(df['Altered votemper'], df['δ18Ow altered'])

    plot_results(df)
    generate_IFA(df, N_IFA)
    ifa_bootstrap_std(df, n=N_IFA, n_iter=1000)

    export_to_excel_multi(df, filename="ORAS5_TempSal_δ18Oc_Series.xlsx")


def convert_lon_to_dataset(lon, lon_data):
    lon_min = lon_data.min().item()
    lon_max = lon_data.max().item()

    if (lon_min >= 0) and (lon_max <= 360):
        if lon < 0:
            lon += 360
    elif (lon_min >= -180) and (lon_max <= 180):
        if lon > 180:
            lon -= 360
    else:
        print("Warning: unexpected longitude format in dataset")

    return lon


def extract_point_data(nc_temp_file, nc_sal_file, lat_site, lon_site, depth_site):
    ds_temp = xr.open_dataset(nc_temp_file)
    ds_sal = xr.open_dataset(nc_sal_file)

    lon_site_temp = convert_lon_to_dataset(lon_site, ds_temp['lon'])
    lon_site_sal = convert_lon_to_dataset(lon_site, ds_sal['lon'])

    temp_point = ds_temp['votemper'].sel(
        lat=lat_site, lon=lon_site_temp, deptht=depth_site, method='nearest')
    sal_point = ds_sal['vosaline'].sel(
        lat=lat_site, lon=lon_site_sal, deptht=depth_site, method='nearest')

    df_temp = temp_point.to_dataframe().reset_index()
    df_sal = sal_point.to_dataframe().reset_index()

    df = pd.merge(df_temp, df_sal, on='time_counter', suffixes=('_temp', '_sal'))
    return df


def alter_seasonality(df, var_name, coef=1.75):
    df['year'] = df['time_counter'].dt.year
    df['month'] = df['time_counter'].dt.month

    annual_mean = df.groupby('year')[var_name].transform('mean')
    monthly_mean = df.groupby(['year', 'month'])[var_name].transform('mean')

    anomaly = monthly_mean - annual_mean
    df[f'Altered {var_name}'] = df[var_name] + anomaly * (coef - 1)

    df.drop(columns=['year', 'month'], inplace=True)
    return df


def compute_d18Ow(sal):
    return eval(EQUATION_D18Ow)


def compute_d18Oc(temp, d18Ow):
    return eval(EQUATION_D18Oc)


def plot_results(df):
    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)

    for ax in axes:
        ax.grid(True, linestyle='--', linewidth=0.5, alpha=1)
        ax.tick_params(axis='x', direction='in')
        ax.tick_params(axis='y', direction='in')

    axes[0].plot(df['time_counter'], df['votemper'], label='Original Temp.', color='blue')
    axes[0].plot(df['time_counter'], df['Altered votemper'], label='Altered Temp.', color='red', alpha=0.6)
    axes[0].set_ylabel('Temperature (°C)')
    axes[0].set_title(f'Temperature at lat={LAT_SITE}, lon={LON_SITE}, depth={DEPTH_SITE} m')
    axes[0].legend()

    axes[1].plot(df['time_counter'], df['vosaline'], label='Original Salinity', color='blue')
    axes[1].plot(df['time_counter'], df['Altered vosaline'], label='Altered Salinity', color='red', alpha=0.6)
    axes[1].set_ylabel('Salinity (PSU)')
    axes[1].legend()

    axes[2].plot(df['time_counter'], df['δ18Oc estimated'], label='Original δ18Oc', color='blue')
    axes[2].plot(df['time_counter'], df['Altered δ18Oc estimated'], label='Altered δ18Oc', color='red', alpha=0.6)
    axes[2].invert_yaxis()
    axes[2].set_ylabel('δ18Oc (‰)')
    axes[2].legend()
    axes[2].set_xlabel('Time')

    plt.tight_layout()
    plt.show()


def generate_IFA(df, n=60):
    original_vals = df['δ18Oc estimated'].dropna().values
    altered_vals = df['Altered δ18Oc estimated'].dropna().values

    # Random draw of IFAs
    ifa_original = np.random.choice(original_vals, size=n, replace=False)
    ifa_altered = np.random.choice(altered_vals, size=n, replace=False)

    # Display means
    mean_orig = np.mean(ifa_original)
    mean_alt = np.mean(ifa_altered)
    print(f"\nMean of {n} original IFAs: {mean_orig:.3f}")
    print(f"Mean of {n} altered IFAs: {mean_alt:.3f}")

    # Display drawn values
    print("\nValues of the 60 original IFAs:")
    print(np.round(ifa_original, 3))
    print("\nValues of the 60 altered IFAs:")
    print(np.round(ifa_altered, 3))

    # Standard deviations
    std_orig = np.std(ifa_original)
    std_alt = np.std(ifa_altered)
    print(f"\nStandard deviation of {n} original IFAs: {std_orig:.3f}")
    print(f"Standard deviation of {n} altered IFAs: {std_alt:.3f}")

    # Normality test
    _, p_orig = shapiro(ifa_original)
    _, p_alt = shapiro(ifa_altered)
    print(f"\nShapiro-Wilk normality test:")
    print(f"  Original : p = {p_orig:.4f} ({'normal' if p_orig > 0.05 else 'non-normal'})")
    print(f"  Altered  : p = {p_alt:.4f} ({'normal' if p_alt > 0.05 else 'non-normal'})")

    # Histograms
    plt.figure(figsize=(12, 5))
    plt.hist(ifa_original, bins=15, alpha=0.6, label='Original', color='blue')
    plt.hist(ifa_altered, bins=15, alpha=0.6, label='Altered', color='red')
    plt.xlabel('Estimated δ18Oc (IFA)')
    plt.ylabel('Frequency')
    plt.title(f'Histogram of randomly drawn IFAs (n={n})')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()


def ifa_bootstrap_std(df, n=60, n_iter=1000):
    """
    Perform n_iter random draws of n IFAs for original and altered δ18Oc.
    Display mean, std of standard deviations, and histograms.
    """
    original_vals = df['δ18Oc estimated'].dropna().values
    altered_vals = df['Altered δ18Oc estimated'].dropna().values

    stds_orig, stds_alt = [], []
    means_orig, means_alt = [], []

    for _ in range(n_iter):
        ifa_o = np.random.choice(original_vals, size=n, replace=False)
        ifa_a = np.random.choice(altered_vals, size=n, replace=False)

        stds_orig.append(np.std(ifa_o))
        stds_alt.append(np.std(ifa_a))
        means_orig.append(np.mean(ifa_o))
        means_alt.append(np.mean(ifa_a))

    mean_std_orig = np.mean(stds_orig)
    mean_std_alt = np.mean(stds_alt)
    std_std_orig = np.std(stds_orig)
    std_std_alt = np.std(stds_alt)

    mean_mean_orig = np.mean(means_orig)
    mean_mean_alt = np.mean(means_alt)
    std_mean_orig = np.std(means_orig)
    std_mean_alt = np.std(means_alt)

    print(f"\n📊 Results from {n_iter} draws of {n} IFAs:")
    print(f"  Mean std original : {mean_std_orig:.4f} ± {std_std_orig:.4f}")
    print(f"  Mean std altered  : {mean_std_alt:.4f} ± {std_std_alt:.4f}")
    print(f"\n  Mean of original means : {mean_mean_orig:.4f} ± {std_mean_orig:.4f}")
    print(f"  Mean of altered means  : {mean_mean_alt:.4f} ± {std_mean_alt:.4f}")

    # Histogram of standard deviations
    plt.figure(figsize=(12, 5))
    plt.hist(stds_orig, bins=30, alpha=0.6, label='Original std', color='blue')
    plt.hist(stds_alt, bins=30, alpha=0.6, label='Altered std', color='red')
    plt.xlabel('Standard deviation of IFAs')
    plt.ylabel('Frequency')
    plt.title(f'Distribution of std over {n_iter} draws (n={n})')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()

    # Histogram of means
    plt.figure(figsize=(12, 5))
    plt.hist(means_orig, bins=30, alpha=0.6, label='Original mean', color='blue')
    plt.hist(means_alt, bins=30, alpha=0.6, label='Altered mean', color='red')
    plt.xlabel('Mean δ18Oc (IFA)')
    plt.ylabel('Frequency')
    plt.title(f'Distribution of means over {n_iter} draws (n={n})')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()


def export_to_excel_multi(df, filename="Output_ORAS5.xlsx"):
    """
    Export original and altered series into a single Excel file with two sheets:
    - 'Original' : unaltered series
    - 'Altered'  : series with modified seasonality
    """
    cols_original = ['time_counter', 'votemper', 'vosaline', 'δ18Ow', 'δ18Oc estimated']
    cols_altered = ['time_counter', 'Altered votemper', 'Altered vosaline', 'δ18Ow altered', 'Altered δ18Oc estimated']

    df_original = df[[col for col in cols_original if col in df.columns]].copy()
    df_altered = df[[col for col in cols_altered if col in df.columns]].copy()

    if 'time_counter' in df_original.columns:
        df_original['time_counter'] = pd.to_datetime(df_original['time_counter'])
    if 'time_counter' in df_altered.columns:
        df_altered['time_counter'] = pd.to_datetime(df_altered['time_counter'])

    output_path = filename if filename.endswith(".xlsx") else filename + ".xlsx"
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        df_original.to_excel(writer, sheet_name="Original", index=False)
        df_altered.to_excel(writer, sheet_name="Altered", index=False)

    print(f"✅ Export completed: {output_path}")
    print(f"   - Sheet 'Original': unaltered series")
    print(f"   - Sheet 'Altered' : series with modified seasonality")


if __name__ == "__main__":
    main()
