#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module 03A: Timing/Shape decomposition (A1 scheme)

Decompose:
  TRc_y - TRc_av = ΔTR_timing + ΔTR_shape

Where:
  TRc_y   = sum_{t=SOS_y..POS_y} B_y(t)
  TRc_av  = sum_{t=SOS_av..POS_av} A(t)
  ΔTR_timing = sum_{t=SOS_y..POS_y} A(t) - sum_{t=SOS_av..POS_av} A(t)
  ΔTR_shape  = sum_{t=SOS_y..POS_y} (B_y(t) - A(t))

Optional split:
  ΔTR_timing = ΔTR_SOS + ΔTR_POS
  ΔTR_SOS = sum_{t=SOS_y..POS_av} A(t) - sum_{t=SOS_av..POS_av} A(t)
  ΔTR_POS = sum_{t=SOS_y..POS_y} A(t) - sum_{t=SOS_y..POS_av} A(t)

Notes:
  - A(t) is the climatology TR_daily_av (365 bands).
  - B_y(t) is the year-specific daily TR; we use TRc_y from 02_TRc_calculation.py.
  - Windows are inclusive (SOS <= t <= POS).
"""

import numpy as np
import rasterio
from pathlib import Path
from tqdm import tqdm

# ==================== Global config ====================
ROOT = Path(r"I:\F\Data4")
PHENO_DIR = ROOT / "Phenology_Output_1" / "GPP_phenology"
TRC_DIR = ROOT / "Wang2025_Analysis" / "TRc_annual"
CLIM_DIR = ROOT / "Wang2025_Analysis" / "Climatology"
OUTPUT_DIR = ROOT / "Wang2025_Analysis" / "Decomposition_TimingShape"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

YEAR_START = 1982
YEAR_END = 2018
NODATA_OUT = -9999.0


def _is_valid_value(value, nodata):
    if nodata is None:
        return np.isfinite(value)
    if np.isnan(nodata):
        return np.isfinite(value)
    return (value != nodata) & np.isfinite(value)


def read_geotiff(file_path):
    with rasterio.open(file_path) as src:
        data = src.read(1)
        profile = src.profile.copy()
        nodata = src.nodata
    return data, profile, nodata


def write_geotiff(file_path, data, profile):
    profile.update(dtype=rasterio.float32, count=1, compress="lzw", nodata=NODATA_OUT)
    with rasterio.open(file_path, "w", **profile) as dst:
        dst.write(data.astype(np.float32), 1)


def _sum_cum_range(cum_arr, start_doy, end_doy):
    """
    Sum over inclusive [start, end], DOY is 1-based.
    """
    flat = cum_arr.reshape(cum_arr.shape[0], -1)
    start = start_doy.ravel()
    end = end_doy.ravel()
    idx = np.arange(start.size)
    end_vals = flat[end - 1, idx]
    start_vals = np.where(start > 1, flat[start - 2, idx], 0.0)
    return (end_vals - start_vals).reshape(start_doy.shape)


def load_climatology():
    tr_clim_file = CLIM_DIR / "TR_daily_climatology.tif"
    sos_av_file = CLIM_DIR / "SOSav.tif"
    pos_av_file = CLIM_DIR / "POSav.tif"

    if not tr_clim_file.exists():
        raise FileNotFoundError(f"Missing climatology: {tr_clim_file}")
    if not sos_av_file.exists():
        raise FileNotFoundError(f"Missing SOSav: {sos_av_file}")
    if not pos_av_file.exists():
        raise FileNotFoundError(f"Missing POSav: {pos_av_file}")

    with rasterio.open(tr_clim_file) as src:
        tr_daily_av = src.read().astype(np.float32)
        profile = src.profile.copy()
        nodata_tr = src.nodata

    if nodata_tr is None or np.isnan(nodata_tr):
        tr_daily_av[~np.isfinite(tr_daily_av)] = np.nan
    else:
        tr_daily_av[tr_daily_av == nodata_tr] = np.nan
    tr_daily_av[tr_daily_av < 0] = np.nan

    sos_av, _, nodata_sos = read_geotiff(sos_av_file)
    pos_av, _, nodata_pos = read_geotiff(pos_av_file)

    return tr_daily_av, sos_av, pos_av, nodata_sos, nodata_pos, profile


def decompose_year(year, tr_cum, trc_av, sos_av, pos_av,
                   nodata_sos, nodata_pos, mask, tr_valid):
    trc_y, _, nodata_trc = read_geotiff(TRC_DIR / f"TRc_{year}.tif")
    sos_y, _, nodata_sos_y = read_geotiff(PHENO_DIR / f"sos_gpp_{year}.tif")
    pos_y, _, nodata_pos_y = read_geotiff(PHENO_DIR / f"pos_doy_gpp_{year}.tif")

    height, width = trc_y.shape

    valid = (
        _is_valid_value(trc_y, nodata_trc) &
        (trc_av != NODATA_OUT) &
        _is_valid_value(sos_y, nodata_sos_y) &
        _is_valid_value(pos_y, nodata_pos_y) &
        _is_valid_value(sos_av, nodata_sos) &
        _is_valid_value(pos_av, nodata_pos) &
        (sos_y > 0) & (pos_y > 0) &
        (sos_av > 0) & (pos_av > 0) &
        (pos_y >= sos_y) &
        (pos_av >= sos_av) &
        mask &
        tr_valid
    )

    # Round and clip to 1..365
    sos_y_i = np.clip(np.rint(sos_y).astype(np.int32), 1, 365)
    pos_y_i = np.clip(np.rint(pos_y).astype(np.int32), 1, 365)
    sos_av_i = np.clip(np.rint(sos_av).astype(np.int32), 1, 365)
    pos_av_i = np.clip(np.rint(pos_av).astype(np.int32), 1, 365)

    # Sum A(t) over SOS_y..POS_y
    sum_a_y = _sum_cum_range(tr_cum, sos_y_i, pos_y_i)

    # Timing and shape components
    tr_timing = np.full((height, width), NODATA_OUT, dtype=np.float32)
    tr_shape = np.full((height, width), NODATA_OUT, dtype=np.float32)
    tr_timing[valid] = (sum_a_y - trc_av)[valid]
    tr_shape[valid] = (trc_y - sum_a_y)[valid]

    # Optional split: SOS and POS timing
    tr_sos = np.full((height, width), NODATA_OUT, dtype=np.float32)
    tr_pos = np.full((height, width), NODATA_OUT, dtype=np.float32)
    valid_sos = valid & (pos_av_i >= sos_y_i)
    sum_a_sos = _sum_cum_range(tr_cum, sos_y_i, pos_av_i)
    tr_sos[valid_sos] = (sum_a_sos - trc_av)[valid_sos]
    tr_pos[valid_sos] = (sum_a_y - sum_a_sos)[valid_sos]

    return tr_timing, tr_shape, tr_sos, tr_pos


def main():
    print("\n" + "=" * 70)
    print("Timing/Shape decomposition (A1 scheme)")
    print("=" * 70)

    years = list(range(YEAR_START, YEAR_END + 1))

    # Mask
    mask_file = ROOT / "Wang2025_Analysis" / "masks" / "combined_mask.tif"
    if not mask_file.exists():
        raise FileNotFoundError(f"Missing mask: {mask_file}")
    mask, profile, mask_nodata = read_geotiff(mask_file)
    mask = _is_valid_value(mask, mask_nodata) & (mask > 0)

    # Climatology
    tr_daily_av, sos_av, pos_av, nodata_sos, nodata_pos, profile = load_climatology()
    tr_valid = np.any(np.isfinite(tr_daily_av), axis=0)
    tr_cum = np.nancumsum(tr_daily_av, axis=0).astype(np.float32)

    # TRc_av (A window over SOSav..POSav)
    sos_av_i = np.clip(np.rint(sos_av).astype(np.int32), 1, 365)
    pos_av_i = np.clip(np.rint(pos_av).astype(np.int32), 1, 365)
    trc_av = np.full(sos_av.shape, NODATA_OUT, dtype=np.float32)
    valid_av = (
        _is_valid_value(sos_av, nodata_sos) &
        _is_valid_value(pos_av, nodata_pos) &
        (pos_av_i >= sos_av_i) &
        mask & tr_valid
    )
    trc_av[valid_av] = _sum_cum_range(tr_cum, sos_av_i, pos_av_i)[valid_av]
    write_geotiff(OUTPUT_DIR / "TRc_av.tif", trc_av, profile)

    # Decompose each year
    for year in tqdm(years, desc="Decomposition"):
        tr_timing, tr_shape, tr_sos, tr_pos = decompose_year(
            year, tr_cum, trc_av, sos_av, pos_av,
            nodata_sos, nodata_pos, mask, tr_valid
        )
        write_geotiff(OUTPUT_DIR / f"TRtiming_{year}.tif", tr_timing, profile)
        write_geotiff(OUTPUT_DIR / f"TRshape_{year}.tif", tr_shape, profile)
        write_geotiff(OUTPUT_DIR / f"TRsos_{year}.tif", tr_sos, profile)
        write_geotiff(OUTPUT_DIR / f"TRpos_{year}.tif", tr_pos, profile)

    print("\n✓ Timing/Shape decomposition complete.")
    print(f"Output: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
