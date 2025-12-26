#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module 05A: ΔSOS regression for timing/shape decomposition outputs

Regressions (pixel-wise):
  - TRshape  ~ ΔSOS
  - TRtiming ~ ΔSOS
  - TRsos    ~ ΔSOS
  - TRpos    ~ ΔSOS
  - Trate    ~ ΔSOS (Trate = TRc / GLS within SOS–POS)

Notes:
  - ΔSOS = SOS_year - SOSav (standard anomaly; advance < 0)
  - Uses outputs from 03_decomposition_timing_shape.py
"""

import numpy as np
import rasterio
from pathlib import Path
from tqdm import tqdm
from scipy import stats

# ==================== Global config ====================
ROOT = Path(r"I:\F\Data4")
ANALYSIS_DIR = ROOT / "Wang2025_Analysis"
DECOMP_DIR = ANALYSIS_DIR / "Decomposition_TimingShape"
TRC_DIR = ANALYSIS_DIR / "TRc_annual"
PHENO_DIR = ROOT / "Phenology_Output_1" / "GPP_phenology"
OUTPUT_DIR = ANALYSIS_DIR / "Statistical_Analysis_TimingShape"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

YEAR_START = 1982
YEAR_END = 2018
NODATA_OUT = -9999.0
NODATA_ABS_MAX = 1e20


def _is_valid_value(value, nodata):
    if nodata is None:
        return np.isfinite(value) & (np.abs(value) < NODATA_ABS_MAX)
    if np.isnan(nodata):
        return np.isfinite(value) & (np.abs(value) < NODATA_ABS_MAX)
    return (value != nodata) & np.isfinite(value) & (np.abs(value) < NODATA_ABS_MAX)


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


def linear_regression_maps(x_stack, y_stack, min_frac=0.6):
    n_years = x_stack.shape[0]
    valid = np.isfinite(x_stack) & np.isfinite(y_stack)
    n_valid = np.sum(valid, axis=0).astype(np.float32)

    x = np.where(valid, x_stack, np.nan)
    y = np.where(valid, y_stack, np.nan)

    mean_x = np.nanmean(x, axis=0)
    mean_y = np.nanmean(y, axis=0)
    dx = x - mean_x
    dy = y - mean_y

    cov_xy = np.nanmean(dx * dy, axis=0)
    var_x = np.nanmean(dx * dx, axis=0)
    var_y = np.nanmean(dy * dy, axis=0)

    with np.errstate(divide="ignore", invalid="ignore"):
        slope = cov_xy / var_x
        r = cov_xy / np.sqrt(var_x * var_y)
        r2 = r ** 2
        df = n_valid - 2
        t_stat = r * np.sqrt(df / (1 - r2))
        p_val = 2 * (1 - stats.t.cdf(np.abs(t_stat), df=df))

    invalid = (
        (n_valid < n_years * min_frac) |
        (n_valid < 3) |
        ~np.isfinite(slope) |
        ~np.isfinite(r2) |
        ~np.isfinite(p_val)
    )

    slope_map = np.full(slope.shape, NODATA_OUT, dtype=np.float32)
    pvalue_map = np.full(p_val.shape, NODATA_OUT, dtype=np.float32)
    r2_map = np.full(r2.shape, NODATA_OUT, dtype=np.float32)

    valid_mask = ~invalid
    slope_map[valid_mask] = slope[valid_mask]
    pvalue_map[valid_mask] = p_val[valid_mask]
    r2_map[valid_mask] = r2[valid_mask]

    return slope_map, pvalue_map, r2_map


def section_3_2_timing_shape(mask):
    print("\n" + "=" * 70)
    print("Section 3.2A: TRshape/TRtiming/TRsos/TRpos/Trate vs ΔSOS")
    print("=" * 70)

    years = list(range(YEAR_START, YEAR_END + 1))
    n_years = len(years)

    # Template
    data_first, profile, _ = read_geotiff(DECOMP_DIR / f"TRshape_{years[0]}.tif")
    height, width = data_first.shape

    # SOS/POS stack and SOSav
    print("\n[Step 1] 读取SOS/POS并计算SOSav...")
    sos_stack = []
    pos_stack = []
    for year in tqdm(years, desc="读取SOS/POS"):
        sos_data, _, sos_nodata = read_geotiff(PHENO_DIR / f"sos_gpp_{year}.tif")
        pos_data, _, pos_nodata = read_geotiff(PHENO_DIR / f"pos_doy_gpp_{year}.tif")
        sos_valid = np.where(_is_valid_value(sos_data, sos_nodata), sos_data, np.nan)
        pos_valid = np.where(_is_valid_value(pos_data, pos_nodata), pos_data, np.nan)
        sos_stack.append(sos_valid)
        pos_stack.append(pos_valid)
    sos_stack = np.stack(sos_stack, axis=0)
    pos_stack = np.stack(pos_stack, axis=0)
    sos_av = np.nanmean(sos_stack, axis=0)

    # ΔSOS (advance < 0)
    print("\n[Step 2] 计算 ΔSOS = SOS_year - SOSav...")
    delta_sos_stack = sos_stack - sos_av

    # GLS for Trate
    print("\n[Step 2b] 计算 GLS = POS - SOS + 1（用于Trate）...")
    gls_stack = np.full_like(sos_stack, np.nan, dtype=np.float32)
    for idx in range(n_years):
        sos_y = sos_stack[idx]
        pos_y = pos_stack[idx]
        sos_int = np.rint(sos_y).astype(np.int32)
        pos_int = np.rint(pos_y).astype(np.int32)
        valid = (
            np.isfinite(sos_y) & np.isfinite(pos_y) &
            (sos_int > 0) & (pos_int > 0) &
            (pos_int >= sos_int)
        )
        gls = np.full((height, width), np.nan, dtype=np.float32)
        gls[valid] = (pos_int[valid] - sos_int[valid] + 1).astype(np.float32)
        gls_stack[idx] = gls

    # Response variables
    print("\n[Step 3] 读取响应变量...")
    response_vars = ["TRshape", "TRtiming", "TRsos", "TRpos", "Trate"]
    response_data = {}
    for resp_var in response_vars:
        if resp_var == "Trate":
            continue
        data_stack = []
        for year in tqdm(years, desc=f"读取{resp_var}", leave=False):
            data, _, nodata = read_geotiff(DECOMP_DIR / f"{resp_var}_{year}.tif")
            data_stack.append(np.where(_is_valid_value(data, nodata), data, np.nan))
        response_data[resp_var] = np.stack(data_stack, axis=0)

    # Trate = TRc / GLS
    print("\n  计算 Trate = TRc / GLS...")
    trc_stack = []
    for year in tqdm(years, desc="读取TRc", leave=False):
        trc_data, _, trc_nodata = read_geotiff(TRC_DIR / f"TRc_{year}.tif")
        trc_stack.append(np.where(_is_valid_value(trc_data, trc_nodata), trc_data, np.nan))
    trc_stack = np.stack(trc_stack, axis=0)
    trate_stack = np.full_like(trc_stack, np.nan, dtype=np.float32)
    valid_trate = np.isfinite(trc_stack) & np.isfinite(gls_stack) & (gls_stack > 0)
    trate_stack[valid_trate] = trc_stack[valid_trate] / gls_stack[valid_trate]
    response_data["Trate"] = trate_stack

    # Regression
    print("\n[Step 4] 像元级线性回归: Response ~ ΔSOS...")
    output_dir = OUTPUT_DIR / "Section_3.2_Phenology_Impact"
    output_dir.mkdir(parents=True, exist_ok=True)

    for resp_var in response_vars:
        print(f"\n  分析: {resp_var} ~ ΔSOS")
        slope_map, pvalue_map, r_squared_map = linear_regression_maps(
            delta_sos_stack, response_data[resp_var], min_frac=0.6
        )

        if mask is not None:
            slope_map[~mask] = NODATA_OUT
            pvalue_map[~mask] = NODATA_OUT
            r_squared_map[~mask] = NODATA_OUT

        write_geotiff(output_dir / f"{resp_var}_vs_deltaSOS_slope.tif", slope_map, profile)
        write_geotiff(output_dir / f"{resp_var}_vs_deltaSOS_pvalue.tif", pvalue_map, profile)
        write_geotiff(output_dir / f"{resp_var}_vs_deltaSOS_R2.tif", r_squared_map, profile)

        sig_mask = (pvalue_map < 0.05) & (pvalue_map != NODATA_OUT)
        n_sig = np.sum(sig_mask)
        mean_slope = np.nanmean(slope_map[sig_mask]) if n_sig > 0 else np.nan
        std_slope = np.nanstd(slope_map[sig_mask]) if n_sig > 0 else np.nan

        print(f"    显著像元数 (p<0.05): {n_sig}")
        print(f"    平均斜率: {mean_slope:.3f} ± {std_slope:.3f} mm/day per day ΔSOS")

    print("\n  ✓ Section 3.2A 分析完成")
    print(f"  输出目录: {output_dir}")


def main():
    print("\n" + "=" * 70)
    print("Timing/Shape ΔSOS regression")
    print("=" * 70)

    mask_file = ANALYSIS_DIR / "masks" / "combined_mask.tif"
    if not mask_file.exists():
        print(f"  ⚠ 掩膜不存在: {mask_file}，使用全有效掩膜")
        data_sample, _, _ = read_geotiff(DECOMP_DIR / f"TRshape_{YEAR_START}.tif")
        mask = np.ones(data_sample.shape, dtype=bool)
    else:
        mask_data, _, mask_nodata = read_geotiff(mask_file)
        mask = _is_valid_value(mask_data, mask_nodata) & (mask_data > 0)
        print(f"  ✓ 掩膜读取成功，有效像元数: {np.sum(mask)}")

    section_3_2_timing_shape(mask)

    print("\n" + "=" * 70)
    print("✓ 执行完成！")
    print(f"输出目录: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
