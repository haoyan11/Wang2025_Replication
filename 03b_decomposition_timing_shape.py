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

# 导入配置
from _config import (
    ROOT, PHENO_DIR, TRC_ANNUAL_DIR, CLIMATOLOGY_DIR, DECOMPOSITION_TIMING_DIR,
    YEAR_START, YEAR_END, NODATA_OUT
)

# 确保输出目录存在
DECOMPOSITION_TIMING_DIR.mkdir(parents=True, exist_ok=True)

# 向后兼容：保留旧变量名
TRC_DIR = TRC_ANNUAL_DIR
CLIM_DIR = CLIMATOLOGY_DIR
OUTPUT_DIR = DECOMPOSITION_TIMING_DIR


def _is_valid_value(value, nodata):
    """
    判断值是否有效（非NODATA）
    兼容各种NODATA类型：-9999, NaN, 3.4e38等
    当nodata元信息缺失时，使用绝对值阈值过滤极大填充值
    """
    # 基础检查：有限性
    valid = np.isfinite(value)

    # 如果有明确的nodata值，检查是否等于nodata
    if nodata is not None and not np.isnan(nodata):
        valid = valid & (value != nodata)

    # 额外保护：过滤常见的极大填充值（1e20, 3.4e38等）
    valid = valid & (np.abs(value) < 1e10)

    return valid


def read_geotiff(file_path):
    with rasterio.open(file_path) as src:
        data = src.read(1)
        profile = src.profile.copy()
        nodata = src.nodata
    return data, profile, nodata


def write_geotiff(file_path, data, profile):
    """写入单波段GeoTIFF（不修改原profile）"""
    out_profile = profile.copy()
    out_profile.update(dtype=rasterio.float32, count=1, compress="lzw", nodata=NODATA_OUT)
    with rasterio.open(file_path, "w", **out_profile) as dst:
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


def check_spatial_consistency(ref_profile, file_path, data_profile, var_name="data"):
    """检查栅格空间一致性（shape/crs/transform）"""
    if ref_profile['width'] != data_profile['width'] or \
       ref_profile['height'] != data_profile['height']:
        raise ValueError(
            f"Shape mismatch for {var_name}: {file_path}\n"
            f"  Expected: {ref_profile['height']}x{ref_profile['width']}\n"
            f"  Got: {data_profile['height']}x{data_profile['width']}"
        )
    if ref_profile['crs'] != data_profile['crs']:
        raise ValueError(
            f"CRS mismatch for {var_name}: {file_path}\n"
            f"  Expected: {ref_profile['crs']}\n"
            f"  Got: {data_profile['crs']}"
        )

    # Transform比较使用容差（避免浮点精度误报）
    ref_transform = ref_profile['transform']
    data_transform = data_profile['transform']
    if ref_transform is not None and data_transform is not None:
        transform_match = all(
            abs(ref_transform[i] - data_transform[i]) < 1e-6
            for i in range(6)
        )
        if not transform_match:
            raise ValueError(
                f"Transform mismatch for {var_name}: {file_path}\n"
                f"  Expected: {ref_transform}\n"
                f"  Got: {data_transform}"
            )


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

    sos_av, sos_profile, nodata_sos = read_geotiff(sos_av_file)
    check_spatial_consistency(profile, sos_av_file, sos_profile, "SOSav")

    pos_av, pos_profile, nodata_pos = read_geotiff(pos_av_file)
    check_spatial_consistency(profile, pos_av_file, pos_profile, "POSav")

    return tr_daily_av, sos_av, pos_av, nodata_sos, nodata_pos, profile


def decompose_year(year, tr_cum, trc_av, sos_av, pos_av,
                   nodata_sos, nodata_pos, mask, tr_valid, ref_profile):
    trc_file = TRC_DIR / f"TRc_{year}.tif"
    sos_file = PHENO_DIR / f"sos_gpp_{year}.tif"
    pos_file = PHENO_DIR / f"pos_doy_gpp_{year}.tif"

    trc_y, trc_profile, nodata_trc = read_geotiff(trc_file)
    sos_y, sos_profile, nodata_sos_y = read_geotiff(sos_file)
    pos_y, pos_profile, nodata_pos_y = read_geotiff(pos_file)

    # 空间一致性检查
    check_spatial_consistency(ref_profile, trc_file, trc_profile, f"TRc_{year}")
    check_spatial_consistency(ref_profile, sos_file, sos_profile, f"sos_gpp_{year}")
    check_spatial_consistency(ref_profile, pos_file, pos_profile, f"pos_doy_gpp_{year}")

    height, width = trc_y.shape

    valid = (
        _is_valid_value(trc_y, nodata_trc) &
        (trc_av != NODATA_OUT) &
        _is_valid_value(sos_y, nodata_sos_y) &
        _is_valid_value(pos_y, nodata_pos_y) &
        _is_valid_value(sos_av, nodata_sos) &
        _is_valid_value(pos_av, nodata_pos) &
        (sos_y > 0) & (sos_y <= 366) &  # 允许闰年DOY=366
        (pos_y > 0) & (pos_y <= 366) &  # 允许闰年DOY=366
        (sos_av > 0) & (sos_av <= 366) &  # 允许闰年DOY=366
        (pos_av > 0) & (pos_av <= 366) &  # 允许闰年DOY=366
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
    mask, mask_profile, mask_nodata = read_geotiff(mask_file)
    mask = _is_valid_value(mask, mask_nodata) & (mask > 0)

    # Climatology
    tr_daily_av, sos_av, pos_av, nodata_sos, nodata_pos, profile = load_climatology()

    # 空间一致性检查
    check_spatial_consistency(profile, mask_file, mask_profile, "combined_mask")
    tr_valid = np.any(np.isfinite(tr_daily_av), axis=0)
    tr_cum = np.nancumsum(tr_daily_av, axis=0).astype(np.float32)

    # TRc_av (A window over SOSav..POSav)
    sos_av_i = np.clip(np.rint(sos_av).astype(np.int32), 1, 365)
    pos_av_i = np.clip(np.rint(pos_av).astype(np.int32), 1, 365)
    trc_av = np.full(sos_av.shape, NODATA_OUT, dtype=np.float32)
    valid_av = (
        _is_valid_value(sos_av, nodata_sos) &
        _is_valid_value(pos_av, nodata_pos) &
        (sos_av > 0) & (sos_av <= 366) &  # 允许闰年DOY=366，排除异常值
        (pos_av > 0) & (pos_av <= 366) &  # 允许闰年DOY=366，排除异常值
        (pos_av_i >= sos_av_i) &
        mask & tr_valid
    )
    trc_av[valid_av] = _sum_cum_range(tr_cum, sos_av_i, pos_av_i)[valid_av]
    write_geotiff(OUTPUT_DIR / "TRc_av.tif", trc_av, profile)

    # Decompose each year
    for year in tqdm(years, desc="Decomposition"):
        tr_timing, tr_shape, tr_sos, tr_pos = decompose_year(
            year, tr_cum, trc_av, sos_av, pos_av,
            nodata_sos, nodata_pos, mask, tr_valid, profile
        )
        write_geotiff(OUTPUT_DIR / f"TRtiming_{year}.tif", tr_timing, profile)
        write_geotiff(OUTPUT_DIR / f"TRshape_{year}.tif", tr_shape, profile)
        write_geotiff(OUTPUT_DIR / f"TRsos_{year}.tif", tr_sos, profile)
        write_geotiff(OUTPUT_DIR / f"TRpos_{year}.tif", tr_pos, profile)

    # Quality check
    print("\n[Quality Check]")
    test_year = 2000
    trc_test, _, nodata_trc_test = read_geotiff(TRC_DIR / f"TRc_{test_year}.tif")
    tr_timing_test, _, _ = read_geotiff(OUTPUT_DIR / f"TRtiming_{test_year}.tif")
    tr_shape_test, _, _ = read_geotiff(OUTPUT_DIR / f"TRshape_{test_year}.tif")

    valid_check = (
        mask &
        _is_valid_value(trc_test, nodata_trc_test) &
        (trc_av != NODATA_OUT) &
        (tr_timing_test != NODATA_OUT) &
        (tr_shape_test != NODATA_OUT)
    )

    # 残差检验: TRc_y - TRc_av - TR_timing - TR_shape
    residual = (trc_test[valid_check] - trc_av[valid_check] -
                tr_timing_test[valid_check] - tr_shape_test[valid_check])

    print(f"\n  Test year: {test_year}")
    print(f"  Residual (TRc - TRc_av - TRtiming - TRshape):")
    print(f"    Mean: {np.mean(residual):.6f} mm")
    print(f"    Std:  {np.std(residual):.6f} mm")
    print(f"    Max:  {np.max(np.abs(residual)):.6f} mm")
    print(f"    95th percentile: {np.percentile(np.abs(residual), 95):.6f} mm")
    print(f"    99th percentile: {np.percentile(np.abs(residual), 99):.6f} mm")

    if np.max(np.abs(residual)) < 1e-2:
        print("  ✓ Quality check passed! Residual negligible.")
    else:
        print("  ⚠ Warning: Large residual detected, please check calculation.")

    print("\n✓ Timing/Shape decomposition complete.")
    print(f"Output: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
