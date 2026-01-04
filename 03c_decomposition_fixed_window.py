#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module 03C: Fixed-Window decomposition (Method 2)

方法2：实际窗口变化 + 固定窗口强度分解

Decompose:
  TRc_y = TRc_av + TR_window_change + TR_fixed_window

Where:
  TRc_y   = sum_{t=SOS_y..POS_y} B_y(t)  # 当年实际累积
  TRc_av  = sum_{t=SOSav..POSav} A(t)    # 气候态基线（固定窗口）

  TR_window_change = 窗口变化带来的实际额外累积
    = TRc_y_estimated - TRc_y_fixed_window
    （窗口外部分按当年强度估算）

  TR_fixed_window = 固定窗口内的强度差异
    = sum_{t=SOSav..POSav} [B_y(t) - A(t)]
    （在固定窗口内，当年与气候态的累积差异）

Optional split:
  TR_window_change = TR_sos_change + TR_pos_change
  TR_sos_change = SOS变化带来的实际累积（如果SOS_y < SOSav，为正）
  TR_pos_change = POS变化带来的实际累积（如果POS_y > POSav，为正）

核心思想：
  - 固定窗口[SOSav, POSav]内的强度变化可以剥离窗口选择效应
  - TR_fixed_window/(POSav-SOSav+1) = 真实的平均速率差异
  - 这是评估"蒸腾速率是否真的变高变低"的最佳指标

实现方法（基于现有数据）：
  由于没有当年每日TR数据B_y(t)，我们使用比例假设：
  1. 计算增强因子 = TRc_y / sum_{SOS_y..POS_y} A(t)
  2. 固定窗口估计 = 增强因子 × TRc_av
  3. 窗口变化部分 = TRc_y - 固定窗口估计

Notes:
  - A(t) is the climatology TR_daily_av (365 bands).
  - B_y(t) is the year-specific daily TR (estimated using proportional scaling).
  - 固定窗口为[SOSav, POSav]，所有年份一致
"""

import numpy as np
import rasterio
from pathlib import Path
from tqdm import tqdm

# 导入配置
from _config import (
    ROOT, PHENO_DIR, TRC_ANNUAL_DIR, CLIMATOLOGY_DIR, DECOMPOSITION_FIXED_DIR,
    YEAR_START, YEAR_END, NODATA_OUT
)

# 确保输出目录存在
DECOMPOSITION_FIXED_DIR.mkdir(parents=True, exist_ok=True)

# 向后兼容：保留旧变量名
TRC_DIR = TRC_ANNUAL_DIR
CLIM_DIR = CLIMATOLOGY_DIR
OUTPUT_DIR = DECOMPOSITION_FIXED_DIR


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
    cum_arr: shape (365, H, W) - cumulative sum along DOY axis
    start_doy, end_doy: shape (H, W) - DOY values (1-365)
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
    """
    方法2分解：固定窗口强度 + 实际窗口变化

    Returns:
    --------
    tr_window_change : 窗口变化带来的实际累积
    tr_fixed_window : 固定窗口内的强度差异
    tr_sos_change : SOS变化的贡献
    tr_pos_change : POS变化的贡献
    """
    # 读取当年数据
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

    # 有效性检查（与03a/03b一致）
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

    # 初始化输出
    tr_window_change = np.full((height, width), NODATA_OUT, dtype=np.float32)
    tr_fixed_window = np.full((height, width), NODATA_OUT, dtype=np.float32)
    tr_sos_change = np.full((height, width), NODATA_OUT, dtype=np.float32)
    tr_pos_change = np.full((height, width), NODATA_OUT, dtype=np.float32)

    if not np.any(valid):
        return tr_window_change, tr_fixed_window, tr_sos_change, tr_pos_change

    # Clip DOY to 1-365
    sos_y_i = np.clip(np.rint(sos_y).astype(np.int32), 1, 365)
    pos_y_i = np.clip(np.rint(pos_y).astype(np.int32), 1, 365)
    sos_av_i = np.clip(np.rint(sos_av).astype(np.int32), 1, 365)
    pos_av_i = np.clip(np.rint(pos_av).astype(np.int32), 1, 365)

    # 进一步检查窗口有效性
    valid_sos = valid & (pos_y_i >= sos_y_i) & (pos_av_i >= sos_av_i)

    if not np.any(valid_sos):
        return tr_window_change, tr_fixed_window, tr_sos_change, tr_pos_change

    # ========== 核心计算 ==========

    # 1. 计算当年窗口在气候态下的理论累积
    #    sum_{SOS_y to POS_y} A(t)
    trc_y_clim = _sum_cum_range(tr_cum, sos_y_i, pos_y_i)

    # 2. 计算增强因子（当年相对于气候态的倍数）
    #    增强因子 = TRc_y / sum_{SOS_y to POS_y} A(t)
    with np.errstate(divide='ignore', invalid='ignore'):
        enhancement_factor = np.where(
            (trc_y_clim > 1e-6) & valid_sos,
            trc_y / trc_y_clim,
            1.0
        )

    # 3. 估算固定窗口内的当年累积
    #    TRc_y_fixed = 增强因子 × TRc_av
    #    假设：相对增强在空间上是均匀的
    trc_y_fixed = enhancement_factor * trc_av

    # 4. 固定窗口强度差异
    #    TR_fixed_window = TRc_y_fixed - TRc_av
    tr_fixed_window[valid_sos] = (trc_y_fixed - trc_av)[valid_sos]

    # 5. 窗口变化部分（残差）
    #    TR_window_change = TRc_y - TRc_y_fixed
    tr_window_change[valid_sos] = (trc_y - trc_y_fixed)[valid_sos]

    # 6. 进一步分解窗口变化为SOS和POS部分
    #    这部分需要更详细的估算
    #    注意：避免np.where包裹_sum_cum_range，防止边界越界

    # SOS变化效应估算
    # 如果SOS_y < SOSav（提前），计算[SOS_y, SOSav-1]的贡献
    # 如果SOS_y > SOSav（推迟），计算[SOSav, SOS_y-1]的贡献（为负）
    sos_before = sos_y_i < sos_av_i  # SOS提前
    sos_after = sos_y_i > sos_av_i   # SOS推迟

    # 初始化SOS变化数组
    sos_advance_clim = np.zeros((height, width), dtype=np.float32)
    sos_delay_clim = np.zeros((height, width), dtype=np.float32)

    # SOS提前：估算[SOS_y, SOSav-1]在当年强度下的累积
    # 确保端点有效：sos_av_i > sos_y_i（即sos_av_i >= sos_y_i + 1，所以sos_av_i - 1 >= sos_y_i）
    valid_advance = valid_sos & sos_before & (sos_av_i > sos_y_i)
    if np.any(valid_advance):
        sos_advance_clim[valid_advance] = _sum_cum_range(
            tr_cum, sos_y_i, sos_av_i - 1
        )[valid_advance]

    # SOS推迟：估算[SOSav, SOS_y-1]在当年强度下的累积（为负贡献）
    # 确保端点有效：sos_y_i > sos_av_i
    valid_delay = valid_sos & sos_after & (sos_y_i > sos_av_i)
    if np.any(valid_delay):
        sos_delay_clim[valid_delay] = _sum_cum_range(
            tr_cum, sos_av_i, sos_y_i - 1
        )[valid_delay]

    # 应用增强因子
    tr_sos_advance = enhancement_factor * sos_advance_clim
    tr_sos_delay = enhancement_factor * sos_delay_clim

    # 根据提前/推迟设置正负号
    tr_sos_change[valid_sos] = np.where(
        sos_before,
        tr_sos_advance,
        -tr_sos_delay
    )[valid_sos]

    # POS变化效应估算（类似SOS）
    pos_after = pos_y_i > pos_av_i   # POS延后
    pos_before = pos_y_i < pos_av_i  # POS提前

    # 初始化POS变化数组
    pos_extend_clim = np.zeros((height, width), dtype=np.float32)
    pos_shorten_clim = np.zeros((height, width), dtype=np.float32)

    # POS延后：估算[POSav+1, POS_y]在当年强度下的累积
    # 确保端点有效：pos_av_i + 1 <= 365 且 pos_y_i >= pos_av_i + 1
    valid_extend = valid_sos & pos_after & (pos_av_i < 365) & (pos_y_i > pos_av_i)
    if np.any(valid_extend):
        pos_extend_clim[valid_extend] = _sum_cum_range(
            tr_cum, pos_av_i + 1, pos_y_i
        )[valid_extend]

    # POS提前：估算[POS_y+1, POSav]在当年强度下的累积（为负贡献）
    # 确保端点有效：pos_y_i + 1 <= 365 且 pos_av_i >= pos_y_i + 1
    valid_shorten = valid_sos & pos_before & (pos_y_i < 365) & (pos_av_i > pos_y_i)
    if np.any(valid_shorten):
        pos_shorten_clim[valid_shorten] = _sum_cum_range(
            tr_cum, pos_y_i + 1, pos_av_i
        )[valid_shorten]

    # 应用增强因子
    tr_pos_extend = enhancement_factor * pos_extend_clim
    tr_pos_shorten = enhancement_factor * pos_shorten_clim

    # 根据延后/提前设置正负号
    tr_pos_change[valid_sos] = np.where(
        pos_after,
        tr_pos_extend,
        -tr_pos_shorten
    )[valid_sos]

    return tr_window_change, tr_fixed_window, tr_sos_change, tr_pos_change


def main():
    print("\n" + "=" * 80)
    print("Fixed-Window Decomposition (Method 2)")
    print("=" * 80)
    print("\n核心思想：")
    print("  - 固定窗口[SOSav, POSav]内的强度变化 → 剥离窗口选择效应")
    print("  - TR_fixed_window/(POSav-SOSav+1) → 真实的平均速率差异")
    print("  - 窗口变化部分 → 窗口延长/缩短带来的实际累积\n")

    years = list(range(YEAR_START, YEAR_END + 1))

    # Load mask
    mask_file = ROOT / "Wang2025_Analysis" / "masks" / "combined_mask.tif"
    if not mask_file.exists():
        raise FileNotFoundError(f"Missing mask: {mask_file}")
    mask, mask_profile, mask_nodata = read_geotiff(mask_file)
    mask = _is_valid_value(mask, mask_nodata) & (mask > 0)

    # Load climatology
    tr_daily_av, sos_av, pos_av, nodata_sos, nodata_pos, profile = load_climatology()

    # 空间一致性检查
    check_spatial_consistency(profile, mask_file, mask_profile, "combined_mask")
    tr_valid = np.any(np.isfinite(tr_daily_av), axis=0)
    tr_cum = np.nancumsum(tr_daily_av, axis=0).astype(np.float32)

    # Compute TRc_av (固定窗口基线)
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

    # 保存固定窗口长度（用于后续速率计算）
    fixed_window_length = np.full(sos_av.shape, NODATA_OUT, dtype=np.float32)
    fixed_window_length[valid_av] = (pos_av_i - sos_av_i + 1)[valid_av]
    write_geotiff(OUTPUT_DIR / "Fixed_Window_Length.tif", fixed_window_length, profile)

    # Decompose each year
    for year in tqdm(years, desc="Fixed-Window Decomposition"):
        tr_window_change, tr_fixed_window, tr_sos_change, tr_pos_change = decompose_year(
            year, tr_cum, trc_av, sos_av, pos_av,
            nodata_sos, nodata_pos, mask, tr_valid, profile
        )
        write_geotiff(OUTPUT_DIR / f"TR_window_change_{year}.tif", tr_window_change, profile)
        write_geotiff(OUTPUT_DIR / f"TR_fixed_window_{year}.tif", tr_fixed_window, profile)
        write_geotiff(OUTPUT_DIR / f"TR_sos_change_{year}.tif", tr_sos_change, profile)
        write_geotiff(OUTPUT_DIR / f"TR_pos_change_{year}.tif", tr_pos_change, profile)

    # Quality check
    print("\n[Quality Check]")
    test_year = 2000
    trc_test, _, nodata_trc_test = read_geotiff(TRC_DIR / f"TRc_{test_year}.tif")
    tr_window_test, _, _ = read_geotiff(OUTPUT_DIR / f"TR_window_change_{test_year}.tif")
    tr_fixed_test, _, _ = read_geotiff(OUTPUT_DIR / f"TR_fixed_window_{test_year}.tif")

    valid_check = (
        mask &
        _is_valid_value(trc_test, nodata_trc_test) &
        (trc_av != NODATA_OUT) &
        (tr_window_test != NODATA_OUT) &
        (tr_fixed_test != NODATA_OUT)
    )

    # 残差检验: TRc_y - TRc_av - TR_window_change - TR_fixed_window
    residual = (trc_test[valid_check] - trc_av[valid_check] -
                tr_window_test[valid_check] - tr_fixed_test[valid_check])

    print(f"\n  Test year: {test_year}")
    print(f"  Residual (TRc - TRc_av - TR_window_change - TR_fixed_window):")
    print(f"    Mean: {np.mean(residual):.6f} mm")
    print(f"    Std:  {np.std(residual):.6f} mm")
    print(f"    Max:  {np.max(np.abs(residual)):.6f} mm")
    print(f"    95th percentile: {np.percentile(np.abs(residual), 95):.6f} mm")
    print(f"    99th percentile: {np.percentile(np.abs(residual), 99):.6f} mm")

    if np.max(np.abs(residual)) < 1e-2:
        print("  ✓ Quality check passed! Residual negligible.")
    else:
        print("  ⚠ Warning: Large residual detected, please check calculation.")

    print("\n✓ Fixed-Window decomposition complete.")
    print(f"Output: {OUTPUT_DIR}")
    print("\n输出文件说明：")
    print("  - TRc_av.tif: 固定窗口基线")
    print("  - Fixed_Window_Length.tif: 固定窗口长度（天数）")
    print("  - TR_fixed_window_YYYY.tif: 固定窗口内的强度差异（核心指标）")
    print("  - TR_window_change_YYYY.tif: 窗口变化带来的实际累积")
    print("  - TR_sos_change_YYYY.tif: SOS变化的贡献")
    print("  - TR_pos_change_YYYY.tif: POS变化的贡献")
    print("\n速率计算：")
    print("  Fixed_Trate = TR_fixed_window / Fixed_Window_Length")


if __name__ == "__main__":
    main()
