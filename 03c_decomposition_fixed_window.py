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

# ==================== Global config ====================
ROOT = Path(r"I:\F\Data4")
PHENO_DIR = ROOT / "Phenology_Output_1" / "GPP_phenology"
TRC_DIR = ROOT / "Wang2025_Analysis" / "TRc_annual"
CLIM_DIR = ROOT / "Wang2025_Analysis" / "Climatology"
OUTPUT_DIR = ROOT / "Wang2025_Analysis" / "Decomposition_FixedWindow"
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
    trc_y, _, nodata_trc = read_geotiff(TRC_DIR / f"TRc_{year}.tif")
    sos_y, _, nodata_sos_y = read_geotiff(PHENO_DIR / f"sos_gpp_{year}.tif")
    pos_y, _, nodata_pos_y = read_geotiff(PHENO_DIR / f"pos_doy_gpp_{year}.tif")

    height, width = trc_y.shape

    # 有效性检查（与03b一致）
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

    # SOS变化效应估算
    # 如果SOS_y < SOSav（提前），计算[SOS_y, SOSav-1]的贡献
    # 如果SOS_y > SOSav（推迟），计算[SOSav, SOS_y-1]的贡献（为负）
    sos_before = sos_y_i < sos_av_i  # SOS提前
    sos_after = sos_y_i > sos_av_i   # SOS推迟

    # SOS提前：估算[SOS_y, SOSav-1]在当年强度下的累积
    # 近似：使用增强因子 × 气候态累积
    sos_advance_clim = np.where(
        sos_before,
        _sum_cum_range(tr_cum, sos_y_i, sos_av_i - 1),
        0.0
    )
    tr_sos_advance = enhancement_factor * sos_advance_clim

    # SOS推迟：估算[SOSav, SOS_y-1]在当年强度下的累积（为负贡献）
    sos_delay_clim = np.where(
        sos_after,
        _sum_cum_range(tr_cum, sos_av_i, sos_y_i - 1),
        0.0
    )
    tr_sos_delay = enhancement_factor * sos_delay_clim

    tr_sos_change[valid_sos] = np.where(
        sos_before,
        tr_sos_advance,
        -tr_sos_delay
    )[valid_sos]

    # POS变化效应估算（类似SOS）
    pos_after = pos_y_i > pos_av_i   # POS延后
    pos_before = pos_y_i < pos_av_i  # POS提前

    pos_extend_clim = np.where(
        pos_after,
        _sum_cum_range(tr_cum, pos_av_i + 1, pos_y_i),
        0.0
    )
    tr_pos_extend = enhancement_factor * pos_extend_clim

    pos_shorten_clim = np.where(
        pos_before,
        _sum_cum_range(tr_cum, pos_y_i + 1, pos_av_i),
        0.0
    )
    tr_pos_shorten = enhancement_factor * pos_shorten_clim

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
    mask, profile, mask_nodata = read_geotiff(mask_file)
    mask = _is_valid_value(mask, mask_nodata) & (mask > 0)

    # Load climatology
    tr_daily_av, sos_av, pos_av, nodata_sos, nodata_pos, profile = load_climatology()
    tr_valid = np.any(np.isfinite(tr_daily_av), axis=0)
    tr_cum = np.nancumsum(tr_daily_av, axis=0).astype(np.float32)

    # Compute TRc_av (固定窗口基线)
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

    # 保存固定窗口长度（用于后续速率计算）
    fixed_window_length = np.full(sos_av.shape, NODATA_OUT, dtype=np.float32)
    fixed_window_length[valid_av] = (pos_av_i - sos_av_i + 1)[valid_av]
    write_geotiff(OUTPUT_DIR / "Fixed_Window_Length.tif", fixed_window_length, profile)

    # Decompose each year
    for year in tqdm(years, desc="Fixed-Window Decomposition"):
        tr_window_change, tr_fixed_window, tr_sos_change, tr_pos_change = decompose_year(
            year, tr_cum, trc_av, sos_av, pos_av,
            nodata_sos, nodata_pos, mask, tr_valid
        )
        write_geotiff(OUTPUT_DIR / f"TR_window_change_{year}.tif", tr_window_change, profile)
        write_geotiff(OUTPUT_DIR / f"TR_fixed_window_{year}.tif", tr_fixed_window, profile)
        write_geotiff(OUTPUT_DIR / f"TR_sos_change_{year}.tif", tr_sos_change, profile)
        write_geotiff(OUTPUT_DIR / f"TR_pos_change_{year}.tif", tr_pos_change, profile)

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
