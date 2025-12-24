#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module 03: 原版TRc分解方法（完全按照Wang 2025论文）

按照Wang (2025)原文方法（真实逐日累加版）：

  1. TRc_y = Σ[SOS_y to POS_y] TR_y(t)  # 当年累积蒸腾

  2. TRc_av = Σ[SOSav to POSav] TR_daily_av(t)  # 多年平均基线 ⭐关键
     使用多年平均日TR在固定窗口累加

  3. TRpheno_y = 物候效应（SOS偏移导致的变化）
     如果 SOS_y < SOSav: TRpheno_y = Σ[SOS_y to SOSav] TR_daily_av(t)  # 提前，增益
     如果 SOS_y > SOSav: TRpheno_y = -Σ[SOSav to SOS_y] TR_daily_av(t) # 推迟，损失
     如果 SOS_y == SOSav: TRpheno_y = 0

  4. TRproduct_y = TRc_y - TRc_av - TRpheno_y  # 剩余项（生产力效应）

重要说明：
  - TR_daily_av: 365天的多年平均日TR气候态
  - SOSav, POSav: 多年平均SOS和POS
  - TRc_av 必须用 TR_daily_av 在固定窗口[SOSav, POSav]累加
  - 不能用 mean(TRc_y)，因为每年窗口不同
"""

import numpy as np
import rasterio
from pathlib import Path
from tqdm import tqdm

# ==================== 全局配置 ====================
ROOT = Path(r"I:\F\Data4")
PHENO_DIR = ROOT / "Phenology_Output_1" / "T_phenology"  # T物候
TRC_DIR = ROOT / "Wang2025_Analysis" / "TRc_annual_T"
CLIMATOLOGY_DIR = ROOT / "Wang2025_Analysis" / "Climatology_T"  # 气候态数据（T物候版）
OUTPUT_DIR = ROOT / "Wang2025_Analysis" / "Decomposition_T"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

YEAR_START = 1982
YEAR_END = 2018
NODATA_OUT = -9999.0  # 仅用于输出，读取时使用真实nodata

# ==================== 辅助函数 ====================
def _is_valid_value(value, nodata):
    """
    判断值是否有效（非NODATA）

    兼容各种NODATA类型：-9999, NaN, 3.4e38等
    """
    if nodata is None:
        return np.isfinite(value)
    elif np.isnan(nodata):
        return np.isfinite(value)
    else:
        return (value != nodata) & np.isfinite(value)

def read_geotiff(file_path):
    """
    读取单波段GeoTIFF

    Returns:
    --------
    data : ndarray
        栅格数据
    profile : dict
        栅格配置
    nodata : float or None
        真实nodata值
    """
    with rasterio.open(file_path) as src:
        data = src.read(1)
        profile = src.profile.copy()
        nodata = src.nodata
    return data, profile, nodata

def write_geotiff(file_path, data, profile):
    """写入单波段GeoTIFF"""
    profile.update(dtype=rasterio.float32, count=1, compress='lzw', nodata=NODATA_OUT)
    with rasterio.open(file_path, 'w', **profile) as dst:
        dst.write(data.astype(np.float32), 1)

def _sum_cum_range(cum_arr, start_doy, end_doy):
    """
    使用累积和快速计算区间和（1-based DOY，区间为 [start, end]）
    start_doy/end_doy 为二维数组，需满足 1 <= start <= end <= 365
    """
    flat = cum_arr.reshape(cum_arr.shape[0], -1)
    start = start_doy.ravel()
    end = end_doy.ravel()
    idx = np.arange(start.size)

    end_vals = flat[end - 1, idx]
    start_vals = np.where(start > 1, flat[start - 2, idx], 0.0)
    return (end_vals - start_vals).reshape(start_doy.shape)

# ==================== 核心函数 ====================
def load_climatology_data():
    """
    加载气候态数据（TR日气候态、SOSav、POSav）

    Returns:
    --------
    TR_daily_av : ndarray, shape (365, H, W)
        365天的多年平均日TR
    SOSav : ndarray, shape (H, W)
        多年平均SOS
    POSav : ndarray, shape (H, W)
        多年平均POS
    nodata_sos : float or None
        SOSav的真实nodata值
    nodata_pos : float or None
        POSav的真实nodata值
    profile : dict
        栅格配置
    """
    print("\n=== 加载气候态数据 ===")

    # 加载TR日气候态（365波段）
    tr_clim_file = CLIMATOLOGY_DIR / "TR_daily_climatology.tif"
    if not tr_clim_file.exists():
        raise FileNotFoundError(
            f"找不到TR日气候态文件: {tr_clim_file}\n"
            f"请先运行 02_TRc_calculation.py 中的 save_climatology_data() 函数"
        )

    with rasterio.open(tr_clim_file) as src:
        # 读取所有365个波段
        TR_daily_av = src.read().astype(np.float32)  # shape: (365, H, W)
        profile = src.profile.copy()
        nodata_tr = src.nodata
        if nodata_tr is None:
            nodata_tr = NODATA_OUT
        if np.isnan(nodata_tr):
            TR_daily_av[~np.isfinite(TR_daily_av)] = np.nan
        else:
            TR_daily_av[TR_daily_av == nodata_tr] = np.nan
        TR_daily_av[TR_daily_av < 0] = np.nan
        print(f"  ✓ TR日气候态: {tr_clim_file.name}, 形状={TR_daily_av.shape}")

    # 加载SOSav
    sos_av_file = CLIMATOLOGY_DIR / "SOSav.tif"
    if not sos_av_file.exists():
        raise FileNotFoundError(
            f"找不到SOSav文件: {sos_av_file}\n"
            f"请先运行 02_TRc_calculation.py 中的 save_climatology_data() 函数"
        )

    with rasterio.open(sos_av_file) as src:
        SOSav = src.read(1)
        nodata_sos = src.nodata
        valid_count = np.sum(_is_valid_value(SOSav, nodata_sos))
        print(f"  ✓ SOSav: {sos_av_file.name}, nodata={nodata_sos}, 有效像元={valid_count}")

    # 加载POSav
    pos_av_file = CLIMATOLOGY_DIR / "POSav.tif"
    if not pos_av_file.exists():
        raise FileNotFoundError(
            f"找不到POSav文件: {pos_av_file}\n"
            f"请先运行 02_TRc_calculation.py 中的 save_climatology_data() 函数"
        )

    with rasterio.open(pos_av_file) as src:
        POSav = src.read(1)
        nodata_pos = src.nodata
        valid_count = np.sum(_is_valid_value(POSav, nodata_pos))
        print(f"  ✓ POSav: {pos_av_file.name}, nodata={nodata_pos}, 有效像元={valid_count}")

    return TR_daily_av, SOSav, POSav, nodata_sos, nodata_pos, profile

def calculate_TRc_av_from_climatology(TR_cum, SOSav, POSav, nodata_sos, nodata_pos, mask, tr_valid):
    """
    按照Wang (2025)论文真实定义计算TRc_av

    TRc_av = Σ[t=SOSav to POSav] TR_daily_av(t)

    使用多年平均日TR在固定窗口[SOSav, POSav]内累加

    Parameters:
    -----------
    TR_cum : ndarray, shape (365, H, W)
        365天的多年平均日TR的累积和
    SOSav : ndarray, shape (H, W)
        多年平均SOS
    POSav : ndarray, shape (H, W)
        多年平均POS
    nodata_sos : float or None
        SOSav的真实nodata值
    nodata_pos : float or None
        POSav的真实nodata值
    mask : ndarray
        有效掩膜

    Returns:
    --------
    TRc_av : ndarray, shape (H, W)
        多年平均TRc（按论文定义）
    """
    print(f"\n=== 计算 TRc_av（论文定义：Σ[SOSav:POSav] TR_daily_av）===")

    height, width = SOSav.shape

    # 构建有效掩膜（使用真实nodata值）
    valid = (
        _is_valid_value(SOSav, nodata_sos) &
        _is_valid_value(POSav, nodata_pos) &
        (SOSav > 0) & (SOSav <= 366) &
        (POSav > 0) & (POSav <= 366) &
        (POSav >= SOSav) &
        mask &
        tr_valid
    )

    # 初始化输出
    TRc_av = np.full((height, width), NODATA_OUT, dtype=np.float32)

    print(f"  向量化累加TR_daily_av[SOSav:POSav]...")

    sos_av_raw = np.rint(SOSav).astype(np.int32)
    pos_av_raw = np.rint(POSav).astype(np.int32)
    sos_av = np.clip(sos_av_raw, 1, 365)
    pos_av = np.clip(pos_av_raw, 1, 365)

    sum_all = _sum_cum_range(TR_cum, sos_av, pos_av)
    TRc_av[valid] = sum_all[valid]

    # 统计（使用valid掩膜，因为输出用NODATA_OUT）
    TRc_av_valid = TRc_av[valid]
    print(f"  ✓ TRc_av统计: 平均={TRc_av_valid.mean():.2f} mm, "
          f"范围=[{TRc_av_valid.min():.2f}, {TRc_av_valid.max():.2f}], "
          f"有效像元={TRc_av_valid.size}")

    return TRc_av

def decompose_TRc_original(year, TRc_av, TR_cum, SOSav, nodata_sos, mask, tr_valid):
    """
    原版分解方法（Wang 2025）- 真实逐日累加版本

    Parameters:
    -----------
    year : int
        年份
    TRc_av : ndarray, shape (H, W)
        多年平均TRc
    TR_cum : ndarray, shape (365, H, W)
        365天的多年平均日TR累积和
    SOSav : ndarray, shape (H, W)
        多年平均SOS
    nodata_sos : float or None
        SOS的真实nodata值
    mask : ndarray
        有效掩膜

    Returns:
    --------
    TRpheno, TRproduct : ndarray
        分解后的两个分量
    """
    # 读取当年数据（使用更新后的read_geotiff，返回3个值）
    TRc_y, _, nodata_trc = read_geotiff(TRC_DIR / f"TRc_{year}.tif")
    sos_y, _, nodata_sos_y = read_geotiff(PHENO_DIR / f"sos_t_{year}.tif")

    height, width = TRc_y.shape

    # 构建有效掩膜（使用真实nodata值）
    valid = (
        _is_valid_value(TRc_y, nodata_trc) &
        (TRc_av != NODATA_OUT) &
        _is_valid_value(sos_y, nodata_sos_y) &
        _is_valid_value(SOSav, nodata_sos) &
        (sos_y > 0) & (sos_y <= 366) &
        (SOSav > 0) & (SOSav <= 366) &
        mask &
        tr_valid
    )

    # 初始化输出
    TRpheno = np.full((height, width), NODATA_OUT, dtype=np.float32)
    TRproduct = np.full((height, width), NODATA_OUT, dtype=np.float32)

    # 原版公式（Wang 2025真实方法）：
    # 如果 SOS_y < SOSav: TRpheno = Σ[SOS_y to SOSav] TR_daily_av(t)
    # 如果 SOS_y > SOSav: TRpheno = -Σ[SOSav to SOS_y] TR_daily_av(t)
    # 如果 SOS_y == SOSav: TRpheno = 0
    # TRproduct = TRc_y - TRc_av - TRpheno

    print(f"  处理年份 {year}: 向量化计算TRpheno...")

    sos_y_raw = np.rint(sos_y).astype(np.int32)
    sos_av_raw = np.rint(SOSav).astype(np.int32)

    sign = np.where(
        sos_y_raw < sos_av_raw, 1.0,
        np.where(sos_y_raw > sos_av_raw, -1.0, 0.0)
    )

    start_raw = np.where(sos_y_raw < sos_av_raw, sos_y_raw, sos_av_raw)
    end_raw = np.where(sos_y_raw < sos_av_raw, sos_av_raw - 1, sos_y_raw - 1)

    valid_range = (
        (start_raw >= 1) & (end_raw >= start_raw) &
        (start_raw <= 366) & (end_raw <= 366)
    )

    start = np.clip(start_raw, 1, 365)
    end = np.clip(end_raw, 1, 365)

    sum_range = _sum_cum_range(TR_cum, start, end)
    sum_range[~valid_range] = 0.0

    TRpheno[valid] = (sign * sum_range)[valid]
    TRproduct[valid] = TRc_y[valid] - TRc_av[valid] - TRpheno[valid]

    return TRpheno, TRproduct

def process_all_years():
    """主处理流程"""
    print("\n" + "="*70)
    print("Wang (2025) 原版TRc分解（真实逐日累加版）")
    print("="*70)

    years = list(range(YEAR_START, YEAR_END + 1))

    # Step 1: 读取掩膜
    print("\n[1/5] 读取掩膜...")
    mask, profile, _ = read_geotiff(OUTPUT_DIR.parent / "masks" / "combined_mask.tif")
    mask = mask.astype(bool)
    print(f"  有效像元数: {np.sum(mask)}")

    # Step 2: 加载气候态数据
    print("\n[2/5] 加载气候态数据...")
    TR_daily_av, SOSav, POSav, nodata_sos, nodata_pos, profile = load_climatology_data()

    tr_valid = np.any(np.isfinite(TR_daily_av), axis=0)
    TR_cum = np.nancumsum(TR_daily_av, axis=0).astype(np.float32)

    # Step 3: 计算多年平均TRc（按论文真实定义）
    print("\n[3/5] 计算多年平均TRc（按论文定义）...")
    TRc_av = calculate_TRc_av_from_climatology(TR_cum, SOSav, POSav, nodata_sos, nodata_pos, mask, tr_valid)

    # 保存多年平均值
    write_geotiff(OUTPUT_DIR / "TRc_av.tif", TRc_av, profile)
    print(f"  ✓ 已保存: TRc_av.tif")

    # Step 4: 逐年分解
    print(f"\n[4/5] 逐年分解 ({YEAR_START}-{YEAR_END})...")
    for year in tqdm(years, desc="分解进度"):
        TRpheno, TRproduct = decompose_TRc_original(year, TRc_av, TR_cum, SOSav, nodata_sos, mask, tr_valid)

        # 保存结果
        write_geotiff(OUTPUT_DIR / f"TRpheno_{year}.tif", TRpheno, profile)
        write_geotiff(OUTPUT_DIR / f"TRproduct_{year}.tif", TRproduct, profile)

    # Step 5: 质量检查
    print("\n[5/5] 质量检查...")
    test_year = 2000
    TRc_test, _, nodata_trc_test = read_geotiff(TRC_DIR / f"TRc_{test_year}.tif")
    TRpheno_test, _, _ = read_geotiff(OUTPUT_DIR / f"TRpheno_{test_year}.tif")
    TRproduct_test, _, _ = read_geotiff(OUTPUT_DIR / f"TRproduct_{test_year}.tif")

    # 构建有效掩膜（使用真实nodata值）
    valid = (
        mask &
        _is_valid_value(TRc_test, nodata_trc_test) &
        (TRpheno_test != NODATA_OUT) &  # 输出文件用NODATA_OUT
        (TRproduct_test != NODATA_OUT)
    )
    residual = TRc_test[valid] - TRc_av[valid] - TRpheno_test[valid] - TRproduct_test[valid]

    print(f"\n  检验年份: {test_year}")
    print(f"  残差统计 (TRc - TRc_av - TRpheno - TRproduct):")
    print(f"    Mean: {np.mean(residual):.6f} mm")
    print(f"    Std:  {np.std(residual):.6f} mm")
    print(f"    Max:  {np.max(np.abs(residual)):.6f} mm")

    if np.max(np.abs(residual)) < 1e-3:
        print("  ✓ 质量检查通过！残差可忽略。")
    else:
        print("  ⚠ 警告：残差较大，请检查计算过程。")

    # 分解效果统计
    print(f"\n  分解效果统计 ({test_year}):")
    print(f"    TRc平均:      {np.mean(TRc_test[valid]):.2f} mm")
    print(f"    TRc_av平均:   {np.mean(TRc_av[valid]):.2f} mm")
    print(f"    TRpheno平均:  {np.mean(TRpheno_test[valid]):.2f} mm")
    print(f"    TRproduct平均:{np.mean(TRproduct_test[valid]):.2f} mm")

    print("\n" + "="*70)
    print("✓ 原版分解完成（真实逐日累加方法）！")
    print(f"输出目录: {OUTPUT_DIR}")
    print("\n方法说明：")
    print("  - 此版本使用真实逐日累加TR气候态计算TRpheno")
    print("  - 替代了简化的线性比例分配方法")
    print("  - 更准确反映物候变化对蒸腾的影响")
    print("="*70)

# ==================== 主程序 ====================
if __name__ == "__main__":
    process_all_years()
