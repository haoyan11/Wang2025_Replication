#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module 02: TRc（累积蒸腾）计算
计算SOS-POS窗口内的累积蒸腾量
TRc_y = Σ[SOS to POS] TR(t)
"""

import numpy as np
import rasterio
from pathlib import Path
from tqdm import tqdm
from datetime import datetime, timedelta
from rasterio.windows import Window
from concurrent.futures import ProcessPoolExecutor, as_completed
import warnings
warnings.filterwarnings('ignore')

# ==================== 全局配置 ====================
ROOT = Path(r"I:\F\Data4")
PHENO_DIR = ROOT / "Phenology_Output_1" / "SIF_phenology"
TR_DAILY_DIR = ROOT / "Meteorological Data" / "GLEAM" / "Daily_0p5deg_TIF" / "Et"
OUTPUT_DIR = ROOT / "Wang2025_Analysis" / "TRc_annual"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

YEAR_START = 1982
YEAR_END = 2018
BLOCK_SIZE = 128
MAX_WORKERS = 2
NODATA_OUT = -9999.0

# ==================== 辅助函数 ====================
def get_TR_file_path(date_obj):
    """
    获取ERA5-Land TR日数据文件路径
    格式: ERA5L_ET_transp_Daily_mm_YYYYMMDD.tif
    """
    yyyymmdd = date_obj.strftime("%Y%m%d")

    # ERA5-Land格式
    file_path = TR_DAILY_DIR / f"ERA5L_ET_transp_Daily_mm_{yyyymmdd}.tif"

    if file_path.exists():
        return file_path

    return None

def is_leap_year(year):
    """判断是否为闰年"""
    return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)

# ==================== 核心计算函数 ====================
def calculate_TRc_for_year(year, mask):
    """
    计算单个年份的TRc

    Parameters:
    -----------
    year : int
        年份
    mask : ndarray
        有效掩膜（布尔数组）

    Returns:
    --------
    TRc : ndarray
        累积蒸腾栅格
    profile : dict
        栅格配置
    """
    print(f"\n=== 计算TRc: {year} ===")

    # 读取物候数据（物候代码输出格式：小写sos_gpp）
    sos_file = PHENO_DIR / f"sos_gpp_{year}.tif"
    pos_file = PHENO_DIR / f"pos_doy_gpp_{year}.tif"  # 注意是pos_doy

    if not sos_file.exists() or not pos_file.exists():
        print(f"  ✗ 错误：缺少物候数据")
        print(f"    SOS文件: {sos_file}")
        print(f"    SOS存在: {sos_file.exists()}")
        print(f"    POS文件: {pos_file}")
        print(f"    POS存在: {pos_file.exists()}")
        print("  请先运行物候计算代码或确认物候数据路径")
        return None, None

    with rasterio.open(sos_file) as src:
        sos = src.read(1)
        profile = src.profile.copy()
        height, width = src.height, src.width

    with rasterio.open(pos_file) as src:
        pos = src.read(1)

    # 初始化TRc数组
    TRc = np.zeros((height, width), dtype=np.float32)
    TRc[~mask] = NODATA_OUT

    # 生成日期列表
    days_in_year = 366 if is_leap_year(year) else 365
    dates_year = [datetime(year, 1, 1) + timedelta(days=i) for i in range(days_in_year)]

    print(f"  年份: {year}, 天数: {days_in_year}")

    # 逐日累加TR
    missing_files = 0
    valid_days = 0

    for date_obj in tqdm(dates_year, desc="累加TR", leave=False):
        doy = date_obj.timetuple().tm_yday
        tr_file = get_TR_file_path(date_obj)

        if tr_file is None or not tr_file.exists():
            missing_files += 1
            continue

        try:
            with rasterio.open(tr_file) as src:
                tr_daily = src.read(1)

            # 累加：仅在SOS ≤ doy ≤ POS的像元
            for i in range(height):
                for j in range(width):
                    if not mask[i, j]:
                        continue

                    sos_ij = sos[i, j]
                    pos_ij = pos[i, j]

                    # 检查有效性
                    if sos_ij == NODATA_OUT or pos_ij == NODATA_OUT:
                        continue

                    if sos_ij <= doy <= pos_ij:
                        tr_val = tr_daily[i, j]
                        if tr_val != NODATA_OUT and not np.isnan(tr_val):
                            TRc[i, j] += tr_val

            valid_days += 1

        except Exception as e:
            tqdm.write(f"  ⚠ 读取失败: {tr_file.name} - {str(e)}")
            missing_files += 1
            continue

    print(f"  有效天数: {valid_days}/{days_in_year}, 缺失: {missing_files}")

    # 质量检查
    TRc_valid = TRc[mask & (TRc != NODATA_OUT)]
    if len(TRc_valid) > 0:
        print(f"  TRc范围: {np.min(TRc_valid):.2f} - {np.max(TRc_valid):.2f} mm")
        print(f"  TRc平均: {np.mean(TRc_valid):.2f} mm")
        print(f"  有效像元: {len(TRc_valid)}")
    else:
        print(f"  ⚠ 警告：没有有效的TRc值")

    return TRc, profile

# ==================== 块处理版本（优化内存） ====================
def calculate_TRc_block_optimized(year, mask):
    """
    块处理版本（更节省内存）

    使用块处理避免一次性读取所有日数据到内存
    """
    print(f"\n=== 计算TRc（块处理）: {year} ===")

    # 读取物候数据（物候代码输出格式）
    sos_file = PHENO_DIR / f"sos_gpp_{year}.tif"
    pos_file = PHENO_DIR / f"pos_doy_gpp_{year}.tif"

    if not sos_file.exists() or not pos_file.exists():
        print(f"  ✗ 错误：缺少物候数据")
        print(f"    SOS: {sos_file} ({'存在' if sos_file.exists() else '不存在'})")
        print(f"    POS: {pos_file} ({'存在' if pos_file.exists() else '不存在'})")
        return None, None

    with rasterio.open(sos_file) as src:
        sos = src.read(1)
        profile = src.profile.copy()
        height, width = src.height, src.width

    with rasterio.open(pos_file) as src:
        pos = src.read(1)

    # 初始化TRc
    TRc = np.full((height, width), NODATA_OUT, dtype=np.float32)

    # 生成块窗口
    block_windows = [
        Window(col_off, row_off,
               min(BLOCK_SIZE, width - col_off),
               min(BLOCK_SIZE, height - row_off))
        for row_off in range(0, height, BLOCK_SIZE)
        for col_off in range(0, width, BLOCK_SIZE)
    ]

    print(f"  总块数: {len(block_windows)}")

    # 日期列表
    days_in_year = 366 if is_leap_year(year) else 365
    dates_year = [datetime(year, 1, 1) + timedelta(days=i) for i in range(days_in_year)]

    # 逐块处理
    for win in tqdm(block_windows, desc="处理块"):
        # 提取块内物候
        sos_block = sos[win.row_off:win.row_off + win.height,
                        win.col_off:win.col_off + win.width]
        pos_block = pos[win.row_off:win.row_off + win.height,
                        win.col_off:win.col_off + win.width]
        mask_block = mask[win.row_off:win.row_off + win.height,
                          win.col_off:win.col_off + win.width]

        # 初始化块内TRc
        TRc_block = np.zeros((win.height, win.width), dtype=np.float32)

        # 逐日累加
        for date_obj in dates_year:
            doy = date_obj.timetuple().tm_yday
            tr_file = get_TR_file_path(date_obj)

            if tr_file is None or not tr_file.exists():
                continue

            try:
                with rasterio.open(tr_file) as src:
                    tr_block = src.read(1, window=win)

                # 累加条件：SOS ≤ doy ≤ POS 且 mask == True
                in_window = (sos_block <= doy) & (doy <= pos_block) & mask_block
                valid_tr = (tr_block != NODATA_OUT) & ~np.isnan(tr_block)
                accumulate_mask = in_window & valid_tr

                TRc_block[accumulate_mask] += tr_block[accumulate_mask]

            except Exception:
                continue

        # 将块写入完整数组
        TRc_block[~mask_block] = NODATA_OUT
        TRc[win.row_off:win.row_off + win.height,
            win.col_off:win.col_off + win.width] = TRc_block

    # 质量检查
    TRc_valid = TRc[mask & (TRc != NODATA_OUT)]
    if len(TRc_valid) > 0:
        print(f"  TRc范围: {np.min(TRc_valid):.2f} - {np.max(TRc_valid):.2f} mm")
        print(f"  TRc平均: {np.mean(TRc_valid):.2f} mm")
        print(f"  有效像元: {len(TRc_valid)}")

    return TRc, profile

# ==================== 主处理流程 ====================
def main(use_block_processing=True):
    """
    主处理流程

    Parameters:
    -----------
    use_block_processing : bool
        是否使用块处理（推荐True以节省内存）
    """
    print("\n" + "="*70)
    print("Module 02: TRc计算")
    print("="*70)

    years = list(range(YEAR_START, YEAR_END + 1))

    # 读取掩膜
    print("\n[1/3] 读取掩膜...")
    mask_file = ROOT / "Wang2025_Analysis" / "masks" / "combined_mask.tif"

    if not mask_file.exists():
        print(f"  ✗ 错误：掩膜文件不存在: {mask_file}")
        print("  请先运行 Module 00: 00_data_preparation.py")
        return

    with rasterio.open(mask_file) as src:
        mask = src.read(1).astype(bool)

    print(f"  ✓ 有效像元数: {np.sum(mask)}")

    # 检查TR数据
    print("\n[2/3] 检查TR数据...")
    test_date = datetime(2000, 1, 15)
    test_file = get_TR_file_path(test_date)

    if test_file is None or not test_file.exists():
        print(f"  ✗ 错误：找不到TR数据")
        print(f"    测试日期: {test_date.strftime('%Y-%m-%d')}")
        print(f"    预期路径示例: {TR_DAILY_DIR / 'Et_20000115.tif'}")
        print(f"\n请检查以下路径是否正确:")
        print(f"  TR_DAILY_DIR = {TR_DAILY_DIR}")
        print(f"\n或修改 get_TR_file_path() 函数以匹配您的文件命名格式")
        return
    else:
        print(f"  ✓ 找到TR数据: {test_file.name}")

    # 逐年计算TRc
    print(f"\n[3/3] 逐年计算TRc ({YEAR_START}-{YEAR_END})...")

    for year in years:
        # 选择处理方法
        if use_block_processing:
            TRc, profile = calculate_TRc_block_optimized(year, mask)
        else:
            TRc, profile = calculate_TRc_for_year(year, mask)

        if TRc is None:
            print(f"  ✗ 跳过年份 {year}")
            continue

        # 保存结果
        output_file = OUTPUT_DIR / f"TRc_{year}.tif"
        profile.update(dtype=rasterio.float32, count=1, compress='lzw', nodata=NODATA_OUT)

        with rasterio.open(output_file, 'w', **profile) as dst:
            dst.write(TRc.astype(np.float32), 1)

        print(f"  ✓ 已保存: {output_file.name}")

    print("\n" + "="*70)
    print("✓ TRc计算完成！")
    print(f"输出目录: {OUTPUT_DIR}")
    print("="*70)

# ==================== 主程序 ====================
if __name__ == "__main__":
    # 使用块处理（节省内存）
    main(use_block_processing=True)

    # 如果内存充足，可使用常规方法（更快）
    # main(use_block_processing=False)
