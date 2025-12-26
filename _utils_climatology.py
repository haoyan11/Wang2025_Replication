#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
工具脚本: 计算日均气候平均态
为端点法分解提供TR_clim数据
"""

import numpy as np
import rasterio
from pathlib import Path
from tqdm import tqdm
from datetime import datetime, timedelta
from concurrent.futures import ProcessPoolExecutor
import warnings
warnings.filterwarnings('ignore')

# ==================== 全局配置 ====================
ROOT = Path(r"I:\F\Data4")
TR_DAILY_DIR = ROOT / "Meteorological Data" / "GLEAM" / "Daily_0p5deg_TIF" / "Et"
OUTPUT_DIR = ROOT / "Wang2025_Analysis" / "climatology"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

YEAR_START = 1982
YEAR_END = 2018
NODATA_OUT = -9999.0

# ==================== 辅助函数 ====================
def get_file_path_with_date(date_obj):
    """获取GLEAM TR日数据文件路径"""
    # 尝试多种命名格式
    yyyymmdd = date_obj.strftime("%Y%m%d")
    year = date_obj.year
    doy = date_obj.timetuple().tm_yday

    # 格式1: Et_YYYYMMDD.tif
    file1 = TR_DAILY_DIR / f"Et_{yyyymmdd}.tif"
    if file1.exists():
        return file1

    # 格式2: Et_YYYY_DOY.tif
    file2 = TR_DAILY_DIR / f"Et_{year}_{doy:03d}.tif"
    if file2.exists():
        return file2

    # 格式3: YYYY/Et_YYYYMMDD.tif
    file3 = TR_DAILY_DIR / str(year) / f"Et_{yyyymmdd}.tif"
    if file3.exists():
        return file3

    return None

def is_leap_year(year):
    """判断是否为闰年"""
    return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)

# ==================== 核心函数 ====================
def calculate_daily_climatology():
    """
    计算366天的日均气候平均态

    输出: TR_daily_climatology.tif (366, H, W)
          每层对应DOY 1-366的气候平均值
    """
    print("\n" + "="*70)
    print("计算日均TR气候平均态")
    print("="*70)

    years = list(range(YEAR_START, YEAR_END + 1))

    # 读取模板获取栅格尺寸
    print("\n[1/3] 初始化...")
    test_date = datetime(2000, 1, 1)
    test_file = get_file_path_with_date(test_date)

    if test_file is None or not test_file.exists():
        print(f"  ✗ 错误：找不到测试文件（{test_date.strftime('%Y-%m-%d')}）")
        print(f"  请检查TR数据路径: {TR_DAILY_DIR}")
        return

    with rasterio.open(test_file) as src:
        profile = src.profile.copy()
        height, width = src.height, src.width

    print(f"  栅格尺寸: {height} x {width}")
    print(f"  年份范围: {YEAR_START}-{YEAR_END}")

    # 初始化累加器（366天 × 高度 × 宽度）
    # 使用366天以包含闰年2月29日
    doy_sum = np.zeros((366, height, width), dtype=np.float64)
    doy_count = np.zeros((366, height, width), dtype=np.int32)

    # 逐年读取并累加
    print("\n[2/3] 逐年读取并累加...")
    for year in tqdm(years, desc="处理年份"):
        days_in_year = 366 if is_leap_year(year) else 365

        for doy in range(1, days_in_year + 1):
            date_obj = datetime(year, 1, 1) + timedelta(days=doy - 1)
            file_path = get_file_path_with_date(date_obj)

            if file_path is None or not file_path.exists():
                continue

            try:
                with rasterio.open(file_path) as src:
                    tr_data = src.read(1)

                # 将有效值累加到对应DOY
                valid = tr_data != NODATA_OUT
                doy_sum[doy - 1][valid] += tr_data[valid]
                doy_count[doy - 1][valid] += 1

            except Exception as e:
                tqdm.write(f"  ⚠ 读取失败: {file_path.name} - {str(e)}")
                continue

    # 计算平均值
    print("\n[3/3] 计算气候平均值...")
    climatology = np.full((366, height, width), NODATA_OUT, dtype=np.float32)

    for doy in tqdm(range(366), desc="计算平均"):
        valid = doy_count[doy] > 0
        climatology[doy][valid] = doy_sum[doy][valid] / doy_count[doy][valid]

    # 对于非闰年没有的2月29日（DOY 60），使用前后平均
    doy_feb28 = 59  # DOY 60 (2月29日，索引59)
    doy_feb29_mask = doy_count[doy_feb28] == 0
    if np.sum(doy_feb29_mask) > 0:
        print(f"  对2月29日（DOY 60）使用前后日期插值...")
        climatology[doy_feb28] = (climatology[doy_feb28 - 1] + climatology[doy_feb28 + 1]) / 2

    # 保存结果
    output_file = OUTPUT_DIR / "TR_daily_climatology.tif"
    profile.update(
        count=366,
        dtype=rasterio.float32,
        compress='lzw',
        nodata=NODATA_OUT
    )

    print(f"\n保存文件: {output_file}")
    with rasterio.open(output_file, 'w', **profile) as dst:
        for doy in range(366):
            dst.write(climatology[doy], doy + 1)

    # 质量检查
    print("\n质量检查:")
    valid_pixels = np.sum(climatology[0] != NODATA_OUT)
    print(f"  有效像元数: {valid_pixels}")
    print(f"  年均TR范围: {np.min(climatology[climatology != NODATA_OUT]):.2f} - "
          f"{np.max(climatology[climatology != NODATA_OUT]):.2f} mm/day")

    # 检查各月代表性
    monthly_coverage = []
    for month in range(1, 13):
        # 每月选择中旬代表日
        mid_day = datetime(2000, month, 15).timetuple().tm_yday - 1
        coverage = np.sum(doy_count[mid_day] > 0) / (height * width) * 100
        monthly_coverage.append(coverage)
        print(f"  {month:2d}月数据覆盖率: {coverage:.1f}%")

    print("\n" + "="*70)
    print("✓ 日均气候平均态计算完成！")
    print(f"输出文件: {output_file}")
    print("="*70)

# ==================== 主程序 ====================
if __name__ == "__main__":
    calculate_daily_climatology()
