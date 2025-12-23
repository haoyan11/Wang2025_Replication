#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module 01: 物候提取
从SIF日尺度数据提取 SOS（春季物候开始）、POS（生长季峰值）、EOS（秋季物候结束）
使用阈值法（20%阈值）
"""

import numpy as np
import rasterio
from pathlib import Path
from tqdm import tqdm
from datetime import datetime, timedelta
from scipy.ndimage import uniform_filter1d
import warnings
warnings.filterwarnings('ignore')

# ==================== 全局配置 ====================
ROOT = Path(r"I:\F\Data4")
SIF_DAILY_DIR = ROOT / "SIF_Data" / "CSIF_daily"
OUTPUT_DIR = ROOT / "Phenology_Output_1" / "SIF_phenology"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

YEAR_START = 1982
YEAR_END = 2018
NODATA_OUT = -9999.0

# 物候提取参数
PHENO_THRESHOLD = 0.20  # 20%阈值
SMOOTH_WINDOW = 7       # 7天移动平均
SOS_MAX_DOY = 150       # SOS最大约束（5月30日）
MIN_VALID_DAYS = 200    # 最少有效天数

# ==================== 辅助函数 ====================
def get_SIF_file_path(date_obj):
    """
    获取SIF日数据文件路径
    支持多种命名格式
    """
    yyyymmdd = date_obj.strftime("%Y%m%d")
    year = date_obj.year
    doy = date_obj.timetuple().tm_yday

    # 格式1: SIF_YYYYMMDD.tif
    file1 = SIF_DAILY_DIR / f"SIF_{yyyymmdd}.tif"
    if file1.exists():
        return file1

    # 格式2: SIF_YYYY_DOY.tif
    file2 = SIF_DAILY_DIR / f"SIF_{year}_{doy:03d}.tif"
    if file2.exists():
        return file2

    # 格式3: YYYY/SIF_YYYYMMDD.tif
    file3 = SIF_DAILY_DIR / str(year) / f"SIF_{yyyymmdd}.tif"
    if file3.exists():
        return file3

    # 格式4: CSIF_YYYYMMDD.tif
    file4 = SIF_DAILY_DIR / f"CSIF_{yyyymmdd}.tif"
    if file4.exists():
        return file4

    return None

def is_leap_year(year):
    """判断是否为闰年"""
    return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)

# ==================== 物候提取算法 ====================
def extract_phenology_threshold(sif_timeseries, smooth=True):
    """
    阈值法提取物候

    Parameters:
    -----------
    sif_timeseries : ndarray
        (365 or 366,) 日尺度SIF时间序列
    smooth : bool
        是否平滑

    Returns:
    --------
    sos, pos, eos : int
        春季开始、峰值、秋季结束的DOY
    """
    n_days = len(sif_timeseries)

    # 检查有效数据
    valid_mask = ~np.isnan(sif_timeseries)
    if np.sum(valid_mask) < MIN_VALID_DAYS:
        return NODATA_OUT, NODATA_OUT, NODATA_OUT

    # 填充缺失值（线性插值）
    sif = sif_timeseries.copy()
    if np.any(~valid_mask):
        x_valid = np.where(valid_mask)[0]
        y_valid = sif[valid_mask]
        x_invalid = np.where(~valid_mask)[0]
        if len(x_valid) > 1:
            sif[x_invalid] = np.interp(x_invalid, x_valid, y_valid)

    # 平滑处理
    if smooth:
        sif_smooth = uniform_filter1d(sif, size=SMOOTH_WINDOW, mode='nearest')
    else:
        sif_smooth = sif

    # 找到峰值（POS）
    pos_idx = np.argmax(sif_smooth)
    pos = pos_idx + 1  # DOY从1开始

    # 计算阈值
    sif_min = np.min(sif_smooth)
    sif_max = sif_smooth[pos_idx]
    threshold = sif_min + PHENO_THRESHOLD * (sif_max - sif_min)

    # 找SOS：POS之前第一个超过阈值的DOY
    sos = NODATA_OUT
    for i in range(pos_idx, -1, -1):
        if sif_smooth[i] <= threshold:
            sos = i + 2  # +2是因为我们要的是超过阈值的第一天
            break

    # 约束SOS ≤ 150 DOY
    if sos != NODATA_OUT and sos > SOS_MAX_DOY:
        sos = SOS_MAX_DOY

    # 找EOS：POS之后第一个降到阈值以下的DOY
    eos = NODATA_OUT
    for i in range(pos_idx, n_days):
        if sif_smooth[i] <= threshold:
            eos = i + 1
            break

    # 如果没找到EOS，使用年末
    if eos == NODATA_OUT:
        eos = n_days

    # 质量检查
    if sos != NODATA_OUT and eos != NODATA_OUT:
        if sos >= pos or pos >= eos or sos >= eos:
            return NODATA_OUT, NODATA_OUT, NODATA_OUT

    return sos, pos, eos

# ==================== 逐年处理 ====================
def process_year(year, mask):
    """
    处理单个年份的物候提取

    Parameters:
    -----------
    year : int
        年份
    mask : ndarray
        有效掩膜（布尔数组）

    Returns:
    --------
    sos, pos, eos : ndarray
        物候栅格数据
    """
    print(f"\n=== 处理年份: {year} ===")

    days_in_year = 366 if is_leap_year(year) else 365
    dates = [datetime(year, 1, 1) + timedelta(days=i) for i in range(days_in_year)]

    # 读取模板获取栅格尺寸
    test_file = get_SIF_file_path(dates[0])
    if test_file is None or not test_file.exists():
        print(f"  ✗ 错误：找不到{year}年的SIF数据")
        print(f"    预期路径示例: {SIF_DAILY_DIR / f'SIF_{year}0101.tif'}")
        return None, None, None

    with rasterio.open(test_file) as src:
        profile = src.profile.copy()
        height, width = src.height, src.width

    # 初始化SIF时间序列数组
    print(f"  读取{days_in_year}天SIF数据...")
    sif_stack = np.full((days_in_year, height, width), np.nan, dtype=np.float32)

    missing_files = 0
    for idx, date_obj in enumerate(tqdm(dates, desc="读取SIF", leave=False)):
        file_path = get_SIF_file_path(date_obj)

        if file_path is None or not file_path.exists():
            missing_files += 1
            continue

        try:
            with rasterio.open(file_path) as src:
                data = src.read(1)
                # 将NODATA转为NaN
                data = np.where(data == NODATA_OUT, np.nan, data)
                sif_stack[idx] = data
        except Exception as e:
            tqdm.write(f"  ⚠ 读取失败: {file_path.name} - {str(e)}")
            continue

    print(f"  有效文件: {days_in_year - missing_files}/{days_in_year}")

    if missing_files > days_in_year * 0.3:
        print(f"  ⚠ 警告：缺失文件过多（{missing_files}/{days_in_year}）")

    # 初始化输出
    sos_map = np.full((height, width), NODATA_OUT, dtype=np.float32)
    pos_map = np.full((height, width), NODATA_OUT, dtype=np.float32)
    eos_map = np.full((height, width), NODATA_OUT, dtype=np.float32)

    # 逐像元提取物候
    print(f"  逐像元提取物候...")
    valid_pixels = 0

    for i in tqdm(range(height), desc="提取物候", leave=False):
        for j in range(width):
            if not mask[i, j]:
                continue

            # 提取时间序列
            sif_ts = sif_stack[:, i, j]

            # 提取物候
            sos, pos, eos = extract_phenology_threshold(sif_ts, smooth=True)

            if sos != NODATA_OUT:
                sos_map[i, j] = sos
                pos_map[i, j] = pos
                eos_map[i, j] = eos
                valid_pixels += 1

    print(f"  ✓ 有效像元数: {valid_pixels}")

    # 质量检查
    if valid_pixels > 0:
        sos_valid = sos_map[sos_map != NODATA_OUT]
        pos_valid = pos_map[pos_map != NODATA_OUT]
        eos_valid = eos_map[eos_map != NODATA_OUT]

        print(f"  SOS范围: {np.min(sos_valid):.0f} - {np.max(sos_valid):.0f} DOY")
        print(f"  POS范围: {np.min(pos_valid):.0f} - {np.max(pos_valid):.0f} DOY")
        print(f"  EOS范围: {np.min(eos_valid):.0f} - {np.max(eos_valid):.0f} DOY")
        print(f"  平均LSP: {np.mean(pos_valid - sos_valid):.1f} 天")

    return sos_map, pos_map, eos_map, profile

def write_geotiff(file_path, data, profile):
    """写入GeoTIFF"""
    profile.update(dtype=rasterio.float32, count=1, compress='lzw', nodata=NODATA_OUT)
    with rasterio.open(file_path, 'w', **profile) as dst:
        dst.write(data.astype(np.float32), 1)

# ==================== 主处理流程 ====================
def main():
    print("\n" + "="*70)
    print("Module 01: SIF物候提取")
    print("="*70)

    years = list(range(YEAR_START, YEAR_END + 1))

    # 读取掩膜
    print("\n[1/2] 读取掩膜...")
    mask_file = ROOT / "Wang2025_Analysis" / "masks" / "combined_mask.tif"

    if not mask_file.exists():
        print(f"  ✗ 错误：掩膜文件不存在: {mask_file}")
        print("  请先运行 Module 00: 00_data_preparation.py")
        return

    with rasterio.open(mask_file) as src:
        mask = src.read(1).astype(bool)

    print(f"  ✓ 有效像元数: {np.sum(mask)}")

    # 逐年提取物候
    print(f"\n[2/2] 逐年提取物候 ({YEAR_START}-{YEAR_END})...")

    for year in years:
        result = process_year(year, mask)

        if result[0] is None:
            print(f"  ✗ 跳过年份 {year}")
            continue

        sos_map, pos_map, eos_map, profile = result

        # 保存结果
        write_geotiff(OUTPUT_DIR / f"SOS_{year}.tif", sos_map, profile)
        write_geotiff(OUTPUT_DIR / f"POS_{year}.tif", pos_map, profile)
        write_geotiff(OUTPUT_DIR / f"EOS_{year}.tif", eos_map, profile)

    print("\n" + "="*70)
    print("✓ 物候提取完成！")
    print(f"输出目录: {OUTPUT_DIR}")
    print("="*70)

# ==================== 主程序 ====================
if __name__ == "__main__":
    main()
