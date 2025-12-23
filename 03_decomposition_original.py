#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module 03: 原版TRc分解方法
按照Wang (2025)原文方法：
  TRc_y = Σ[SOS to POS] TR(t)
  TRpheno_y = TRc_av × (LSP_y - LSP_av) / LSP_av
  TRproduct_y = TRc_y - TRc_av - TRpheno_y
"""

import numpy as np
import rasterio
from pathlib import Path
from rasterio.windows import Window
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timedelta

# ==================== 全局配置 ====================
ROOT = Path(r"I:\F\Data4")
PHENO_DIR = ROOT / "Phenology_Output_1" / "GPP_phenology"  # GPP物候（物候代码输出）
TRC_DIR = ROOT / "Wang2025_Analysis" / "TRc_annual"
OUTPUT_DIR = ROOT / "Wang2025_Analysis" / "Decomposition"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

YEAR_START = 1982
YEAR_END = 2018
BLOCK_SIZE = 128
MAX_WORKERS = 2
NODATA_OUT = -9999.0

# ==================== 辅助函数 ====================
def read_geotiff(file_path):
    """读取单波段GeoTIFF"""
    with rasterio.open(file_path) as src:
        data = src.read(1)
        profile = src.profile.copy()
    return data, profile

def write_geotiff(file_path, data, profile):
    """写入单波段GeoTIFF"""
    profile.update(dtype=rasterio.float32, count=1, compress='lzw', nodata=NODATA_OUT)
    with rasterio.open(file_path, 'w', **profile) as dst:
        dst.write(data.astype(np.float32), 1)

def calculate_LSP(sos, pos):
    """计算生长季长度 LSP = POS - SOS"""
    valid = (sos > 0) & (pos > 0) & (pos >= sos)
    lsp = np.where(valid, pos - sos, NODATA_OUT)
    return lsp

# ==================== 核心函数 ====================
def calculate_multiyear_mean(years, var_name, input_dir):
    """
    计算多年平均值（TRc_av, LSP_av）

    Parameters:
    -----------
    years : list
        年份列表
    var_name : str
        变量名 ('TRc' or 'LSP')
    input_dir : Path
        输入目录

    Returns:
    --------
    mean_array : ndarray
        多年平均值
    profile : dict
        栅格配置
    """
    print(f"\n=== 计算 {var_name}_av (多年平均) ===")

    stack = []
    for year in tqdm(years, desc=f"读取 {var_name}"):
        if var_name == 'LSP':
            # 物候代码输出格式：小写sos_gpp, pos_doy_gpp
            sos, _ = read_geotiff(PHENO_DIR / f"sos_gpp_{year}.tif")
            pos, profile = read_geotiff(PHENO_DIR / f"pos_doy_gpp_{year}.tif")
            data = calculate_LSP(sos, pos)
        else:  # TRc
            data, profile = read_geotiff(input_dir / f"{var_name}_{year}.tif")

        # 将NODATA替换为NaN用于计算
        data_masked = np.where(data == NODATA_OUT, np.nan, data)
        stack.append(data_masked)

    stack = np.stack(stack, axis=0)  # (n_years, H, W)

    # 沿时间轴计算平均值（忽略NaN）
    with np.errstate(invalid='ignore', divide='ignore'):
        mean_array = np.nanmean(stack, axis=0)
        mean_array = np.where(np.isnan(mean_array), NODATA_OUT, mean_array)

    print(f"✓ {var_name}_av 有效像元数: {np.sum(mean_array != NODATA_OUT)}")
    return mean_array, profile

def decompose_TRc_original(year, TRc_av, LSP_av, mask):
    """
    原版分解方法（Wang 2025）

    Parameters:
    -----------
    year : int
        年份
    TRc_av : ndarray
        多年平均TRc
    LSP_av : ndarray
        多年平均LSP
    mask : ndarray
        有效掩膜（≥30°N + 森林）

    Returns:
    --------
    TRpheno, TRproduct : ndarray
        分解后的两个分量
    """
    # 读取当年数据
    TRc_y, _ = read_geotiff(TRC_DIR / f"TRc_{year}.tif")
    # 物候代码输出格式：sos_gpp, pos_doy_gpp
    sos_y, _ = read_geotiff(PHENO_DIR / f"sos_gpp_{year}.tif")
    pos_y, _ = read_geotiff(PHENO_DIR / f"pos_doy_gpp_{year}.tif")
    LSP_y = calculate_LSP(sos_y, pos_y)

    # 构建有效掩膜
    valid = (
        (TRc_y != NODATA_OUT) &
        (TRc_av != NODATA_OUT) &
        (LSP_y != NODATA_OUT) &
        (LSP_av != NODATA_OUT) &
        (LSP_av > 0) &  # 避免除零
        mask
    )

    # 初始化输出
    TRpheno = np.full_like(TRc_y, NODATA_OUT, dtype=np.float32)
    TRproduct = np.full_like(TRc_y, NODATA_OUT, dtype=np.float32)

    # 原版公式（Wang 2025 Eq. 3-4）
    # TRpheno_y = TRc_av × (LSP_y - LSP_av) / LSP_av
    # TRproduct_y = TRc_y - TRc_av - TRpheno_y
    with np.errstate(divide='ignore', invalid='ignore'):
        TRpheno[valid] = TRc_av[valid] * (LSP_y[valid] - LSP_av[valid]) / LSP_av[valid]
        TRproduct[valid] = TRc_y[valid] - TRc_av[valid] - TRpheno[valid]

    return TRpheno, TRproduct

def process_all_years():
    """主处理流程"""
    print("\n" + "="*70)
    print("Wang (2025) 原版TRc分解")
    print("="*70)

    years = list(range(YEAR_START, YEAR_END + 1))

    # Step 1: 读取掩膜
    print("\n[1/4] 读取掩膜...")
    mask, profile = read_geotiff(OUTPUT_DIR.parent / "masks" / "combined_mask.tif")
    mask = mask.astype(bool)
    print(f"  有效像元数: {np.sum(mask)}")

    # Step 2: 计算多年平均值
    print("\n[2/4] 计算多年平均值...")
    TRc_av, _ = calculate_multiyear_mean(years, 'TRc', TRC_DIR)
    LSP_av, _ = calculate_multiyear_mean(years, 'LSP', PHENO_DIR)

    # 保存多年平均值
    write_geotiff(OUTPUT_DIR / "TRc_av.tif", TRc_av, profile)
    write_geotiff(OUTPUT_DIR / "LSP_av.tif", LSP_av, profile)
    print(f"  ✓ 已保存: TRc_av.tif, LSP_av.tif")

    # Step 3: 逐年分解
    print(f"\n[3/4] 逐年分解 ({YEAR_START}-{YEAR_END})...")
    for year in tqdm(years, desc="分解进度"):
        TRpheno, TRproduct = decompose_TRc_original(year, TRc_av, LSP_av, mask)

        # 保存结果
        write_geotiff(OUTPUT_DIR / f"TRpheno_{year}.tif", TRpheno, profile)
        write_geotiff(OUTPUT_DIR / f"TRproduct_{year}.tif", TRproduct, profile)

    # Step 4: 质量检查
    print("\n[4/4] 质量检查...")
    test_year = 2000
    TRc_test, _ = read_geotiff(TRC_DIR / f"TRc_{test_year}.tif")
    TRpheno_test, _ = read_geotiff(OUTPUT_DIR / f"TRpheno_{test_year}.tif")
    TRproduct_test, _ = read_geotiff(OUTPUT_DIR / f"TRproduct_{test_year}.tif")

    valid = mask & (TRc_test != NODATA_OUT)
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

    print("\n" + "="*70)
    print("✓ 原版分解完成！")
    print(f"输出目录: {OUTPUT_DIR}")
    print("="*70)

# ==================== 主程序 ====================
if __name__ == "__main__":
    process_all_years()
