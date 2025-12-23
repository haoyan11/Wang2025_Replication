#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module 04: 改进版TRc分解方法（补充材料）
方法1：端点法（Endpoint Decomposition）
方法2：速率法（Rate-based Decomposition）
方法3：SIF窗口替代法
"""

import numpy as np
import rasterio
from pathlib import Path
from tqdm import tqdm
from datetime import datetime, timedelta

# ==================== 全局配置 ====================
ROOT = Path(r"I:\F\Data4")
PHENO_DIR = ROOT / "Phenology_Output_1" / "SIF_phenology"
TRC_DIR = ROOT / "Wang2025_Analysis" / "TRc_annual"
TR_DAILY_DIR = ROOT / "Meteorological Data" / "GLEAM" / "Daily_0p5deg_TIF" / "Et"
SIF_DIR = ROOT / "SIF_Data" / "CSIF_daily"  # 假设路径
OUTPUT_DIR = ROOT / "Wang2025_Analysis" / "Decomposition_Improved"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

YEAR_START = 1982
YEAR_END = 2018
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
    """计算生长季长度"""
    valid = (sos > 0) & (pos > 0) & (pos >= sos)
    lsp = np.where(valid, pos - sos, NODATA_OUT)
    return lsp

def get_file_path_with_date(base_dir, prefix, date_obj):
    """
    获取日期文件路径（灵活命名）
    支持格式: YYYYMMDD.tif 或 YYYY_DOY.tif
    """
    # 优先尝试 YYYYMMDD 格式
    yyyymmdd = date_obj.strftime("%Y%m%d")
    file1 = base_dir / f"{prefix}_{yyyymmdd}.tif"
    if file1.exists():
        return file1

    # 尝试 YYYY_DOY 格式
    year = date_obj.year
    doy = date_obj.timetuple().tm_yday
    file2 = base_dir / f"{prefix}_{year}_{doy:03d}.tif"
    if file2.exists():
        return file2

    return None

# ==================== 方法1: 端点法 ====================
def decompose_endpoint(year, TRc_av, LSP_av, mask):
    """
    端点法（ChatGPT建议）：
    TRc_y^clim = Σ[SOS_y to POS_y] TR_clim(t)
      其中 TR_clim(t) 是气候平均态日均蒸腾
    TRpheno_endpoint = TRc_y^clim - TRc_av
    TRproduct_endpoint = TRc_y - TRc_y^clim
    """
    print(f"\n  处理年份: {year}")

    # 读取当年物候
    sos_y, _ = read_geotiff(PHENO_DIR / f"SOS_{year}.tif")
    pos_y, _ = read_geotiff(PHENO_DIR / f"POS_{year}.tif")
    TRc_y, profile = read_geotiff(TRC_DIR / f"TRc_{year}.tif")

    # 读取气候平均态日均TR（需提前计算）
    TR_clim_file = OUTPUT_DIR.parent / "climatology" / "TR_daily_climatology.tif"
    if not TR_clim_file.exists():
        print("  ⚠ 警告：缺少TR气候平均态文件，跳过端点法")
        return None, None

    with rasterio.open(TR_clim_file) as src:
        TR_clim = src.read()  # (366, H, W) 或 (365, H, W)

    # 计算TRc_y^clim
    TRc_y_clim = np.full_like(TRc_y, NODATA_OUT, dtype=np.float32)

    for i in range(TRc_y.shape[0]):
        for j in range(TRc_y.shape[1]):
            if not mask[i, j] or sos_y[i, j] == NODATA_OUT or pos_y[i, j] == NODATA_OUT:
                continue

            sos_doy = int(sos_y[i, j])
            pos_doy = int(pos_y[i, j])

            if sos_doy <= 0 or pos_doy < sos_doy:
                continue

            # 累加气候平均态
            doy_range = range(sos_doy - 1, pos_doy)  # DOY从1开始，数组从0开始
            TRc_y_clim[i, j] = np.sum(TR_clim[doy_range, i, j])

    # 分解
    valid = (TRc_y != NODATA_OUT) & (TRc_y_clim != NODATA_OUT) & (TRc_av != NODATA_OUT) & mask
    TRpheno_endpoint = np.full_like(TRc_y, NODATA_OUT)
    TRproduct_endpoint = np.full_like(TRc_y, NODATA_OUT)

    TRpheno_endpoint[valid] = TRc_y_clim[valid] - TRc_av[valid]
    TRproduct_endpoint[valid] = TRc_y[valid] - TRc_y_clim[valid]

    return TRpheno_endpoint, TRproduct_endpoint

# ==================== 方法2: 速率法 ====================
def decompose_rate_based(year, mask):
    """
    速率法（基于日均TR）：
    TRmean_y = TRc_y / LSP_y
    ΔTRmean = TRmean_y - TRmean_av
    ΔLSP = LSP_y - LSP_av
    TRpheno_rate = TRmean_av × ΔLSP
    TRproduct_rate = ΔTRmean × LSP_y
    """
    # 读取数据
    TRc_y, _ = read_geotiff(TRC_DIR / f"TRc_{year}.tif")
    sos_y, _ = read_geotiff(PHENO_DIR / f"SOS_{year}.tif")
    pos_y, profile = read_geotiff(PHENO_DIR / f"POS_{year}.tif")
    LSP_y = calculate_LSP(sos_y, pos_y)

    TRc_av, _ = read_geotiff(OUTPUT_DIR.parent / "Decomposition" / "TRc_av.tif")
    LSP_av, _ = read_geotiff(OUTPUT_DIR.parent / "Decomposition" / "LSP_av.tif")

    # 计算TRmean
    valid = (TRc_y != NODATA_OUT) & (LSP_y != NODATA_OUT) & (LSP_y > 0) & mask
    TRmean_y = np.full_like(TRc_y, NODATA_OUT)
    TRmean_y[valid] = TRc_y[valid] / LSP_y[valid]

    TRmean_av = np.full_like(TRc_av, NODATA_OUT)
    valid_av = (TRc_av != NODATA_OUT) & (LSP_av != NODATA_OUT) & (LSP_av > 0) & mask
    TRmean_av[valid_av] = TRc_av[valid_av] / LSP_av[valid_av]

    # 分解
    valid_final = valid & valid_av
    TRpheno_rate = np.full_like(TRc_y, NODATA_OUT)
    TRproduct_rate = np.full_like(TRc_y, NODATA_OUT)

    with np.errstate(divide='ignore', invalid='ignore'):
        delta_LSP = LSP_y[valid_final] - LSP_av[valid_final]
        delta_TRmean = TRmean_y[valid_final] - TRmean_av[valid_final]

        TRpheno_rate[valid_final] = TRmean_av[valid_final] * delta_LSP
        TRproduct_rate[valid_final] = delta_TRmean * LSP_y[valid_final]

    return TRpheno_rate, TRproduct_rate

# ==================== 方法3: SIF窗口替代 ====================
def calculate_SIF_alternative_window(year, window_type, mask):
    """
    SIF窗口替代法：
    window_type = 'early' (SOS+0~30天) 或 'mid' (SOS+30~90天)
    计算该窗口的SIF累积，作为生产力代理变量
    """
    sos_y, _ = read_geotiff(PHENO_DIR / f"SOS_{year}.tif")
    profile = rasterio.open(PHENO_DIR / f"SOS_{year}.tif").profile

    SIF_window = np.zeros_like(sos_y, dtype=np.float32)

    for doy in range(1, 367):
        date_obj = datetime(year, 1, 1) + timedelta(days=doy - 1)
        sif_file = get_file_path_with_date(SIF_DIR, "SIF", date_obj)

        if sif_file is None or not sif_file.exists():
            continue

        with rasterio.open(sif_file) as src:
            sif_daily = src.read(1)

        for i in range(sos_y.shape[0]):
            for j in range(sos_y.shape[1]):
                if not mask[i, j] or sos_y[i, j] == NODATA_OUT:
                    continue

                sos_doy = int(sos_y[i, j])

                if window_type == 'early':
                    if sos_doy <= doy <= sos_doy + 30:
                        SIF_window[i, j] += sif_daily[i, j]
                elif window_type == 'mid':
                    if sos_doy + 30 <= doy <= sos_doy + 90:
                        SIF_window[i, j] += sif_daily[i, j]

    SIF_window[~mask] = NODATA_OUT
    return SIF_window, profile

# ==================== 主处理流程 ====================
def process_improved_methods():
    """执行所有改进方法"""
    print("\n" + "="*70)
    print("改进版TRc分解方法（补充材料）")
    print("="*70)

    years = list(range(YEAR_START, YEAR_END + 1))

    # 读取掩膜
    print("\n读取掩膜...")
    mask, profile = read_geotiff(OUTPUT_DIR.parent / "masks" / "combined_mask.tif")
    mask = mask.astype(bool)

    TRc_av, _ = read_geotiff(OUTPUT_DIR.parent / "Decomposition" / "TRc_av.tif")
    LSP_av, _ = read_geotiff(OUTPUT_DIR.parent / "Decomposition" / "LSP_av.tif")

    # 方法1: 端点法
    print("\n[方法1] 端点法...")
    endpoint_dir = OUTPUT_DIR / "Method1_Endpoint"
    endpoint_dir.mkdir(exist_ok=True)

    for year in tqdm(years, desc="端点法"):
        TRpheno_ep, TRproduct_ep = decompose_endpoint(year, TRc_av, LSP_av, mask)
        if TRpheno_ep is not None:
            write_geotiff(endpoint_dir / f"TRpheno_endpoint_{year}.tif", TRpheno_ep, profile)
            write_geotiff(endpoint_dir / f"TRproduct_endpoint_{year}.tif", TRproduct_ep, profile)

    # 方法2: 速率法
    print("\n[方法2] 速率法...")
    rate_dir = OUTPUT_DIR / "Method2_Rate"
    rate_dir.mkdir(exist_ok=True)

    for year in tqdm(years, desc="速率法"):
        TRpheno_rate, TRproduct_rate = decompose_rate_based(year, mask)
        write_geotiff(rate_dir / f"TRpheno_rate_{year}.tif", TRpheno_rate, profile)
        write_geotiff(rate_dir / f"TRproduct_rate_{year}.tif", TRproduct_rate, profile)

    # 方法3: SIF窗口替代
    print("\n[方法3] SIF窗口替代法...")
    sif_window_dir = OUTPUT_DIR / "Method3_SIF_Windows"
    sif_window_dir.mkdir(exist_ok=True)

    for year in tqdm(years, desc="SIF窗口"):
        for window_type in ['early', 'mid']:
            SIF_win, _ = calculate_SIF_alternative_window(year, window_type, mask)
            write_geotiff(sif_window_dir / f"SIF_{window_type}_{year}.tif", SIF_win, profile)

    print("\n" + "="*70)
    print("✓ 改进方法完成！")
    print(f"输出目录: {OUTPUT_DIR}")
    print("="*70)

if __name__ == "__main__":
    process_improved_methods()
