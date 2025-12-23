#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module 05: 统计分析
1. Sen趋势 + Mann-Kendall显著性检验
2. 归因分析（标准化多元回归）
3. 15年滑动窗口敏感性分析
"""

import numpy as np
import rasterio
from pathlib import Path
from tqdm import tqdm
from scipy import stats
from scipy.stats import mannwhitall
import warnings
warnings.filterwarnings('ignore')

# ==================== 全局配置 ====================
ROOT = Path(r"I:\F\Data4")
DECOMP_DIR = ROOT / "Wang2025_Analysis" / "Decomposition"
PHENO_DIR = ROOT / "Phenology_Output_1" / "SIF_phenology"
SIF_DIR = ROOT / "SIF_Data" / "CSIF_annual"  # 年度SIF总量
SM_DIR = ROOT / "Meteorological Data" / "GLEAM" / "Annual" / "SMrz"
OUTPUT_DIR = ROOT / "Wang2025_Analysis" / "Statistical_Analysis"
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

def sen_slope(x, y):
    """
    Sen斜率估计
    x: 自变量（年份）
    y: 因变量（观测值）
    返回: (slope, intercept)
    """
    n = len(y)
    if n < 3:
        return np.nan, np.nan

    slopes = []
    for i in range(n):
        for j in range(i + 1, n):
            if x[j] != x[i]:
                slopes.append((y[j] - y[i]) / (x[j] - x[i]))

    if len(slopes) == 0:
        return np.nan, np.nan

    slope = np.median(slopes)
    intercept = np.median(y - slope * x)
    return slope, intercept

def mann_kendall_test(y):
    """
    Mann-Kendall趋势检验
    返回: (Z统计量, p值)
    """
    n = len(y)
    if n < 3:
        return np.nan, np.nan

    s = 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            s += np.sign(y[j] - y[i])

    # 计算方差
    unique, counts = np.unique(y, return_counts=True)
    tie_sum = np.sum(counts * (counts - 1) * (2 * counts + 5))
    var_s = (n * (n - 1) * (2 * n + 5) - tie_sum) / 18.0

    if var_s <= 0:
        return np.nan, np.nan

    # Z统计量
    if s > 0:
        z = (s - 1) / np.sqrt(var_s)
    elif s < 0:
        z = (s + 1) / np.sqrt(var_s)
    else:
        z = 0

    # 双侧检验p值
    p_value = 2 * (1 - stats.norm.cdf(abs(z)))

    return z, p_value

# ==================== 1. 趋势分析 ====================
def calculate_trends(var_name, input_dir, mask):
    """
    计算Sen趋势 + MK检验

    Parameters:
    -----------
    var_name : str
        变量名（如 'TRpheno', 'TRproduct', 'SOS'）
    input_dir : Path
        输入目录
    mask : ndarray
        有效掩膜

    Returns:
    --------
    slope, p_value : ndarray
        Sen斜率和MK检验p值
    """
    print(f"\n  计算 {var_name} 趋势...")

    years = list(range(YEAR_START, YEAR_END + 1))
    n_years = len(years)

    # 读取多年数据
    first_file = input_dir / f"{var_name}_{years[0]}.tif"
    with rasterio.open(first_file) as src:
        profile = src.profile.copy()
        height, width = src.height, src.width

    data_stack = np.full((n_years, height, width), np.nan, dtype=np.float32)

    for idx, year in enumerate(tqdm(years, desc=f"读取{var_name}", leave=False)):
        file_path = input_dir / f"{var_name}_{year}.tif"
        with rasterio.open(file_path) as src:
            data = src.read(1)
            data_stack[idx] = np.where(data == NODATA_OUT, np.nan, data)

    # 初始化输出
    slope_map = np.full((height, width), NODATA_OUT, dtype=np.float32)
    p_value_map = np.full((height, width), NODATA_OUT, dtype=np.float32)

    x = np.array(years, dtype=np.float32)

    # 逐像元计算
    for i in tqdm(range(height), desc=f"{var_name} 趋势分析", leave=False):
        for j in range(width):
            if not mask[i, j]:
                continue

            y = data_stack[:, i, j]
            valid_mask = ~np.isnan(y)
            valid_count = np.sum(valid_mask)

            if valid_count < n_years * 0.6:  # 至少60%有效数据
                continue

            y_valid = y[valid_mask]
            x_valid = x[valid_mask]

            # Sen斜率
            slope, _ = sen_slope(x_valid, y_valid)

            # MK检验
            z, p_val = mann_kendall_test(y_valid)

            slope_map[i, j] = slope if not np.isnan(slope) else NODATA_OUT
            p_value_map[i, j] = p_val if not np.isnan(p_val) else NODATA_OUT

    return slope_map, p_value_map, profile

# ==================== 2. 归因分析 ====================
def standardize(data):
    """标准化（z-score）"""
    mean = np.nanmean(data, axis=0)
    std = np.nanstd(data, axis=0)
    std[std == 0] = 1  # 避免除零
    return (data - mean) / std

def attribution_analysis(mask):
    """
    标准化多元回归归因分析（Wang 2025 Eq. 5-6）

    Y = TRc (or TRpheno, TRproduct)
    X = [SOS, LSP, SIF, SM]

    归因系数：α_i = β_i × std(X_i) / std(Y)
    """
    print("\n[2] 归因分析...")

    years = list(range(YEAR_START, YEAR_END + 1))
    n_years = len(years)

    # 读取模板
    with rasterio.open(DECOMP_DIR / f"TRc_{years[0]}.tif") as src:
        profile = src.profile.copy()
        height, width = src.height, src.width

    # 变量列表
    response_vars = ['TRc', 'TRpheno', 'TRproduct']
    predictor_vars = ['SOS', 'LSP', 'SIF', 'SM']  # 简化版，实际应包括更多气象变量

    results = {}

    for response_var in response_vars:
        print(f"\n  因变量: {response_var}")

        # 读取因变量数据
        Y_stack = []
        for year in years:
            if response_var in ['TRpheno', 'TRproduct']:
                data, _ = read_geotiff(DECOMP_DIR / f"{response_var}_{year}.tif")
            else:  # TRc
                data, _ = read_geotiff(DECOMP_DIR.parent / "TRc_annual" / f"TRc_{year}.tif")
            Y_stack.append(np.where(data == NODATA_OUT, np.nan, data))
        Y_stack = np.stack(Y_stack, axis=0)  # (n_years, H, W)

        # 读取自变量数据
        X_stack = {}
        for pred_var in predictor_vars:
            data_list = []
            for year in years:
                if pred_var == 'LSP':
                    sos, _ = read_geotiff(PHENO_DIR / f"SOS_{year}.tif")
                    pos, _ = read_geotiff(PHENO_DIR / f"POS_{year}.tif")
                    data = np.where((sos > 0) & (pos > sos), pos - sos, NODATA_OUT)
                elif pred_var in ['SOS', 'POS']:
                    data, _ = read_geotiff(PHENO_DIR / f"{pred_var}_{year}.tif")
                elif pred_var == 'SIF':
                    data, _ = read_geotiff(SIF_DIR / f"SIF_annual_{year}.tif")
                elif pred_var == 'SM':
                    data, _ = read_geotiff(SM_DIR / f"SMrz_{year}.tif")
                else:
                    data = np.full((height, width), NODATA_OUT)

                data_list.append(np.where(data == NODATA_OUT, np.nan, data))
            X_stack[pred_var] = np.stack(data_list, axis=0)  # (n_years, H, W)

        # 逐像元回归
        attribution_maps = {var: np.full((height, width), NODATA_OUT, dtype=np.float32)
                           for var in predictor_vars}

        for i in tqdm(range(height), desc=f"归因分析 {response_var}", leave=False):
            for j in range(width):
                if not mask[i, j]:
                    continue

                y = Y_stack[:, i, j]
                if np.sum(~np.isnan(y)) < n_years * 0.6:
                    continue

                # 构建设计矩阵
                X_matrix = np.column_stack([X_stack[var][:, i, j] for var in predictor_vars])

                # 移除含NaN的行
                valid_rows = ~np.isnan(X_matrix).any(axis=1) & ~np.isnan(y)
                if np.sum(valid_rows) < 10:
                    continue

                X_valid = X_matrix[valid_rows]
                y_valid = y[valid_rows]

                # 标准化
                X_std = standardize(X_valid)
                y_std = standardize(y_valid.reshape(-1, 1)).flatten()

                # 多元线性回归
                try:
                    beta, _, _, _ = np.linalg.lstsq(X_std, y_std, rcond=None)

                    # 计算归因系数（使用标准化系数）
                    for idx, var in enumerate(predictor_vars):
                        attribution_maps[var][i, j] = beta[idx]

                except np.linalg.LinAlgError:
                    continue

        # 保存结果
        attr_dir = OUTPUT_DIR / "Attribution" / response_var
        attr_dir.mkdir(parents=True, exist_ok=True)

        for var in predictor_vars:
            write_geotiff(attr_dir / f"attribution_{var}.tif", attribution_maps[var], profile)

        results[response_var] = attribution_maps

    print("\n  ✓ 归因分析完成")
    return results

# ==================== 3. 滑动窗口敏感性 ====================
def moving_window_sensitivity(window_size=15, mask=None):
    """
    15年滑动窗口敏感性分析
    计算每个窗口内 TRc ~ SOS 的偏相关系数
    """
    print(f"\n[3] {window_size}年滑动窗口敏感性分析...")

    years = list(range(YEAR_START, YEAR_END + 1))
    n_years = len(years)

    if n_years < window_size:
        print(f"  ⚠ 警告：总年数({n_years})小于窗口大小({window_size})，跳过")
        return

    # 读取模板
    with rasterio.open(DECOMP_DIR / f"TRc_{years[0]}.tif") as src:
        profile = src.profile.copy()
        height, width = src.height, src.width

    # 读取所有年份数据
    print("  读取数据...")
    TRc_stack = []
    SOS_stack = []

    for year in tqdm(years, desc="读取数据", leave=False):
        trc, _ = read_geotiff(DECOMP_DIR.parent / "TRc_annual" / f"TRc_{year}.tif")
        sos, _ = read_geotiff(PHENO_DIR / f"SOS_{year}.tif")

        TRc_stack.append(np.where(trc == NODATA_OUT, np.nan, trc))
        SOS_stack.append(np.where(sos == NODATA_OUT, np.nan, sos))

    TRc_stack = np.stack(TRc_stack, axis=0)  # (n_years, H, W)
    SOS_stack = np.stack(SOS_stack, axis=0)

    # 计算滑动窗口
    n_windows = n_years - window_size + 1
    sensitivity_maps = np.full((n_windows, height, width), NODATA_OUT, dtype=np.float32)

    for win_idx in tqdm(range(n_windows), desc="滑动窗口"):
        start_year = YEAR_START + win_idx
        end_year = start_year + window_size - 1

        for i in range(height):
            for j in range(width):
                if mask is not None and not mask[i, j]:
                    continue

                trc_win = TRc_stack[win_idx:win_idx + window_size, i, j]
                sos_win = SOS_stack[win_idx:win_idx + window_size, i, j]

                valid = ~np.isnan(trc_win) & ~np.isnan(sos_win)
                if np.sum(valid) < window_size * 0.6:
                    continue

                # Pearson相关系数
                r, p = stats.pearsonr(sos_win[valid], trc_win[valid])

                if p < 0.05:  # 仅保留显著相关
                    sensitivity_maps[win_idx, i, j] = r

    # 保存每个窗口的结果
    sens_dir = OUTPUT_DIR / "Moving_Window_Sensitivity"
    sens_dir.mkdir(parents=True, exist_ok=True)

    for win_idx in range(n_windows):
        start_year = YEAR_START + win_idx
        end_year = start_year + window_size - 1
        filename = f"sensitivity_TRc_SOS_{start_year}-{end_year}.tif"
        write_geotiff(sens_dir / filename, sensitivity_maps[win_idx], profile)

    print(f"  ✓ 滑动窗口分析完成，共{n_windows}个窗口")

# ==================== 主程序 ====================
def main():
    print("\n" + "="*70)
    print("统计分析模块")
    print("="*70)

    # 读取掩膜
    print("\n[0] 读取掩膜...")
    mask, profile = read_geotiff(DECOMP_DIR.parent / "masks" / "combined_mask.tif")
    mask = mask.astype(bool)

    # 1. 趋势分析
    print("\n[1] 趋势分析...")
    trend_dir = OUTPUT_DIR / "Trends"
    trend_dir.mkdir(parents=True, exist_ok=True)

    variables = {
        'SOS': PHENO_DIR,
        'LSP': PHENO_DIR,  # 需要特殊处理
        'TRpheno': DECOMP_DIR,
        'TRproduct': DECOMP_DIR
    }

    for var_name, input_dir in variables.items():
        if var_name == 'LSP':
            # 特殊处理LSP（需从SOS和POS计算）
            continue

        slope, p_value, _ = calculate_trends(var_name, input_dir, mask)

        write_geotiff(trend_dir / f"{var_name}_slope.tif", slope, profile)
        write_geotiff(trend_dir / f"{var_name}_pvalue.tif", p_value, profile)

    # 2. 归因分析
    attribution_results = attribution_analysis(mask)

    # 3. 滑动窗口敏感性
    moving_window_sensitivity(window_size=15, mask=mask)

    print("\n" + "="*70)
    print("✓ 统计分析完成！")
    print(f"输出目录: {OUTPUT_DIR}")
    print("="*70)

if __name__ == "__main__":
    main()
