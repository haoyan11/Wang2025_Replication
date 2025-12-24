#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module 05: 统计分析 - Wang (2025) Sections 3.2 & 3.3

Section 3.2: Increasing TRpheno but decreasing TRproduct with earlier spring phenology
    - 像元级线性回归: TRc/TRpheno/TRproduct vs dSOS (advance)

Section 3.3: Drivers of TRproduct decrease with spring phenology change
    - 偏相关归因分析（控制其他变量）
    - 15年滑动窗口偏相关演变
    - Theil-Sen趋势 + Mann-Kendall检验

Version: 2.0.0
Author: Wang2025 Replication Project

核心方法：
- ΔSOS = SOS_year - SOSav (标准异常定义，advance < 0)
- 偏相关（控制其他变量，Z-score，Wang 2025 Eq. 3）
- VIF > 10 的变量剔除
- 主驱动因子判定: |R| > 0.1
"""

import os
import numpy as np
import rasterio
from pathlib import Path
from tqdm import tqdm
from scipy import stats
from datetime import datetime, timedelta
from functools import lru_cache
from concurrent.futures import ProcessPoolExecutor, as_completed
import warnings
warnings.filterwarnings('ignore')

# 性能优化：Numba加速
try:
    from numba import jit, prange
    NUMBA_AVAILABLE = True
except ImportError:
    print("⚠️ 警告：未安装numba，将使用numpy版本（较慢）")
    print("  建议安装：pip install numba")
    NUMBA_AVAILABLE = False
    # 定义dummy装饰器
    def jit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator
    prange = range

# ==================== 全局配置 ====================
ROOT = Path(r"I:\F\Data4")
ANALYSIS_DIR = ROOT / "Wang2025_Analysis"
DECOMP_DIR = ANALYSIS_DIR / "Decomposition_T"
TRC_DIR = ANALYSIS_DIR / "TRc_annual_T"
PHENO_DIR = ROOT / "Phenology_Output_1" / "T_phenology"
GPP_DAILY_DIR = ROOT / "GLASS_GPP" / "GLASS_GPP_daily_interpolated"  # 保持不变，用于SIF代理
GPP_DAILY_FORMAT = "GPP_{date}.tif"  # {date} = YYYYMMDD
METEO_DIR = ROOT / "Meteorological Data"
OUTPUT_DIR = ANALYSIS_DIR / "Statistical_Analysis_T"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

YEAR_START = 1982
YEAR_END = 2018
NODATA_OUT = -9999.0
NODATA_ABS_MAX = 1e20
BLOCK_SIZE_3_3 = 64
USE_BLOCK_PARALLEL = True
MAX_WORKERS_3_3 = min(4, os.cpu_count() or 1)

# 季节定义
SPRING_MONTHS = [3, 4, 5]  # 3-5月
SUMMER_MONTHS = [6, 7, 8]  # 6-8月

DAILY_VAR_SPECS = {
    'Ta': {
        'dir': METEO_DIR / "ERA5_Land" / "Tem" / "Tem_Daily" / "Tem_Daily_2",
        'pattern': "ERA5L_T2mDaily_C_{date}.tif"
    },
    'Rs': {
        'dir': METEO_DIR / "ERA5_Land" / "DSW" / "DSW_Daily" / "DSW_Daily_2",
        'pattern': "ERA5L_SWDaily_MJ_{date}.tif"
    },
    'P': {
        'dir': METEO_DIR / "ERA5_Land" / "Pre" / "Pre_Daily" / "Pre_Daily_2",
        'pattern': "ERA5L_PrecipDaily_mm_{date}.tif"
    },
    'SMsurf': {
        'dir': METEO_DIR / "GLEAM" / "SMs" / "SMs_Daily_1",
        'pattern': "SMs_{date}.tif"
    },
    'SMroot': {
        'dir': METEO_DIR / "GLEAM" / "SMrz" / "SMrz_Daily_1",
        'pattern': "SMrz_{date}.tif"
    }
}

# ==================== 辅助函数 ====================
def _is_valid_value(value, nodata):
    """检查值是否有效（非NODATA）"""
    if nodata is None:
        return np.isfinite(value) & (np.abs(value) < NODATA_ABS_MAX)
    if np.isnan(nodata):
        return np.isfinite(value) & (np.abs(value) < NODATA_ABS_MAX)
    return (value != nodata) & np.isfinite(value) & (np.abs(value) < NODATA_ABS_MAX)

@lru_cache(maxsize=None)
def _has_daily_files(var_name, year):
    """快速检测某年日尺度文件是否存在（抽样几个DOY）"""
    spec = DAILY_VAR_SPECS.get(var_name)
    if spec is None:
        return False
    var_dir = spec['dir']
    pattern = spec['pattern']
    for doy in (1, 120, 240):
        date_str = (datetime(year, 1, 1) + timedelta(days=doy - 1)).strftime("%Y%m%d")
        if (var_dir / pattern.format(date=date_str)).exists():
            return True
    return False

@lru_cache(maxsize=None)
def _has_gpp_files(year):
    """快速检测某年日尺度GPP文件是否存在"""
    for doy in (90, 120, 180, 220):
        date_str = (datetime(year, 1, 1) + timedelta(days=doy - 1)).strftime("%Y%m%d")
        if (GPP_DAILY_DIR / GPP_DAILY_FORMAT.format(date=date_str)).exists():
            return True
    return False

def iter_valid_blocks(mask, block_size):
    """生成仅包含有效像元的块索引"""
    height, width = mask.shape
    for r0 in range(0, height, block_size):
        r1 = min(height, r0 + block_size)
        for c0 in range(0, width, block_size):
            c1 = min(width, c0 + block_size)
            block_mask = mask[r0:r1, c0:c1]
            if np.any(block_mask):
                yield r0, r1, c0, c1, block_mask

def read_geotiff(file_path):
    """读取单波段GeoTIFF，返回数据和profile（保持原始NODATA）"""
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

def get_doy_from_date(year, month, day):
    """计算儒略日（DOY）"""
    return datetime(year, month, day).timetuple().tm_yday

def sen_slope(x, y):
    """
    Theil-Sen斜率估计（非参数趋势分析）

    Parameters:
    -----------
    x : array_like
        自变量（通常为时间）
    y : array_like
        因变量

    Returns:
    --------
    slope, intercept : float
        斜率和截距
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
    Mann-Kendall趋势显著性检验

    Parameters:
    -----------
    y : array_like
        时间序列数据

    Returns:
    --------
    z, p_value : float
        Z统计量和双侧p值
    """
    n = len(y)
    if n < 3:
        return np.nan, np.nan

    s = 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            s += np.sign(y[j] - y[i])

    # 计算方差（考虑并列值）
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

def _sum_cum_range(cum_arr, start_doy, end_doy):
    """
    使用累积和快速计算区间和（1-based DOY，区间为 [start, end]）
    start_doy/end_doy 为二维数组
    """
    flat = cum_arr.reshape(cum_arr.shape[0], -1)
    start = start_doy.ravel()
    end = end_doy.ravel()
    idx = np.arange(start.size)

    end_vals = flat[end - 1, idx]
    start_vals = np.where(start > 1, flat[start - 2, idx], 0)
    return (end_vals - start_vals).reshape(start_doy.shape)

def linear_regression_maps(x_stack, y_stack, min_frac=0.6):
    """
    向量化线性回归：y ~ x

    Parameters:
    -----------
    x_stack, y_stack : ndarray
        形状 (n_years, H, W)
    min_frac : float
        最低有效样本比例

    Returns:
    --------
    slope_map, pvalue_map, r2_map : ndarray
    """
    n_years = x_stack.shape[0]
    valid = np.isfinite(x_stack) & np.isfinite(y_stack)
    n_valid = np.sum(valid, axis=0).astype(np.float32)

    x = np.where(valid, x_stack, np.nan)
    y = np.where(valid, y_stack, np.nan)

    mean_x = np.nanmean(x, axis=0)
    mean_y = np.nanmean(y, axis=0)

    dx = x - mean_x
    dy = y - mean_y

    cov_xy = np.nanmean(dx * dy, axis=0)
    var_x = np.nanmean(dx * dx, axis=0)
    var_y = np.nanmean(dy * dy, axis=0)

    with np.errstate(divide='ignore', invalid='ignore'):
        slope = cov_xy / var_x
        r = cov_xy / np.sqrt(var_x * var_y)
        r2 = r ** 2

        df = n_valid - 2
        t_stat = r * np.sqrt(df / (1 - r2))
        p_val = 2 * (1 - stats.t.cdf(np.abs(t_stat), df=df))

    invalid = (
        (n_valid < n_years * min_frac) |
        (n_valid < 3) |
        ~np.isfinite(slope) |
        ~np.isfinite(r2) |
        ~np.isfinite(p_val)
    )

    slope_map = np.full(slope.shape, NODATA_OUT, dtype=np.float32)
    pvalue_map = np.full(p_val.shape, NODATA_OUT, dtype=np.float32)
    r2_map = np.full(r2.shape, NODATA_OUT, dtype=np.float32)

    valid_mask = ~invalid
    slope_map[valid_mask] = slope[valid_mask]
    pvalue_map[valid_mask] = p_val[valid_mask]
    r2_map[valid_mask] = r2[valid_mask]

    return slope_map, pvalue_map, r2_map

def partial_corr_from_std(y_std, X_std):
    """
    基于相关矩阵的偏相关系数计算（输入已标准化）

    Returns:
    --------
    r : ndarray (p,)
    p_vals : ndarray (p,)
    """
    if X_std.ndim == 1:
        X_std = X_std.reshape(-1, 1)
    n, p = X_std.shape
    if p == 0:
        return None, None

    data = np.column_stack([y_std, X_std])
    if not np.isfinite(data).all():
        return None, None

    try:
        corr = np.corrcoef(data, rowvar=False)
        prec = np.linalg.inv(corr)
    except np.linalg.LinAlgError:
        return None, None

    denom = np.sqrt(prec[0, 0] * np.diag(prec)[1:])
    with np.errstate(divide='ignore', invalid='ignore'):
        r = -prec[0, 1:] / denom

    r = np.clip(r, -0.999999, 0.999999)
    df = n - p - 1
    if df <= 0:
        p_vals = np.full(p, np.nan, dtype=np.float32)
    else:
        t_stat = r * np.sqrt(df / (1 - r ** 2))
        p_vals = 2 * (1 - stats.t.cdf(np.abs(t_stat), df=df))

    return r.astype(np.float32), p_vals.astype(np.float32)

def _partial_corr_block_worker(args):
    """
    处理一个块的偏相关与VIF过滤（用于并行/分块）
    """
    (r0, r1, c0, c1, block_mask, predictor_vars, X_block, Y_block, min_rows) = args
    block_h = r1 - r0
    block_w = c1 - c0

    partial_r_block = {var: np.full((block_h, block_w), NODATA_OUT, dtype=np.float32)
                       for var in predictor_vars}
    partial_p_block = {var: np.full((block_h, block_w), NODATA_OUT, dtype=np.float32)
                       for var in predictor_vars}
    vif_block = {var: np.zeros((block_h, block_w), dtype=bool)
                 for var in predictor_vars}
    r2_block = np.full((block_h, block_w), NODATA_OUT, dtype=np.float32)

    for i_rel, j_rel in np.argwhere(block_mask):
        y = Y_block[:, i_rel, j_rel]
        if np.sum(~np.isnan(y)) < min_rows:
            continue

        X_matrix = np.column_stack([X_block[var][:, i_rel, j_rel] for var in predictor_vars])
        valid_rows = ~np.isnan(X_matrix).any(axis=1) & ~np.isnan(y)
        if np.sum(valid_rows) < min_rows:
            continue

        X_valid = X_matrix[valid_rows]
        y_valid = y[valid_rows]

        X_std = standardize(X_valid)
        y_std = standardize(y_valid.reshape(-1, 1)).flatten()

        current_vars = list(range(len(predictor_vars)))
        X_filtered = X_std.copy()

        while len(current_vars) > 1:
            vif_values = calculate_vif(X_filtered)
            max_vif_idx = np.argmax(vif_values)
            if vif_values[max_vif_idx] > 10:
                X_filtered = np.delete(X_filtered, max_vif_idx, axis=1)
                current_vars.pop(max_vif_idx)
            else:
                break

        try:
            r_vals, p_vals = partial_corr_from_std(y_std, X_filtered)
            if r_vals is None:
                continue

            for idx, var_idx in enumerate(current_vars):
                var_name = predictor_vars[var_idx]
                partial_r_block[var_name][i_rel, j_rel] = r_vals[idx]
                partial_p_block[var_name][i_rel, j_rel] = p_vals[idx]
                vif_block[var_name][i_rel, j_rel] = True

            beta = np.linalg.lstsq(X_filtered, y_std, rcond=None)[0]
            ss_res = np.sum((y_std - X_filtered @ beta) ** 2)
            ss_tot = np.sum(y_std ** 2)
            r2_block[i_rel, j_rel] = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

        except np.linalg.LinAlgError:
            continue

    return r0, r1, c0, c1, partial_r_block, partial_p_block, vif_block, r2_block

def _partial_corr_window_block_worker(args):
    """
    处理一个块的滑动窗口偏相关（不做VIF）
    """
    (r0, r1, c0, c1, block_mask, predictor_vars, X_block, Y_block, min_rows) = args
    block_h = r1 - r0
    block_w = c1 - c0

    partial_r_block = {var: np.full((block_h, block_w), NODATA_OUT, dtype=np.float32)
                       for var in predictor_vars}

    for i_rel, j_rel in np.argwhere(block_mask):
        y = Y_block[:, i_rel, j_rel]
        if np.sum(~np.isnan(y)) < min_rows:
            continue

        X_matrix = np.column_stack([X_block[var][:, i_rel, j_rel] for var in predictor_vars])
        valid_rows = ~np.isnan(X_matrix).any(axis=1) & ~np.isnan(y)
        if np.sum(valid_rows) < min_rows:
            continue

        X_valid = X_matrix[valid_rows]
        y_valid = y[valid_rows]

        X_std = standardize(X_valid)
        y_std = standardize(y_valid.reshape(-1, 1)).flatten()

        try:
            r_vals, _ = partial_corr_from_std(y_std, X_std)
            if r_vals is None:
                continue

            for idx, var in enumerate(predictor_vars):
                partial_r_block[var][i_rel, j_rel] = r_vals[idx]

        except np.linalg.LinAlgError:
            continue

    return r0, r1, c0, c1, partial_r_block

def _trend_block_worker(args):
    """
    处理一个块的趋势分析（Theil-Sen + MK）
    """
    (r0, r1, c0, c1, block_mask, sens_block, window_years, min_frac) = args
    block_h = r1 - r0
    block_w = c1 - c0
    slope_block = np.full((block_h, block_w), NODATA_OUT, dtype=np.float32)
    p_block = np.full((block_h, block_w), NODATA_OUT, dtype=np.float32)

    n_windows = sens_block.shape[0]

    for i_rel, j_rel in np.argwhere(block_mask):
        sens_series = sens_block[:, i_rel, j_rel]
        valid = (sens_series != NODATA_OUT) & np.isfinite(sens_series)
        if np.sum(valid) < n_windows * min_frac:
            continue

        sens_valid = sens_series[valid]
        window_years_valid = window_years[valid]

        slope, _ = sen_slope(window_years_valid, sens_valid)
        _, p_val = mann_kendall_test(sens_valid)

        if np.isfinite(slope):
            slope_block[i_rel, j_rel] = slope
        if np.isfinite(p_val):
            p_block[i_rel, j_rel] = p_val

    return r0, r1, c0, c1, slope_block, p_block

def standardize(data):
    """
    Z-score标准化

    Parameters:
    -----------
    data : ndarray
        输入数据 (n_samples, n_features) 或 (n_samples,)

    Returns:
    --------
    data_std : ndarray
        标准化后的数据
    """
    mean = np.nanmean(data, axis=0)
    std = np.nanstd(data, axis=0)
    std[std == 0] = 1  # 避免除零
    return (data - mean) / std

def calculate_vif(X):
    """
    计算方差膨胀因子（Variance Inflation Factor, VIF）

    VIF_i = 1 / (1 - R_i²)
    其中 R_i² 是第i个变量对其他变量回归的决定系数

    注：输入变量已做Z-score标准化，因此不显式加入截距项。

    Parameters:
    -----------
    X : ndarray
        设计矩阵 (n_samples, n_features)

    Returns:
    --------
    vif : ndarray
        每个特征的VIF值 (n_features,)
    """
    n_features = X.shape[1]
    vif = np.zeros(n_features)

    for i in range(n_features):
        # 将第i个变量作为因变量
        y = X[:, i]
        # 其他变量作为自变量
        X_others = np.delete(X, i, axis=1)

        # 线性回归计算R²
        try:
            beta = np.linalg.lstsq(X_others, y, rcond=None)[0]
            y_pred = X_others @ beta
            ss_res = np.sum((y - y_pred) ** 2)
            ss_tot = np.sum((y - np.mean(y)) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

            # VIF计算
            vif[i] = 1 / (1 - r_squared) if r_squared < 0.9999 else np.inf
        except:
            vif[i] = np.inf

    return vif

def calculate_seasonal_gpp(year, season='spring'):
    """
    计算季节平均GPP（从日GPP数据）

    Parameters:
    -----------
    year : int
        年份
    season : str
        'spring' (3-5月) 或 'summer' (6-8月)

    Returns:
    --------
    sif_seasonal : ndarray
        季节平均GPP (H, W)
    """
    if not _has_gpp_files(year):
        return None

    months = SPRING_MONTHS if season == 'spring' else SUMMER_MONTHS

    gpp_list = []
    for month in months:
        for day in range(1, 32):
            try:
                date_str = datetime(year, month, day).strftime("%Y%m%d")
                gpp_file = GPP_DAILY_DIR / GPP_DAILY_FORMAT.format(date=date_str)

                if gpp_file.exists():
                    data, _, nodata = read_geotiff(gpp_file)
                    valid_data = np.where(_is_valid_value(data, nodata), data, np.nan)
                    gpp_list.append(valid_data)
            except ValueError:
                continue

    if len(gpp_list) == 0:
        return None

    gpp_stack = np.stack(gpp_list, axis=0)
    gpp_seasonal = np.nanmean(gpp_stack, axis=0)

    return gpp_seasonal

def calculate_lsp_period_average(var_name, year, sos_map, pos_map):
    """
    计算LSP期间的变量平均值（Wang 2025 Section 2.2.2）

    Parameters:
    -----------
    var_name : str
        变量名 ('Ta', 'Rs', 'P', 'SMsurf', 'SMroot')
    year : int
        年份
    sos_map : ndarray
        SOS地图 (H, W)，单位：DOY
    pos_map : ndarray
        POS地图 (H, W)，单位：DOY

    Returns:
    --------
    lsp_avg : ndarray
        LSP期间平均值 (H, W)
    """
    spec = DAILY_VAR_SPECS.get(var_name)
    if spec is None:
        raise ValueError(f"Unknown variable: {var_name}")

    var_dir = spec['dir']
    pattern = spec['pattern']
    height, width = sos_map.shape

    if not _has_daily_files(var_name, year):
        return np.full((height, width), np.nan, dtype=np.float32)

    sos_int = np.rint(sos_map).astype(np.int32)
    pos_int = np.rint(pos_map).astype(np.int32)

    valid = (
        np.isfinite(sos_map) & np.isfinite(pos_map) &
        (sos_int > 0) & (pos_int > 0) &
        (pos_int >= sos_int) &
        (sos_int <= 365) & (pos_int <= 365)
    )

    if not np.any(valid):
        return np.full((height, width), np.nan, dtype=np.float32)

    doy_start = int(np.nanmin(sos_int[valid]))
    doy_end = int(np.nanmax(pos_int[valid]))
    doy_start = max(1, min(365, doy_start))
    doy_end = max(1, min(365, doy_end))

    total_sum = np.zeros((height, width), dtype=np.float32)
    total_cnt = np.zeros((height, width), dtype=np.int16)

    for doy in range(doy_start, doy_end + 1):
        in_window = valid & (sos_int <= doy) & (pos_int >= doy)
        if not np.any(in_window):
            continue

        date_str = (datetime(year, 1, 1) + timedelta(days=doy - 1)).strftime("%Y%m%d")
        var_file = var_dir / pattern.format(date=date_str)
        if not var_file.exists():
            continue

        data, _, nodata = read_geotiff(var_file)
        valid_data = _is_valid_value(data, nodata)
        use_mask = in_window & valid_data

        if not np.any(use_mask):
            continue

        total_sum[use_mask] += data[use_mask].astype(np.float32)
        total_cnt[use_mask] += 1

    window_len = pos_int - sos_int + 1
    lsp_avg = np.full((height, width), np.nan, dtype=np.float32)
    good = valid & (total_cnt >= 0.6 * window_len)
    lsp_avg[good] = total_sum[good] / total_cnt[good]

    return lsp_avg

# ==================== Section 3.2: ΔSOS回归分析 ====================
def section_3_2_phenology_impact(mask):
    """
    Section 3.2: Increasing TRpheno but decreasing TRproduct with earlier spring phenology

    方法：
    1. 计算 ΔSOS = SOS_year - SOSav (标准异常定义，advance < 0)
    2. 像元级线性回归：
       - TRc ~ ΔSOS
       - TRpheno ~ ΔSOS
       - TRproduct ~ ΔSOS
    3. 输出回归斜率和显著性

    预期结果（论文 Figure 4，advance < 0）：
    - TRc: -0.76 ± 0.97 mm/day per day of ΔSOS (负值表示SOS提前时TRc增加)
    - TRpheno: -1.80 ± 0.44 mm/day per day of ΔSOS (负值表示SOS提前时TRpheno增加)
    - TRproduct: +1.04 ± 0.80 mm/day per day of ΔSOS (正值表示SOS提前时TRproduct减少)

    Parameters:
    -----------
    mask : ndarray
        有效掩膜 (H, W)

    Returns:
    --------
    None (结果保存到文件)
    """
    print("\n" + "="*70)
    print("Section 3.2: TRc Components vs ΔSOS 回归分析")
    print("="*70)

    years = list(range(YEAR_START, YEAR_END + 1))
    n_years = len(years)

    # 读取模板
    data_first, profile, _ = read_geotiff(TRC_DIR / f"TRc_{years[0]}.tif")
    height, width = data_first.shape

    # Step 1: 计算多年平均SOS（SOSav），同步读取POS用于Trate
    print("\n[Step 1] 计算多年平均 SOS (SOSav)...")
    sos_stack = []
    pos_stack = []

    for year in tqdm(years, desc="读取SOS/POS"):
        sos_data, _, sos_nodata = read_geotiff(PHENO_DIR / f"sos_t_{year}.tif")
        pos_data, _, pos_nodata = read_geotiff(PHENO_DIR / f"pos_doy_t_{year}.tif")
        sos_valid = np.where(_is_valid_value(sos_data, sos_nodata), sos_data, np.nan)
        pos_valid = np.where(_is_valid_value(pos_data, pos_nodata), pos_data, np.nan)
        sos_stack.append(sos_valid)
        pos_stack.append(pos_valid)

    sos_stack = np.stack(sos_stack, axis=0)  # (n_years, H, W)
    sos_av = np.nanmean(sos_stack, axis=0)  # (H, W)

    # Step 2: 计算ΔSOS（标准异常定义，advance < 0）
    print("\n[Step 2] 计算 ΔSOS = SOS_year - SOSav (标准异常定义)...")
    delta_sos_stack = sos_stack - sos_av  # (n_years, H, W)

    # Step 2b: 计算GLS（用于Trate）
    print("\n[Step 2b] 计算 GLS = POS - SOS + 1（用于Trate）...")
    gls_stack = np.full_like(sos_stack, np.nan, dtype=np.float32)
    for idx in range(n_years):
        sos_y = sos_stack[idx]
        pos_y = pos_stack[idx]
        sos_int = np.rint(sos_y).astype(np.int32)
        pos_int = np.rint(pos_y).astype(np.int32)
        valid = (
            np.isfinite(sos_y) & np.isfinite(pos_y) &
            (sos_int > 0) & (pos_int > 0) &
            (pos_int >= sos_int)
        )
        gls = np.full((height, width), np.nan, dtype=np.float32)
        gls[valid] = (pos_int[valid] - sos_int[valid] + 1).astype(np.float32)
        gls_stack[idx] = gls

    # Step 3: 读取响应变量（TRc, TRpheno, TRproduct, Trate）
    print("\n[Step 3] 读取响应变量...")
    response_vars = ['TRc', 'TRpheno', 'TRproduct', 'Trate']
    response_data = {}

    for resp_var in ['TRc', 'TRpheno', 'TRproduct']:
        data_stack = []
        for year in tqdm(years, desc=f"读取{resp_var}", leave=False):
            if resp_var == 'TRc':
                data, _, nodata = read_geotiff(TRC_DIR / f"TRc_{year}.tif")
            else:
                data, _, nodata = read_geotiff(DECOMP_DIR / f"{resp_var}_{year}.tif")

            data_stack.append(np.where(_is_valid_value(data, nodata), data, np.nan))

        response_data[resp_var] = np.stack(data_stack, axis=0)  # (n_years, H, W)

    # 计算 Trate = TRc / GLS（mm/day）
    print("\n  计算 Trate = TRc / GLS...")
    trc_stack = response_data['TRc']
    trate_stack = np.full_like(trc_stack, np.nan, dtype=np.float32)
    valid_trate = np.isfinite(trc_stack) & np.isfinite(gls_stack) & (gls_stack > 0)
    trate_stack[valid_trate] = trc_stack[valid_trate] / gls_stack[valid_trate]
    response_data['Trate'] = trate_stack

    # Step 4: 像元级线性回归
    print("\n[Step 4] 像元级线性回归: Response ~ ΔSOS...")

    output_dir = OUTPUT_DIR / "Section_3.2_Phenology_Impact"
    output_dir.mkdir(parents=True, exist_ok=True)

    for resp_var in response_vars:
        print(f"\n  分析: {resp_var} ~ ΔSOS")

        slope_map, pvalue_map, r_squared_map = linear_regression_maps(
            delta_sos_stack, response_data[resp_var], min_frac=0.6
        )

        if mask is not None:
            slope_map[~mask] = NODATA_OUT
            pvalue_map[~mask] = NODATA_OUT
            r_squared_map[~mask] = NODATA_OUT

        # 保存结果
        write_geotiff(output_dir / f"{resp_var}_vs_deltaSOS_slope.tif", slope_map, profile)
        write_geotiff(output_dir / f"{resp_var}_vs_deltaSOS_pvalue.tif", pvalue_map, profile)
        write_geotiff(output_dir / f"{resp_var}_vs_deltaSOS_R2.tif", r_squared_map, profile)

        # 统计显著像元（p < 0.05）
        sig_mask = (pvalue_map < 0.05) & (pvalue_map != NODATA_OUT)
        n_sig = np.sum(sig_mask)
        mean_slope = np.nanmean(slope_map[sig_mask]) if n_sig > 0 else np.nan
        std_slope = np.nanstd(slope_map[sig_mask]) if n_sig > 0 else np.nan

        print(f"    显著像元数 (p<0.05): {n_sig}")
        print(f"    平均斜率: {mean_slope:.3f} ± {std_slope:.3f} mm/day per day ΔSOS")

    print("\n  ✓ Section 3.2 分析完成")
    print(f"  输出目录: {output_dir}")

# ==================== Section 3.3: 驱动因子分析 ====================
def section_3_3_driver_analysis(mask):
    """
    Section 3.3: Drivers of TRproduct decrease with spring phenology change

    方法：
    1. 全时段（1982-2018）偏相关分析 + VIF过滤
       控制其余变量后计算 TR 与各因子的偏相关系数 R

    2. 15年滑动窗口偏相关演变
       - 计算每个窗口的偏相关系数
       - Theil-Sen趋势分析
       - Mann-Kendall显著性检验

    3. 主驱动因子识别（|R| > 0.1）

    Parameters:
    -----------
    mask : ndarray
        有效掩膜 (H, W)

    Returns:
    --------
    None (结果保存到文件)
    """
    print("\n" + "="*70)
    print("Section 3.3: TRproduct 驱动因子分析")
    print("="*70)

    years = list(range(YEAR_START, YEAR_END + 1))
    n_years = len(years)

    # 读取模板
    data_first, profile, _ = read_geotiff(TRC_DIR / f"TRc_{years[0]}.tif")
    height, width = data_first.shape

    # 定义变量（Wang 2025 Eq. 3）
    response_vars = ['TRc', 'TRpheno', 'TRproduct']
    # 注意：SIFspr/SIFsum 实际使用日尺度 GPP 的季节均值作为代理
    predictor_vars = ['LSP', 'SMroot', 'Ta', 'Rs', 'P', 'SIFspr', 'SIFsum']
    available_vars = []
    missing_vars = []
    for var in predictor_vars:
        if var == 'LSP':
            available_vars.append(var)
            continue
        if var in ('SIFspr', 'SIFsum'):
            if _has_gpp_files(years[0]) or _has_gpp_files(years[-1]):
                available_vars.append(var)
            else:
                missing_vars.append(var)
            continue
        if var not in DAILY_VAR_SPECS:
            missing_vars.append(var)
            continue
        if _has_daily_files(var, years[0]) or _has_daily_files(var, years[-1]):
            available_vars.append(var)
        else:
            missing_vars.append(var)

    if missing_vars:
        print(f"  ⚠ 缺少日尺度输入，跳过变量: {', '.join(missing_vars)}")
    predictor_vars = available_vars
    if not predictor_vars:
        print("  ✗ 没有可用的预测变量，跳过 Section 3.3")
        return
    blocks = list(iter_valid_blocks(mask, BLOCK_SIZE_3_3))
    if not blocks:
        print("  ✗ 有效像元为空，跳过 Section 3.3")
        return

    # ========== Part 1: 全时段偏相关归因分析（带VIF过滤）==========
    print("\n[Part 1] 全时段偏相关归因分析（1982-2018，VIF过滤）...")

    # 预计算所有预测变量
    print("\n  预计算预测变量...")
    X_all_years = {var: [] for var in predictor_vars}

    for year in tqdm(years, desc="预计算自变量"):
        # 读取SOS和POS
        sos_data, _, sos_nodata = read_geotiff(PHENO_DIR / f"sos_t_{year}.tif")
        pos_data, _, pos_nodata = read_geotiff(PHENO_DIR / f"pos_doy_t_{year}.tif")

        sos_map = np.where(_is_valid_value(sos_data, sos_nodata), sos_data, np.nan)
        pos_map = np.where(_is_valid_value(pos_data, pos_nodata), pos_data, np.nan)

        if 'LSP' in predictor_vars:
            lsp_map = np.where((~np.isnan(sos_map)) & (~np.isnan(pos_map)) & (pos_map >= sos_map),
                               pos_map - sos_map + 1, np.nan)
            X_all_years['LSP'].append(lsp_map)

        # 计算LSP期间气象变量平均值
        for var in ['SMroot', 'Ta', 'Rs', 'P']:
            if var not in predictor_vars:
                continue
            lsp_avg = calculate_lsp_period_average(var, year, sos_map, pos_map)
            X_all_years[var].append(lsp_avg)

        if 'SIFspr' in predictor_vars:
            gpp_spring = calculate_seasonal_gpp(year, season='spring')
            if gpp_spring is None:
                gpp_spring = np.full((height, width), np.nan)
            X_all_years['SIFspr'].append(gpp_spring)

        if 'SIFsum' in predictor_vars:
            gpp_summer = calculate_seasonal_gpp(year, season='summer')
            if gpp_summer is None:
                gpp_summer = np.full((height, width), np.nan)
            X_all_years['SIFsum'].append(gpp_summer)

    # 转换为numpy数组
    for var in predictor_vars:
        X_all_years[var] = np.stack(X_all_years[var], axis=0)  # (n_years, H, W)

    # 对每个响应变量进行归因分析
    for response_var in response_vars:
        print(f"\n  响应变量: {response_var}")

        # 读取响应变量
        Y_stack = []
        for year in years:
            if response_var == 'TRc':
                data, _, nodata = read_geotiff(TRC_DIR / f"TRc_{year}.tif")
            else:
                data, _, nodata = read_geotiff(DECOMP_DIR / f"{response_var}_{year}.tif")

            Y_stack.append(np.where(_is_valid_value(data, nodata), data, np.nan))

        Y_stack = np.stack(Y_stack, axis=0)  # (n_years, H, W)

        # 逐像元VIF过滤 + 偏相关（分块）
        partial_r_maps = {var: np.full((height, width), NODATA_OUT, dtype=np.float32)
                          for var in predictor_vars}
        partial_p_maps = {var: np.full((height, width), NODATA_OUT, dtype=np.float32)
                          for var in predictor_vars}
        vif_filtered_vars = {var: np.zeros((height, width), dtype=bool)
                            for var in predictor_vars}
        r_squared_map = np.full((height, width), NODATA_OUT, dtype=np.float32)
        p = len(predictor_vars)
        min_rows = max(15, p + 8)
        if USE_BLOCK_PARALLEL:
            with ProcessPoolExecutor(max_workers=MAX_WORKERS_3_3) as executor:
                futures = []
                for r0, r1, c0, c1, block_mask in blocks:
                    X_block = {var: X_all_years[var][:, r0:r1, c0:c1] for var in predictor_vars}
                    Y_block = Y_stack[:, r0:r1, c0:c1]
                    args = (r0, r1, c0, c1, block_mask, predictor_vars, X_block, Y_block, min_rows)
                    futures.append(executor.submit(_partial_corr_block_worker, args))

                for fut in tqdm(as_completed(futures), total=len(futures),
                                desc=f"归因分析 {response_var}", leave=False):
                    r0, r1, c0, c1, pr_block, pp_block, vif_block, r2_block = fut.result()
                    for var in predictor_vars:
                        partial_r_maps[var][r0:r1, c0:c1] = pr_block[var]
                        partial_p_maps[var][r0:r1, c0:c1] = pp_block[var]
                        vif_filtered_vars[var][r0:r1, c0:c1] = vif_block[var]
                    r_squared_map[r0:r1, c0:c1] = r2_block
        else:
            for r0, r1, c0, c1, block_mask in tqdm(blocks, desc=f"归因分析 {response_var}", leave=False):
                X_block = {var: X_all_years[var][:, r0:r1, c0:c1] for var in predictor_vars}
                Y_block = Y_stack[:, r0:r1, c0:c1]
                args = (r0, r1, c0, c1, block_mask, predictor_vars, X_block, Y_block, min_rows)
                r0, r1, c0, c1, pr_block, pp_block, vif_block, r2_block = _partial_corr_block_worker(args)
                for var in predictor_vars:
                    partial_r_maps[var][r0:r1, c0:c1] = pr_block[var]
                    partial_p_maps[var][r0:r1, c0:c1] = pp_block[var]
                    vif_filtered_vars[var][r0:r1, c0:c1] = vif_block[var]
                r_squared_map[r0:r1, c0:c1] = r2_block

        # 保存全时段归因结果
        output_dir_full = OUTPUT_DIR / "Section_3.3_Drivers" / "Full_Period" / response_var
        output_dir_full.mkdir(parents=True, exist_ok=True)

        for var in predictor_vars:
            write_geotiff(output_dir_full / f"partial_r_{var}.tif",
                          partial_r_maps[var], profile)
            write_geotiff(output_dir_full / f"partial_p_{var}.tif",
                          partial_p_maps[var], profile)
            vif_out = vif_filtered_vars[var].astype(np.float32)
            vif_out[~mask] = NODATA_OUT
            write_geotiff(output_dir_full / f"vif_retained_{var}.tif",
                          vif_out, profile)

        write_geotiff(output_dir_full / f"R_squared.tif", r_squared_map, profile)

        # 统计主驱动因子（|β| > 0.1）
        print(f"\n    主驱动因子（|R| > 0.1）统计:")
        for var in predictor_vars:
            attr = partial_r_maps[var]
            main_driver_mask = (np.abs(attr) > 0.1) & (attr != NODATA_OUT)
            n_main = np.sum(main_driver_mask)
            print(f"      {var}: {n_main} 像元")

    # ========== Part 2: 15年滑动窗口偏相关演变 ==========
    print("\n[Part 2] 15年滑动窗口偏相关演变...")

    window_size = 15
    n_windows = n_years - window_size + 1

    if n_windows < 2:
        print("  ⚠ 警告：窗口数不足，跳过滑动窗口分析")
        return

    # 对TRproduct进行滑动窗口分析（论文重点关注）
    response_var = 'TRproduct'
    print(f"\n  响应变量: {response_var}")

    # 读取响应变量
    Y_stack = []
    for year in years:
        data, _, nodata = read_geotiff(DECOMP_DIR / f"{response_var}_{year}.tif")
        Y_stack.append(np.where(_is_valid_value(data, nodata), data, np.nan))

    Y_stack = np.stack(Y_stack, axis=0)  # (n_years, H, W)

    # 存储每个窗口的偏相关系数
    partial_r_evolution = {var: np.full((n_windows, height, width), NODATA_OUT, dtype=np.float32)
                           for var in predictor_vars}

    for win_idx in tqdm(range(n_windows), desc="滑动窗口"):
        start_idx = win_idx
        end_idx = win_idx + window_size
        start_year = years[start_idx]
        end_year = years[end_idx - 1]

        # 提取窗口数据
        Y_window = Y_stack[start_idx:end_idx]  # (window_size, H, W)
        X_window = {var: X_all_years[var][start_idx:end_idx] for var in predictor_vars}

        # 逐像元回归（分块/并行）
        p = len(predictor_vars)
        min_rows = max(10, p + 5)
        if USE_BLOCK_PARALLEL:
            with ProcessPoolExecutor(max_workers=MAX_WORKERS_3_3) as executor:
                futures = []
                for r0, r1, c0, c1, block_mask in blocks:
                    X_block = {var: X_window[var][:, r0:r1, c0:c1] for var in predictor_vars}
                    Y_block = Y_window[:, r0:r1, c0:c1]
                    args = (r0, r1, c0, c1, block_mask, predictor_vars, X_block, Y_block, min_rows)
                    futures.append(executor.submit(_partial_corr_window_block_worker, args))

                for fut in tqdm(as_completed(futures), total=len(futures),
                                desc="滑动窗口像元", leave=False):
                    r0, r1, c0, c1, pr_block = fut.result()
                    for var in predictor_vars:
                        partial_r_evolution[var][win_idx, r0:r1, c0:c1] = pr_block[var]
        else:
            for r0, r1, c0, c1, block_mask in tqdm(blocks, desc="滑动窗口像元", leave=False):
                X_block = {var: X_window[var][:, r0:r1, c0:c1] for var in predictor_vars}
                Y_block = Y_window[:, r0:r1, c0:c1]
                args = (r0, r1, c0, c1, block_mask, predictor_vars, X_block, Y_block, min_rows)
                r0, r1, c0, c1, pr_block = _partial_corr_window_block_worker(args)
                for var in predictor_vars:
                    partial_r_evolution[var][win_idx, r0:r1, c0:c1] = pr_block[var]

    # 保存每个窗口的偏相关系数
    output_dir_window = OUTPUT_DIR / "Section_3.3_Drivers" / "Moving_Window" / response_var
    output_dir_window.mkdir(parents=True, exist_ok=True)

    for win_idx in range(n_windows):
        start_year = years[win_idx]
        end_year = years[win_idx + window_size - 1]

        for var in predictor_vars:
            filename = f"partial_r_{var}_{start_year}-{end_year}.tif"
            write_geotiff(output_dir_window / filename,
                         partial_r_evolution[var][win_idx], profile)

    # ========== Part 3: 偏相关趋势分析 ==========
    print("\n[Part 3] 偏相关趋势分析（Theil-Sen + Mann-Kendall）...")

    output_dir_trend = OUTPUT_DIR / "Section_3.3_Drivers" / "Sensitivity_Trends" / response_var
    output_dir_trend.mkdir(parents=True, exist_ok=True)

    window_years = np.array(range(n_windows), dtype=np.float32)

    for var in predictor_vars:
        print(f"\n  分析: {var} 偏相关趋势")

        trend_slope_map = np.full((height, width), NODATA_OUT, dtype=np.float32)
        trend_pvalue_map = np.full((height, width), NODATA_OUT, dtype=np.float32)

        if USE_BLOCK_PARALLEL:
            with ProcessPoolExecutor(max_workers=MAX_WORKERS_3_3) as executor:
                futures = []
                for r0, r1, c0, c1, block_mask in blocks:
                    sens_block = partial_r_evolution[var][:, r0:r1, c0:c1]
                    args = (r0, r1, c0, c1, block_mask, sens_block, window_years, 0.6)
                    futures.append(executor.submit(_trend_block_worker, args))

                for fut in tqdm(as_completed(futures), total=len(futures),
                                desc=f"{var}趋势", leave=False):
                    r0, r1, c0, c1, slope_block, p_block = fut.result()
                    trend_slope_map[r0:r1, c0:c1] = slope_block
                    trend_pvalue_map[r0:r1, c0:c1] = p_block
        else:
            for r0, r1, c0, c1, block_mask in tqdm(blocks, desc=f"{var}趋势", leave=False):
                sens_block = partial_r_evolution[var][:, r0:r1, c0:c1]
                args = (r0, r1, c0, c1, block_mask, sens_block, window_years, 0.6)
                r0, r1, c0, c1, slope_block, p_block = _trend_block_worker(args)
                trend_slope_map[r0:r1, c0:c1] = slope_block
                trend_pvalue_map[r0:r1, c0:c1] = p_block

        # 保存趋势结果
        write_geotiff(output_dir_trend / f"{var}_trend_slope.tif", trend_slope_map, profile)
        write_geotiff(output_dir_trend / f"{var}_trend_pvalue.tif", trend_pvalue_map, profile)

        # 统计显著趋势（p < 0.05）
        sig_mask = (trend_pvalue_map < 0.05) & (trend_pvalue_map != NODATA_OUT)
        n_sig = np.sum(sig_mask)
        mean_trend = np.nanmean(trend_slope_map[sig_mask]) if n_sig > 0 else np.nan

        print(f"    显著趋势像元数 (p<0.05): {n_sig}")
        print(f"    平均趋势斜率: {mean_trend:.6f} per window")

    print("\n  ✓ Section 3.3 分析完成")
    print(f"  输出目录: {OUTPUT_DIR / 'Section_3.3_Drivers'}")

# ==================== 主程序 ====================
def main():
    print("\n" + "="*70)
    print("统计分析模块 - Wang (2025) Sections 3.2 & 3.3")
    print("="*70)

    # 读取掩膜
    print("\n[0] 读取掩膜...")
    mask_file = ANALYSIS_DIR / "masks" / "combined_mask.tif"

    if not mask_file.exists():
        print(f"  ⚠ 警告：掩膜文件不存在: {mask_file}")
        print("  使用默认掩膜（全部有效）")
        data_sample, profile_sample, _ = read_geotiff(TRC_DIR / f"TRc_{YEAR_START}.tif")
        mask = np.ones(data_sample.shape, dtype=bool)
    else:
        mask_data, profile, mask_nodata = read_geotiff(mask_file)
        mask = _is_valid_value(mask_data, mask_nodata) & (mask_data > 0)
        print(f"  ✓ 掩膜读取成功，有效像元数: {np.sum(mask)}")

    # Section 3.2: TRc组分 vs dSOS 回归分析
    try:
        section_3_2_phenology_impact(mask)
    except Exception as e:
        print(f"\n  ✗ Section 3.2 分析失败: {str(e)}")
        import traceback
        traceback.print_exc()

    # Section 3.3: TRproduct驱动因子分析
    try:
        section_3_3_driver_analysis(mask)
    except Exception as e:
        print(f"\n  ✗ Section 3.3 分析失败: {str(e)}")
        import traceback
        traceback.print_exc()

    print("\n" + "="*70)
    print("✓ 统计分析模块执行完成！")
    print(f"输出目录: {OUTPUT_DIR}")
    print("\n输出文件结构：")
    print("  ├── Section_3.2_Phenology_Impact/")
    print("  │   ├── TRc_vs_deltaSOS_slope.tif")
    print("  │   ├── TRc_vs_deltaSOS_pvalue.tif")
    print("  │   ├── TRc_vs_deltaSOS_R2.tif")
    print("  │   ├── TRpheno_vs_deltaSOS_slope.tif")
    print("  │   ├── TRpheno_vs_deltaSOS_pvalue.tif")
    print("  │   ├── TRpheno_vs_deltaSOS_R2.tif")
    print("  │   ├── TRproduct_vs_deltaSOS_slope.tif")
    print("  │   ├── TRproduct_vs_deltaSOS_pvalue.tif")
    print("  │   ├── TRproduct_vs_deltaSOS_R2.tif")
    print("  │   ├── Trate_vs_deltaSOS_slope.tif")
    print("  │   ├── Trate_vs_deltaSOS_pvalue.tif")
    print("  │   └── Trate_vs_deltaSOS_R2.tif")
    print("  ├── Section_3.3_Drivers/")
    print("  │   ├── Full_Period/")
    print("  │   │   ├── TRc/")
    print("  │   │   ├── TRpheno/")
    print("  │   │   └── TRproduct/")
    print("  │   │       ├── partial_r_{var}.tif (LSP, SMroot, Ta, Rs, P, SIFspr, SIFsum; SIF=GPP proxy)")
    print("  │   │       ├── partial_p_{var}.tif")
    print("  │   │       ├── vif_retained_{var}.tif")
    print("  │   │       └── R_squared.tif")
    print("  │   ├── Moving_Window/")
    print("  │   │   └── TRproduct/")
    print("  │   │       └── partial_r_{var}_{start}-{end}.tif")
    print("  │   └── Sensitivity_Trends/")
    print("  │       └── TRproduct/")
    print("  │           ├── {var}_trend_slope.tif")
    print("  │           └── {var}_trend_pvalue.tif")
    print("="*70)

if __name__ == "__main__":
    main()
