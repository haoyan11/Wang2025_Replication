#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module 04B: 统计分析 - Timing/Shape Decomposition (Sections 3.2 & 3.3)

Section 3.2: Timing/Shape/SOS/POS vs ΔSOS
    - 像元级线性回归

Section 3.3: Drivers of TRshape decrease with spring phenology change
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
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
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

# 导入配置
from _config import (
    ROOT, OUTPUT_ROOT, PHENO_DIR, TRC_ANNUAL_DIR, DECOMPOSITION_TIMING_DIR,
    GPP_DAILY_DIR, STATISTICAL_TIMING_DIR, YEAR_START, YEAR_END, NODATA_OUT
)

# 确保输出目录存在
STATISTICAL_TIMING_DIR.mkdir(parents=True, exist_ok=True)

# 向后兼容：保留旧变量名
ANALYSIS_DIR = OUTPUT_ROOT
DECOMP_DIR = DECOMPOSITION_TIMING_DIR
TRC_DIR = TRC_ANNUAL_DIR
OUTPUT_DIR = STATISTICAL_TIMING_DIR

# 脚本特定配置
GPP_DAILY_FORMAT = "GPP_{date}.tif"  # {date} = YYYYMMDD
METEO_DIR = ROOT / "Meteorological Data"
NODATA_ABS_MAX = 1e20
BLOCK_SIZE_3_3 = 128  # 优化：增大块减少进程创建次数（原64->128）

# ==================== 并行配置（性能优化）====================
# 性能优化策略（2-3x加速）：
#   1. BLOCK_SIZE增大（128）减少块数量和进程创建开销
#   2. MAX_WORKERS提升（8）充分利用多核CPU
#   3. MAX_IO_WORKERS提升（16）加速文件读取
#
# Windows注意事项：
#   - 如遇内存不足：降低MAX_WORKERS_3_3到4
#   - 如遇进程启动慢：设USE_BLOCK_PARALLEL = False
#   - WSL/Linux推荐（fork模式比spawn快）
# ==================================================================
USE_BLOCK_PARALLEL = True  # 块级并行（推荐保持开启）
MAX_WORKERS_3_3 = min(8, os.cpu_count() or 1)  # 优化：从4提升到8

# ==================== 偏相关计算优化配置 ====================
# 选择偏相关计算模式：
#   True  - 批量向量化模式（性能优先，10-50x加速，无VIF过滤）
#   False - 原版VIF过滤模式（统计严谨性优先，逐像元VIF过滤）
#
# 建议：
#   - 对于探索性分析和大规模数据，使用 True（快速）
#   - 对于最终发表结果，使用 False（统计严谨）
USE_BATCH_VECTORIZED = True

# ==================== I/O优化配置 ====================
# 1. 多线程并行读取（治本方案）
#    - I/O密集型任务，线程池能显著加速文件读取
#    - Windows/Linux均有效，无spawn模式序列化开销
MAX_IO_WORKERS = min(16, os.cpu_count() or 1)  # 优化：从8提升到16加速文件读取

# 2. 缓存配置（辅助方案）
#    - LSP/GPP计算需要读取大量小文件（每年数百个日尺度文件）
#    - 启用缓存可显著加速重复运行（首次运行会慢，后续快）
#
# 注意事项：
#   - 缓存文件存储在 ANALYSIS_DIR / "Cache_Statistical_TimingShape"
#   - 如果输入数据更新，需手动删除缓存目录
#   - 缓存文件较大（每个变量/年约几MB），注意磁盘空间
USE_LSP_CACHE = False   # 启用LSP期间气象变量均值缓存（推荐开启）
USE_GPP_CACHE = False   # 启用季节GPP均值缓存（推荐开启）
CACHE_DIR = ANALYSIS_DIR / "Cache_Statistical_TimingShape"
if USE_LSP_CACHE or USE_GPP_CACHE:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

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


def linear_regression_maps(x_stack, y_stack, min_frac=0.6):
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

    with np.errstate(divide="ignore", invalid="ignore"):
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

    # 优先使用Numba加速版本（约5-20x加速）
    if NUMBA_AVAILABLE:
        if not np.isfinite(y_std).all() or not np.isfinite(X_std).all():
            return None, None
        r, p_vals = partial_corr_from_std_numba(y_std, X_std)
        return r.astype(np.float32), p_vals.astype(np.float32)

    # 原版NumPy/SciPy实现
    data = np.column_stack([y_std, X_std])
    if not np.isfinite(data).all():
        return None, None

    try:
        corr = np.corrcoef(data, rowvar=False)
        prec = np.linalg.inv(corr)
    except np.linalg.LinAlgError:
        return None, None

    denom = np.sqrt(prec[0, 0] * np.diag(prec)[1:])
    with np.errstate(divide="ignore", invalid="ignore"):
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


def filter_statistical_outliers(slope_map, r2_map=None, pvalue_map=None, mask=None,
                                 slope_threshold=100.0, print_details=True):
    """
    过滤统计结果异常值并打印详细信息

    Parameters:
    -----------
    slope_map : ndarray (H, W)
        回归斜率地图
    r2_map : ndarray (H, W), optional
        R²地图
    pvalue_map : ndarray (H, W), optional
        p值地图
    mask : ndarray (H, W), optional
        空间掩膜
    slope_threshold : float
        斜率绝对值阈值（默认100）
    print_details : bool
        是否打印详细异常值信息（默认True）

    Returns:
    --------
    n_filtered : int
        被过滤的像元数量
    valid_mask : ndarray (H, W)
        过滤后的有效像元掩膜
    outlier_info : dict
        异常值详细信息字典
    """
    # 确定有效像元总数
    if mask is not None:
        base_mask = mask
        n_total = np.sum(mask)
    else:
        base_mask = np.ones_like(slope_map, dtype=bool)
        n_total = slope_map.size

    # 初始化异常值统计
    outlier_info = {
        'n_total': n_total,
        'n_valid': 0,
        'n_filtered': 0,
        'pct_filtered': 0.0,
        'slope_nan_inf': 0,
        'slope_extreme': 0,
        'r2_invalid': 0,
        'pvalue_invalid': 0,
        'outlier_values': {}
    }

    # 从base_mask开始，逐步应用过滤条件
    valid = base_mask.copy()

    # 1. 斜率异常值过滤
    # 1a. NaN/Inf值
    slope_finite = np.isfinite(slope_map)
    slope_nan_inf_mask = base_mask & ~slope_finite
    outlier_info['slope_nan_inf'] = np.sum(slope_nan_inf_mask)

    # 1b. 极端值
    slope_extreme_mask = base_mask & slope_finite & (np.abs(slope_map) >= slope_threshold)
    outlier_info['slope_extreme'] = np.sum(slope_extreme_mask)
    if outlier_info['slope_extreme'] > 0:
        extreme_values = slope_map[slope_extreme_mask]
        outlier_info['outlier_values']['slope_range'] = (
            f"[{np.min(extreme_values):.2f}, {np.max(extreme_values):.2f}]"
        )

    # 应用斜率过滤
    valid &= slope_finite & (np.abs(slope_map) < slope_threshold)

    # 2. R²范围检查（0到1）
    if r2_map is not None:
        r2_invalid_mask = base_mask & valid & (
            ~np.isfinite(r2_map) | (r2_map < 0) | (r2_map > 1)
        )
        outlier_info['r2_invalid'] = np.sum(r2_invalid_mask)
        if outlier_info['r2_invalid'] > 0:
            invalid_r2 = r2_map[r2_invalid_mask]
            outlier_info['outlier_values']['r2_range'] = (
                f"[{np.nanmin(invalid_r2):.3f}, {np.nanmax(invalid_r2):.3f}]"
            )
        valid &= np.isfinite(r2_map) & (r2_map >= 0) & (r2_map <= 1)

    # 3. p值范围检查（0到1）
    if pvalue_map is not None:
        pvalue_invalid_mask = base_mask & valid & (
            ~np.isfinite(pvalue_map) | (pvalue_map < 0) | (pvalue_map > 1)
        )
        outlier_info['pvalue_invalid'] = np.sum(pvalue_invalid_mask)
        if outlier_info['pvalue_invalid'] > 0:
            invalid_p = pvalue_map[pvalue_invalid_mask]
            outlier_info['outlier_values']['pvalue_range'] = (
                f"[{np.nanmin(invalid_p):.3e}, {np.nanmax(invalid_p):.3e}]"
            )
        valid &= np.isfinite(pvalue_map) & (pvalue_map >= 0) & (pvalue_map <= 1)

    # 统计最终结果
    outlier_info['n_valid'] = np.sum(valid)
    outlier_info['n_filtered'] = n_total - outlier_info['n_valid']
    outlier_info['pct_filtered'] = outlier_info['n_filtered'] / n_total * 100 if n_total > 0 else 0

    # 打印详细信息
    if print_details and outlier_info['n_filtered'] > 0:
        print(f"\n    === 异常值过滤详情 ===")
        print(f"      总像元数: {n_total}")
        print(f"      有效像元: {outlier_info['n_valid']} ({100-outlier_info['pct_filtered']:.1f}%)")
        print(f"      过滤像元: {outlier_info['n_filtered']} ({outlier_info['pct_filtered']:.1f}%)")
        print(f"\n      过滤原因分解:")

        if outlier_info['slope_nan_inf'] > 0:
            print(f"        斜率NaN/Inf:   {outlier_info['slope_nan_inf']:>6} "
                  f"({outlier_info['slope_nan_inf']/n_total*100:>5.1f}%)")

        if outlier_info['slope_extreme'] > 0:
            print(f"        斜率极端值:    {outlier_info['slope_extreme']:>6} "
                  f"({outlier_info['slope_extreme']/n_total*100:>5.1f}%) "
                  f"[|slope| > {slope_threshold}]")
            if 'slope_range' in outlier_info['outlier_values']:
                print(f"          极端值范围: {outlier_info['outlier_values']['slope_range']}")

        if outlier_info['r2_invalid'] > 0:
            print(f"        R²异常值:      {outlier_info['r2_invalid']:>6} "
                  f"({outlier_info['r2_invalid']/n_total*100:>5.1f}%) "
                  f"[R²<0 或 R²>1 或 NaN]")
            if 'r2_range' in outlier_info['outlier_values']:
                print(f"          异常R²范围: {outlier_info['outlier_values']['r2_range']}")

        if outlier_info['pvalue_invalid'] > 0:
            print(f"        p值异常值:     {outlier_info['pvalue_invalid']:>6} "
                  f"({outlier_info['pvalue_invalid']/n_total*100:>5.1f}%) "
                  f"[p<0 或 p>1 或 NaN]")
            if 'pvalue_range' in outlier_info['outlier_values']:
                print(f"          异常p值范围: {outlier_info['outlier_values']['pvalue_range']}")

        print()

    elif print_details and outlier_info['n_filtered'] == 0:
        print(f"    ✓ 无异常值（全部{n_total}个像元均有效）")

    return outlier_info['n_filtered'], valid, outlier_info


def print_comprehensive_statistics(data_map, pvalue_map=None, mask=None,
                                   var_name="", unit="",
                                   print_percentiles=True, print_sign_split=True):
    """
    打印全面统计（全部有效像元 + 显著性像元）

    Parameters:
    -----------
    data_map : ndarray (H, W)
        统计量地图（如slope, R, R²）
    pvalue_map : ndarray (H, W), optional
        p值地图（用于显著性筛选），如为None则不区分显著性
    mask : ndarray (H, W), optional
        有效像元掩膜
    var_name : str
        变量名
    unit : str
        单位
    print_percentiles : bool
        是否打印百分位数（默认True）
    print_sign_split : bool
        是否打印正负值统计（默认True）

    Returns:
    --------
    stats_dict : dict
        包含统计量的字典
    """
    # 提取所有有效像元
    if mask is not None:
        valid_all = mask & np.isfinite(data_map)
    else:
        valid_all = np.isfinite(data_map)

    values_all = data_map[valid_all]
    n_all = len(values_all)

    if n_all == 0:
        print(f"\n  === {var_name} 统计 ===")
        print(f"    ⚠️ 无有效数据")
        return {}

    # 计算全部有效像元统计量
    stats_all = {
        'mean': np.mean(values_all),
        'std': np.std(values_all),
        'median': np.median(values_all),
        'n': n_all
    }

    if print_percentiles:
        stats_all.update({
            'p5': np.percentile(values_all, 5),
            'p25': np.percentile(values_all, 25),
            'p75': np.percentile(values_all, 75),
            'p95': np.percentile(values_all, 95)
        })

    if print_sign_split:
        n_pos = np.sum(values_all > 0)
        n_neg = np.sum(values_all < 0)
        n_zero = np.sum(values_all == 0)
        stats_all.update({
            'n_pos': n_pos,
            'n_neg': n_neg,
            'n_zero': n_zero,
            'pct_pos': n_pos / n_all * 100 if n_all > 0 else 0,
            'pct_neg': n_neg / n_all * 100 if n_all > 0 else 0
        })

    # 打印全部有效像元统计
    print(f"\n  === {var_name} 统计 ===")
    print(f"\n    全部有效像元（N={n_all}）:")
    print(f"      平均值: {stats_all['mean']:+.4f} ± {stats_all['std']:.4f} {unit}")
    print(f"      中位数: {stats_all['median']:+.4f} {unit}")

    if print_percentiles:
        print(f"      百分位: [P5: {stats_all['p5']:+.4f}, P25: {stats_all['p25']:+.4f}, "
              f"P75: {stats_all['p75']:+.4f}, P95: {stats_all['p95']:+.4f}] {unit}")

    if print_sign_split:
        print(f"      正值像元: {stats_all['n_pos']} ({stats_all['pct_pos']:.1f}%)")
        print(f"      负值像元: {stats_all['n_neg']} ({stats_all['pct_neg']:.1f}%)")

    # 如果提供了p值地图，计算显著性像元统计
    stats_sig = None
    if pvalue_map is not None:
        valid_sig = valid_all & (pvalue_map < 0.05)
        values_sig = data_map[valid_sig]
        n_sig = len(values_sig)

        if n_sig > 0:
            sig_pct = n_sig / n_all * 100

            stats_sig = {
                'mean': np.mean(values_sig),
                'std': np.std(values_sig),
                'median': np.median(values_sig),
                'n': n_sig,
                'pct': sig_pct
            }

            if print_percentiles:
                stats_sig.update({
                    'p5': np.percentile(values_sig, 5),
                    'p25': np.percentile(values_sig, 25),
                    'p75': np.percentile(values_sig, 75),
                    'p95': np.percentile(values_sig, 95)
                })

            # 打印显著性像元统计
            print(f"\n    显著性像元（p<0.05, N={n_sig}, {sig_pct:.1f}%）:")
            print(f"      平均值: {stats_sig['mean']:+.4f} ± {stats_sig['std']:.4f} {unit}")
            print(f"      中位数: {stats_sig['median']:+.4f} {unit}")

            if print_percentiles:
                print(f"      百分位: [P5: {stats_sig['p5']:+.4f}, P25: {stats_sig['p25']:+.4f}, "
                      f"P75: {stats_sig['p75']:+.4f}, P95: {stats_sig['p95']:+.4f}] {unit}")

    return {'all': stats_all, 'significant': stats_sig}


# ==================== Numba加速版本（性能优化）====================
if NUMBA_AVAILABLE:
    @jit(nopython=True, cache=True)
    def _lstsq_simple_numba(X, y):
        """Numba优化的简单最小二乘求解（用于VIF计算）"""
        XtX = X.T @ X
        Xty = X.T @ y
        try:
            beta = np.linalg.solve(XtX, Xty)
            return beta, True
        except:
            return np.zeros(X.shape[1], dtype=X.dtype), False

    @jit(nopython=True, cache=True)
    def calculate_vif_numba(X):
        """Numba加速的VIF计算"""
        n_samples, n_features = X.shape
        vif = np.zeros(n_features, dtype=X.dtype)
        for i in range(n_features):
            y = X[:, i]
            X_others = np.empty((n_samples, n_features - 1), dtype=X.dtype)
            col = 0
            for j in range(n_features):
                if j != i:
                    X_others[:, col] = X[:, j]
                    col += 1
            beta, success = _lstsq_simple_numba(X_others, y)
            if not success:
                vif[i] = np.inf
                continue
            y_pred = X_others @ beta
            ss_res = np.sum((y - y_pred) ** 2)
            ss_tot = np.sum((y - np.mean(y)) ** 2)
            if ss_tot > 0:
                r_squared = 1.0 - (ss_res / ss_tot)
                if r_squared < 0.9999:
                    vif[i] = 1.0 / (1.0 - r_squared)
                else:
                    vif[i] = np.inf
            else:
                vif[i] = np.inf
        return vif

    @jit(nopython=True, cache=True)
    def partial_corr_from_std_numba(y_std, X_std):
        """Numba加速的偏相关系数计算（基于精度矩阵）"""
        n = len(y_std)
        p = X_std.shape[1]
        data = np.empty((n, p + 1), dtype=y_std.dtype)
        data[:, 0] = y_std
        data[:, 1:] = X_std
        corr = np.zeros((p + 1, p + 1), dtype=y_std.dtype)
        for i in range(p + 1):
            for j in range(i, p + 1):
                c = np.mean(data[:, i] * data[:, j])
                corr[i, j] = c
                corr[j, i] = c
        try:
            prec = np.linalg.inv(corr)
        except:
            return np.full(p, np.nan, dtype=y_std.dtype), np.full(p, np.nan, dtype=y_std.dtype)
        r = np.empty(p, dtype=y_std.dtype)
        for i in range(p):
            denom = np.sqrt(prec[0, 0] * prec[i + 1, i + 1])
            if denom > 0:
                r[i] = -prec[0, i + 1] / denom
            else:
                r[i] = np.nan
        r = np.clip(r, -0.999999, 0.999999)
        df = n - p - 1
        p_vals = np.empty(p, dtype=y_std.dtype)
        if df > 0:
            for i in range(p):
                if not np.isnan(r[i]):
                    t_stat = r[i] * np.sqrt(df / (1.0 - r[i] ** 2))
                    abs_t = abs(t_stat)
                    if abs_t < 10:
                        p_vals[i] = 2.0 * (1.0 - 0.5 * (1.0 + np.tanh(abs_t * np.sqrt(2.0 / np.pi))))
                    else:
                        p_vals[i] = 0.0
                else:
                    p_vals[i] = np.nan
        else:
            p_vals[:] = np.nan
        return r, p_vals

# ==================== 批量向量化版本（超高性能，参考用户代码）====================
def partial_corr_batch_vectorized(Y_block, X_block, predictor_vars, min_rows, enable_vif=False):
    """
    批量向量化偏相关计算（参考用户代码优化架构）

    性能优势：批量矩阵操作 >> 逐像元循环
    - 使用einsum批量计算协方差矩阵
    - 使用pinv批量求逆所有像元的矩阵
    - 向量化提取偏相关系数

    Parameters:
    -----------
    Y_block : ndarray, shape (n_years, block_h, block_w)
        响应变量块
    X_block : dict of ndarray, shape (n_years, block_h, block_w)
        预测变量块字典
    predictor_vars : list of str
        预测变量名列表
    min_rows : int
        最小有效样本数
    enable_vif : bool
        是否启用VIF过滤（默认False，因为VIF难以批量化）

    Returns:
    --------
    partial_r_block : dict of ndarray (block_h, block_w)
    partial_p_block : dict of ndarray (block_h, block_w)
    vif_block : dict of ndarray (block_h, block_w)
    r2_block : ndarray (block_h, block_w)
    """
    n_years, block_h, block_w = Y_block.shape
    N = block_h * block_w
    P = len(predictor_vars) + 1  # Y + X变量数

    # 初始化输出
    partial_r_block = {var: np.full((block_h, block_w), NODATA_OUT, dtype=np.float32)
                       for var in predictor_vars}
    partial_p_block = {var: np.full((block_h, block_w), NODATA_OUT, dtype=np.float32)
                       for var in predictor_vars}
    vif_block = {var: np.zeros((block_h, block_w), dtype=bool) for var in predictor_vars}
    r2_block = np.full((block_h, block_w), NODATA_OUT, dtype=np.float32)

    # 构建数据立方体 (N, T, P)
    cube = np.empty((N, n_years, P), dtype=np.float32)
    cube[:, :, 0] = Y_block.reshape(n_years, -1).T  # Y变量
    for j, var in enumerate(predictor_vars, start=1):
        cube[:, :, j] = X_block[var].reshape(n_years, -1).T  # X变量

    # 有效性掩膜 (N, T)
    valid_year = ~np.isnan(cube).any(axis=2)
    n_eff_all = valid_year.sum(axis=1).astype(np.int16)
    good_idx = np.where(n_eff_all >= min_rows)[0]

    if good_idx.size == 0:
        return partial_r_block, partial_p_block, vif_block, r2_block

    # 提取有效像元
    data = cube[good_idx]  # (Ng, T, P)
    vmask = valid_year[good_idx]  # (Ng, T)
    n_eff = vmask.sum(axis=1).astype(np.float32)  # (Ng,)

    # 批量标准化（Z-score）
    sums = np.sum(np.where(vmask[..., None], data, 0.0), axis=1, keepdims=True)  # (Ng, 1, P)
    means = sums / vmask.sum(axis=1, keepdims=True)[..., None]
    dmean = np.where(vmask[..., None], data - means, 0.0)  # (Ng, T, P)

    # 批量计算协方差矩阵（einsum核心优化！）
    denom = (n_eff[:, None, None] - 1.0)
    cov = np.einsum('ntk,ntl->nkl', dmean, dmean) / denom  # (Ng, P, P)

    # 批量矩阵求逆（pinv核心优化！）
    inv_cov = np.linalg.pinv(cov.astype(np.float64)).astype(np.float32)  # (Ng, P, P)

    # 向量化提取偏相关系数
    num = -inv_cov[:, 0, 1:]  # (Ng, P-1)
    den = np.sqrt(inv_cov[:, 0, 0, None] * inv_cov[:, 1:, 1:].diagonal(axis1=1, axis2=2))
    r_partial = np.divide(num, den, out=np.full_like(num, np.nan), where=den > 0)
    r_partial = np.clip(r_partial, -0.999999, 0.999999)

    # 向量化计算p值
    m = len(predictor_vars)  # 控制变量数
    df = (n_eff.astype(np.int32) - m - 2)[:, None]  # (Ng, 1)
    ok_df = (df > 3) & np.isfinite(r_partial)

    # 利用广播计算 t 统计量（避免布尔索引维度不匹配）
    t_stat_full = r_partial * np.sqrt(df / (1.0 - r_partial**2))
    t_stat = np.where(ok_df, t_stat_full, np.nan)

    # 计算p值（使用广播的df）
    p_vals_full = 2.0 * (1.0 - stats.t.cdf(np.abs(t_stat), df))
    p_vals = np.where(ok_df, p_vals_full, np.nan)

    # 计算R²（多元回归）- 批量向量化版本
    # 使用协方差矩阵直接计算R²（避免lstsq循环）
    #
    # 理论公式：R² = Var(ŷ) / Var(y)
    # 其中：
    #   β = Cov_XX^-1 @ Cov_yX  (回归系数)
    #   Var(ŷ) = Cov_yX^T @ Cov_XX^-1 @ Cov_yX  (预测值方差)
    #   Var(y) = cov[:, 0, 0]  (y的方差)
    #
    # 因此：R² = (Cov_yX^T @ Cov_XX^-1 @ Cov_yX) / Var(y)
    #
    # 关键修复：必须先提取Cov_XX子矩阵，再对它求逆
    #          不能从inv_cov中提取子矩阵（矩阵逆的子矩阵 ≠ 子矩阵的逆）

    # 提取协方差向量和矩阵
    Cov_yX = cov[:, 0, 1:]  # (Ng, m) - y与X的协方差向量
    Cov_XX = cov[:, 1:, 1:]  # (Ng, m, m) - X的协方差矩阵
    var_y_total = cov[:, 0, 0]  # (Ng,) - y的总方差

    # 关键步骤：对Cov_XX求逆（不是从inv_cov中提取！）
    Cov_XX_inv = np.linalg.pinv(Cov_XX.astype(np.float64)).astype(np.float32)  # (Ng, m, m)

    # 批量计算预测值方差：Var(ŷ) = Cov_yX^T @ Cov_XX^-1 @ Cov_yX
    var_y_pred = np.einsum('ni,nij,nj->n', Cov_yX, Cov_XX_inv, Cov_yX)  # (Ng,)

    # 计算R²：R² = Var(ŷ) / Var(y)
    r2_vals = var_y_pred / (var_y_total + 1e-10)
    r2_vals = np.clip(r2_vals, 0.0, 1.0)  # 确保R²在 [0, 1] 范围内

    # 样本数不足的像元设为nan
    r2_vals[n_eff < min_rows] = np.nan

    # 填充结果到块
    r_blk = {var: np.full((N,), np.nan, dtype=np.float32) for var in predictor_vars}
    p_blk = {var: np.full((N,), np.nan, dtype=np.float32) for var in predictor_vars}
    r2_blk = np.full((N,), np.nan, dtype=np.float32)

    for j, var in enumerate(predictor_vars):
        r_blk[var][good_idx] = r_partial[:, j]
        p_blk[var][good_idx] = p_vals[:, j]
        vif_block[var].reshape(-1)[good_idx] = True  # 标记为已处理

    r2_blk[good_idx] = r2_vals

    # Reshape回块形状
    for var in predictor_vars:
        partial_r_block[var] = r_blk[var].reshape(block_h, block_w)
        partial_p_block[var] = p_blk[var].reshape(block_h, block_w)
    r2_block = r2_blk.reshape(block_h, block_w)

    return partial_r_block, partial_p_block, vif_block, r2_block


def _partial_corr_block_worker_batch(args):
    """
    批量向量化版本的块处理器（包装器）

    与 _partial_corr_block_worker 接口兼容，但使用批量向量化计算
    性能提升：10-50x（相比原版），5-20x（相比Numba版）

    注意：此版本不执行VIF过滤（为了批量化性能）
    """
    (r0, r1, c0, c1, block_mask, predictor_vars, X_block, Y_block, min_rows) = args

    # 调用批量向量化函数
    partial_r_block, partial_p_block, vif_block, r2_block = \
        partial_corr_batch_vectorized(Y_block, X_block, predictor_vars, min_rows, enable_vif=False)

    return (r0, r1, c0, c1, partial_r_block, partial_p_block, vif_block, r2_block)


def _partial_corr_window_block_worker_batch(args):
    """
    批量向量化版本的滑动窗口块处理器

    与 _partial_corr_window_block_worker 接口兼容，使用批量向量化计算
    """
    (r0, r1, c0, c1, block_mask, predictor_vars, X_block, Y_block, min_rows) = args

    # 调用批量向量化函数（只需要partial_r）
    partial_r_block, _, _, _ = \
        partial_corr_batch_vectorized(Y_block, X_block, predictor_vars, min_rows, enable_vif=False)

    return r0, r1, c0, c1, partial_r_block


# ==================== 原版函数（兼容性）====================
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
    # 优先使用Numba加速版本
    if NUMBA_AVAILABLE:
        return calculate_vif_numba(X)

    # 原版NumPy实现
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

    优化策略：多线程并行读取日文件，显著加速I/O

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
    # 缓存检查（可选优化）
    if USE_GPP_CACHE:
        cache_file = CACHE_DIR / f"GPP_{season}_{year}.tif"
        if cache_file.exists():
            data, _, _ = read_geotiff(cache_file)
            return data

    if not _has_gpp_files(year):
        return None

    months = SPRING_MONTHS if season == 'spring' else SUMMER_MONTHS

    # ========== 多线程并行读取优化 ==========
    # 步骤1：收集所有需要的文件路径
    file_paths = []
    for month in months:
        for day in range(1, 32):
            try:
                date_str = datetime(year, month, day).strftime("%Y%m%d")
                gpp_file = GPP_DAILY_DIR / GPP_DAILY_FORMAT.format(date=date_str)
                if gpp_file.exists():
                    file_paths.append(gpp_file)
            except ValueError:
                continue

    if len(file_paths) == 0:
        return None

    # 步骤2：多线程并行读取所有文件
    def _read_and_process(file_path):
        """工作线程：读取单个文件并处理"""
        data, profile, nodata = read_geotiff(file_path)
        valid_data = np.where(_is_valid_value(data, nodata), data, np.nan)
        return valid_data, profile

    gpp_list = []
    profile_template = None

    with ThreadPoolExecutor(max_workers=MAX_IO_WORKERS) as executor:
        # 提交所有读取任务
        futures = [executor.submit(_read_and_process, fp) for fp in file_paths]

        # 收集结果（保持顺序）
        for future in futures:
            valid_data, profile = future.result()
            if profile_template is None:
                profile_template = profile
            gpp_list.append(valid_data)

    # 步骤3：内存计算
    gpp_stack = np.stack(gpp_list, axis=0)
    gpp_seasonal = np.nanmean(gpp_stack, axis=0)

    # 保存缓存（可选优化）
    if USE_GPP_CACHE and profile_template is not None:
        cache_file = CACHE_DIR / f"GPP_{season}_{year}.tif"
        if not cache_file.exists():
            write_geotiff(cache_file, gpp_seasonal, profile_template)

    return gpp_seasonal


def calculate_lsp_period_average(var_name, year, sos_map, pos_map):
    """
    计算LSP期间的变量平均值（Wang 2025 Section 2.2.2）

    优化策略：多线程并行读取日文件，显著加速I/O

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
    # 缓存检查（可选优化）
    # 注意：缓存键包含年份，但不包含SOS/POS（假设同一年的SOS/POS不变）
    if USE_LSP_CACHE:
        cache_file = CACHE_DIR / f"LSP_{var_name}_{year}.tif"
        if cache_file.exists():
            data, _, _ = read_geotiff(cache_file)
            return data

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
        (sos_int <= 366) & (pos_int <= 366)  # 允许闰年DOY=366
    )

    # Clip到365以匹配365天气候态（闰年DOY=366映射到平年DOY=365，即12月31日）
    sos_int = np.clip(sos_int, 1, 365)
    pos_int = np.clip(pos_int, 1, 365)

    if not np.any(valid):
        return np.full((height, width), np.nan, dtype=np.float32)

    doy_start = int(np.nanmin(sos_int[valid]))
    doy_end = int(np.nanmax(pos_int[valid]))
    doy_start = max(1, min(365, doy_start))
    doy_end = max(1, min(365, doy_end))

    # ========== 多线程并行读取优化 ==========
    # 步骤1：准备文件路径列表
    file_doy_pairs = []
    for doy in range(doy_start, doy_end + 1):
        date_str = (datetime(year, 1, 1) + timedelta(days=doy - 1)).strftime("%Y%m%d")
        var_file = var_dir / pattern.format(date=date_str)
        if var_file.exists():
            file_doy_pairs.append((var_file, doy))

    # 步骤2：多线程并行读取所有文件
    def _read_and_process_lsp(args):
        """工作线程：读取单个文件并处理"""
        file_path, doy = args
        data, profile, nodata = read_geotiff(file_path)
        valid_data = _is_valid_value(data, nodata)
        return doy, data.astype(np.float32), valid_data, profile

    daily_data_dict = {}  # {doy: (data_array, valid_mask)}
    profile_template = None

    with ThreadPoolExecutor(max_workers=MAX_IO_WORKERS) as executor:
        # 提交所有读取任务
        futures = [executor.submit(_read_and_process_lsp, pair) for pair in file_doy_pairs]

        # 收集结果
        for future in futures:
            doy, data, valid_data, profile = future.result()
            if profile_template is None:
                profile_template = profile
            daily_data_dict[doy] = (data, valid_data)

    # 步骤3：从内存字典中提取数据并计算（无文件I/O）
    total_sum = np.zeros((height, width), dtype=np.float32)
    total_cnt = np.zeros((height, width), dtype=np.int16)

    for doy in range(doy_start, doy_end + 1):
        if doy not in daily_data_dict:
            continue

        in_window = valid & (sos_int <= doy) & (pos_int >= doy)
        if not np.any(in_window):
            continue

        data, valid_data = daily_data_dict[doy]
        use_mask = in_window & valid_data

        if not np.any(use_mask):
            continue

        total_sum[use_mask] += data[use_mask]
        total_cnt[use_mask] += 1

    window_len = pos_int - sos_int + 1
    lsp_avg = np.full((height, width), np.nan, dtype=np.float32)
    good = valid & (total_cnt >= 0.6 * window_len)
    lsp_avg[good] = total_sum[good] / total_cnt[good]

    # 保存缓存（可选优化）
    if USE_LSP_CACHE and profile_template is not None:
        cache_file = CACHE_DIR / f"LSP_{var_name}_{year}.tif"
        if not cache_file.exists():
            write_geotiff(cache_file, lsp_avg, profile_template)

    return lsp_avg



def section_3_2_timing_shape(mask):
    print("\n" + "=" * 70)
    print("Section 3.2A: TRshape/TRtiming/TRsos/TRpos/Trate vs ΔSOS")
    print("=" * 70)

    years = list(range(YEAR_START, YEAR_END + 1))
    n_years = len(years)

    # Template (使用TRc以确保与04a/04c空间一致)
    data_first, profile, _ = read_geotiff(TRC_DIR / f"TRc_{years[0]}.tif")
    height, width = data_first.shape

    # SOS/POS stack and SOSav
    print("\n[Step 1] 读取SOS/POS并计算SOSav...")
    sos_stack = []
    pos_stack = []
    for year in tqdm(years, desc="读取SOS/POS"):
        sos_data, _, sos_nodata = read_geotiff(PHENO_DIR / f"sos_gpp_{year}.tif")
        pos_data, _, pos_nodata = read_geotiff(PHENO_DIR / f"pos_doy_gpp_{year}.tif")
        sos_valid = np.where(_is_valid_value(sos_data, sos_nodata), sos_data, np.nan)
        pos_valid = np.where(_is_valid_value(pos_data, pos_nodata), pos_data, np.nan)
        sos_stack.append(sos_valid)
        pos_stack.append(pos_valid)
    sos_stack = np.stack(sos_stack, axis=0)
    pos_stack = np.stack(pos_stack, axis=0)
    sos_av = np.nanmean(sos_stack, axis=0)

    # ΔSOS (advance < 0)
    print("\n[Step 2] 计算 ΔSOS = SOS_year - SOSav...")
    delta_sos_stack = sos_stack - sos_av

    # GLS for Trate
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

    # Response variables
    print("\n[Step 3] 读取响应变量...")
    # 核心分解组分
    decomp_vars = ["TRshape", "TRtiming", "TRsos", "TRpos"]
    response_data = {}

    for resp_var in decomp_vars:
        data_stack = []
        for year in tqdm(years, desc=f"读取{resp_var}", leave=False):
            data, _, nodata = read_geotiff(DECOMP_DIR / f"{resp_var}_{year}.tif")
            data_stack.append(np.where(_is_valid_value(data, nodata), data, np.nan))
        response_data[resp_var] = np.stack(data_stack, axis=0)

    # TRc（用于对比和计算Trate）
    print("\n  读取 TRc...")
    trc_stack = []
    for year in tqdm(years, desc="读取TRc", leave=False):
        trc_data, _, trc_nodata = read_geotiff(TRC_DIR / f"TRc_{year}.tif")
        trc_stack.append(np.where(_is_valid_value(trc_data, trc_nodata), trc_data, np.nan))
    trc_stack = np.stack(trc_stack, axis=0)
    response_data["TRc"] = trc_stack

    # Trate = TRc / GLS（变化窗口归一化）
    print("\n  计算 Trate = TRc / GLS（变化窗口）...")
    trate_stack = np.full_like(trc_stack, np.nan, dtype=np.float32)
    valid_trate = np.isfinite(trc_stack) & np.isfinite(gls_stack) & (gls_stack > 0)
    trate_stack[valid_trate] = trc_stack[valid_trate] / gls_stack[valid_trate]
    response_data["Trate"] = trate_stack

    # 完整的响应变量列表（按重要性排序）
    response_vars = ["TRc", "TRshape", "TRtiming", "TRsos", "TRpos", "Trate"]

    # Regression
    print("\n[Step 4] 像元级线性回归: Response ~ ΔSOS...")
    output_dir = OUTPUT_DIR / "Section_3.2_Phenology_Impact"
    output_dir.mkdir(parents=True, exist_ok=True)

    for resp_var in response_vars:
        print(f"\n  分析: {resp_var} ~ ΔSOS")
        slope_map, pvalue_map, r_squared_map = linear_regression_maps(
            delta_sos_stack, response_data[resp_var], min_frac=0.6
        )

        # 异常值过滤（斜率阈值根据变量调整）
        slope_threshold = 10.0  # TRshape/TRtiming/TRsos/TRpos的合理范围
        if resp_var == 'Trate':
            slope_threshold = 1.0  # Trate的斜率范围更小

        n_filtered, valid_mask, outlier_info = filter_statistical_outliers(
            slope_map, r2_map=r_squared_map, pvalue_map=pvalue_map,
            mask=mask, slope_threshold=slope_threshold, print_details=True
        )

        # 应用过滤掩膜
        slope_map_clean = np.where(valid_mask, slope_map, np.nan)
        pvalue_map_clean = np.where(valid_mask, pvalue_map, np.nan)
        r_squared_map_clean = np.where(valid_mask, r_squared_map, np.nan)

        # 将NaN转换为NODATA_OUT用于保存
        slope_map_save = np.where(np.isfinite(slope_map_clean), slope_map_clean, NODATA_OUT)
        pvalue_map_save = np.where(np.isfinite(pvalue_map_clean), pvalue_map_clean, NODATA_OUT)
        r2_map_save = np.where(np.isfinite(r_squared_map_clean), r_squared_map_clean, NODATA_OUT)

        if mask is not None:
            slope_map_save[~mask] = NODATA_OUT
            pvalue_map_save[~mask] = NODATA_OUT
            r2_map_save[~mask] = NODATA_OUT

        # 保存结果
        write_geotiff(output_dir / f"{resp_var}_vs_deltaSOS_slope.tif", slope_map_save, profile)
        write_geotiff(output_dir / f"{resp_var}_vs_deltaSOS_pvalue.tif", pvalue_map_save, profile)
        write_geotiff(output_dir / f"{resp_var}_vs_deltaSOS_R2.tif", r2_map_save, profile)

        # 全面统计输出（全部有效像元 + 显著性像元）
        unit = "mm/day per day ΔSOS" if resp_var == 'Trate' else "mm per day ΔSOS"
        print_comprehensive_statistics(
            slope_map_clean, pvalue_map_clean, valid_mask,
            var_name=f"{resp_var} ~ ΔSOS 回归斜率",
            unit=unit,
            print_percentiles=True,
            print_sign_split=True
        )

    print("\n  ✓ Section 3.2A 分析完成")
    print(f"  输出目录: {output_dir}")


# ==================== Section 3.3: 驱动因子分析 ====================
def section_3_3_driver_analysis(mask):
    """
    Section 3.3: Drivers of TRshape decrease with spring phenology change

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
    print("Section 3.3: TRshape 驱动因子分析")
    print("="*70)

    years = list(range(YEAR_START, YEAR_END + 1))
    n_years = len(years)

    # 读取模板（使用TRc以确保与04a/04c空间一致）
    data_first, profile, _ = read_geotiff(TRC_DIR / f"TRc_{years[0]}.tif")
    height, width = data_first.shape

    # 定义变量（Wang 2025 Eq. 3）
    # 统一因变量：移除TRtiming（时序效应），仅保留TRc和TRshape（形状效应）
    response_vars = ['TRc', 'TRshape']
    # 注意：SIFspr/SIFsum 实际使用日尺度 GPP 的季节均值作为代理
    predictor_vars = ['SOS', 'SMroot', 'Ta', 'Rs', 'P', 'SIFspr', 'SIFsum']
    available_vars = []
    missing_vars = []
    for var in predictor_vars:
        if var == 'SOS':
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
    mode_desc = "批量向量化模式" if USE_BATCH_VECTORIZED else "VIF过滤模式"
    print(f"\n[Part 1] 全时段偏相关归因分析（1982-2018，{mode_desc}）...")

    # 预计算所有预测变量
    print("\n  预计算预测变量...")
    X_all_years = {var: [] for var in predictor_vars}

    for year in tqdm(years, desc="预计算自变量"):
        # 读取SOS和POS
        sos_data, _, sos_nodata = read_geotiff(PHENO_DIR / f"sos_gpp_{year}.tif")
        pos_data, _, pos_nodata = read_geotiff(PHENO_DIR / f"pos_doy_gpp_{year}.tif")

        sos_map = np.where(_is_valid_value(sos_data, sos_nodata), sos_data, np.nan)
        pos_map = np.where(_is_valid_value(pos_data, pos_nodata), pos_data, np.nan)

        if 'SOS' in predictor_vars:
            X_all_years['SOS'].append(sos_map)

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

        # 选择worker函数（批量向量化 vs VIF过滤）
        worker_func = _partial_corr_block_worker_batch if USE_BATCH_VECTORIZED else _partial_corr_block_worker

        if USE_BLOCK_PARALLEL:
            with ProcessPoolExecutor(max_workers=MAX_WORKERS_3_3) as executor:
                futures = []
                for r0, r1, c0, c1, block_mask in blocks:
                    X_block = {var: X_all_years[var][:, r0:r1, c0:c1] for var in predictor_vars}
                    Y_block = Y_stack[:, r0:r1, c0:c1]
                    args = (r0, r1, c0, c1, block_mask, predictor_vars, X_block, Y_block, min_rows)
                    futures.append(executor.submit(worker_func, args))

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
                r0, r1, c0, c1, pr_block, pp_block, vif_block, r2_block = worker_func(args)
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

        # ========== 详细统计汇总 ==========
        print(f"\n    === {response_var} 偏相关详细统计 ===")

        # 1. 全局统计（每个预测变量）
        stats_table = []
        for var in predictor_vars:
            r_map = partial_r_maps[var]
            p_map = partial_p_maps[var]
            valid = (r_map != NODATA_OUT) & np.isfinite(r_map)

            if np.sum(valid) > 0:
                r_valid = r_map[valid]
                p_valid = p_map[valid]

                mean_r = np.mean(r_valid)
                std_r = np.std(r_valid)
                median_r = np.median(r_valid)
                sig_pct = np.sum(p_valid < 0.05) / len(p_valid) * 100
                main_driver_pct = np.sum(np.abs(r_valid) > 0.1) / len(r_valid) * 100
                n_valid = len(r_valid)

                stats_table.append({
                    'var': var,
                    'mean_R': mean_r,
                    'std_R': std_r,
                    'median_R': median_r,
                    'sig_pct': sig_pct,
                    'main_pct': main_driver_pct,
                    'n_valid': n_valid
                })

        # 按平均|R|降序排列
        stats_table.sort(key=lambda x: abs(x['mean_R']), reverse=True)

        print("\n    主导因子排名（按平均|R|降序）:")
        print(f"      {'排名':<4} {'变量':<10} {'平均R':<12} {'中位R':<10} {'显著%':<8} {'主导%':<8} {'有效像元':<10}")
        print(f"      {'-'*4} {'-'*10} {'-'*12} {'-'*10} {'-'*8} {'-'*8} {'-'*10}")
        for i, s in enumerate(stats_table, 1):
            print(f"      {i:<4} {s['var']:<10} "
                  f"{s['mean_R']:+.3f}±{s['std_R']:.3f}  "
                  f"{s['median_R']:+.3f}    "
                  f"{s['sig_pct']:>6.1f}%  "
                  f"{s['main_pct']:>6.1f}%  "
                  f"{s['n_valid']:>10}")

        # 2. R² 拟合统计 - 使用全面统计函数
        r2_mask = (r_squared_map != NODATA_OUT) & np.isfinite(r_squared_map)
        if mask is not None:
            r2_mask &= mask

        if np.sum(r2_mask) > 0:
            print(f"\n    === 模型拟合统计（R²）===")

            # 全面统计（全部像元）
            r2_stats = print_comprehensive_statistics(
                r_squared_map, pvalue_map=None, mask=r2_mask,
                var_name="R²",
                unit="",
                print_percentiles=True,
                print_sign_split=False  # R²不需要正负值统计
            )

            # 额外统计：R²阈值分布
            r2_valid = r_squared_map[r2_mask]
            print(f"\n    R²阈值分布:")
            print(f"      R² > 0.1: {np.sum(r2_valid > 0.1) / len(r2_valid) * 100:.1f}% ({np.sum(r2_valid > 0.1)} 像元)")
            print(f"      R² > 0.3: {np.sum(r2_valid > 0.3) / len(r2_valid) * 100:.1f}% ({np.sum(r2_valid > 0.3)} 像元)")
            print(f"      R² > 0.5: {np.sum(r2_valid > 0.5) / len(r2_valid) * 100:.1f}% ({np.sum(r2_valid > 0.5)} 像元)")
            print(f"      有效预测像元总数: {len(r2_valid)}")

    # ========== Part 2: 15年滑动窗口偏相关演变 ==========
    print("\n[Part 2] 15年滑动窗口偏相关演变...")

    window_size = 15
    n_windows = n_years - window_size + 1

    if n_windows < 2:
        print("  ⚠ 警告：窗口数不足，跳过滑动窗口分析")
        return

    # 对TRshape进行滑动窗口分析（论文重点关注，修改为TRshape）
    response_var = 'TRshape'
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

    # 创建进程池（复用以避免重复创建开销）
    p = len(predictor_vars)
    min_rows = max(10, p + 5)
    executor = ProcessPoolExecutor(max_workers=MAX_WORKERS_3_3) if USE_BLOCK_PARALLEL else None

    # 选择worker函数（批量向量化 vs 原版）
    window_worker_func = _partial_corr_window_block_worker_batch if USE_BATCH_VECTORIZED else _partial_corr_window_block_worker

    try:
        for win_idx in tqdm(range(n_windows), desc="滑动窗口"):
            start_idx = win_idx
            end_idx = win_idx + window_size
            start_year = years[start_idx]
            end_year = years[end_idx - 1]

            # 提取窗口数据
            Y_window = Y_stack[start_idx:end_idx]  # (window_size, H, W)
            X_window = {var: X_all_years[var][start_idx:end_idx] for var in predictor_vars}

            # 逐像元回归（分块/并行）
            if USE_BLOCK_PARALLEL:
                futures = []
                for r0, r1, c0, c1, block_mask in blocks:
                    X_block = {var: X_window[var][:, r0:r1, c0:c1] for var in predictor_vars}
                    Y_block = Y_window[:, r0:r1, c0:c1]
                    args = (r0, r1, c0, c1, block_mask, predictor_vars, X_block, Y_block, min_rows)
                    futures.append(executor.submit(window_worker_func, args))

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
                    r0, r1, c0, c1, pr_block = window_worker_func(args)
                    for var in predictor_vars:
                        partial_r_evolution[var][win_idx, r0:r1, c0:c1] = pr_block[var]
    finally:
        if executor is not None:
            executor.shutdown(wait=True)

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


def main():
    print("\n" + "="*70)
    print("统计分析模块 - Timing/Shape Decomposition (Sections 3.2 & 3.3)")
    print("="*70)

    # 性能优化状态
    if NUMBA_AVAILABLE:
        print("\n⚡ Numba JIT加速已启用（VIF + 偏相关计算约10-50x加速）")
    else:
        print("\n⚠️  Numba未安装，使用NumPy/SciPy版本（较慢）")
        print("   建议安装以获得10-50x性能提升: pip install numba")

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

    # Section 3.2: Timing/Shape组分 vs dSOS 回归分析
    try:
        section_3_2_timing_shape(mask)
    except Exception as e:
        print(f"\n  ✗ Section 3.2 分析失败: {str(e)}")
        import traceback
        traceback.print_exc()

    # Section 3.3: TRshape驱动因子分析
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
    print("  │   ├── TRshape_vs_deltaSOS_slope.tif")
    print("  │   ├── TRtiming_vs_deltaSOS_slope.tif")
    print("  │   ├── TRsos_vs_deltaSOS_slope.tif")
    print("  │   ├── TRpos_vs_deltaSOS_slope.tif")
    print("  │   └── Trate_vs_deltaSOS_slope.tif")
    print("  ├── Section_3.3_Drivers/")
    print("  │   ├── Full_Period/")
    print("  │   │   ├── TRc/")
    print("  │   │   ├── TRshape/")
    print("  │   │   └── TRtiming/")
    print("  │   │       ├── partial_r_{var}.tif (SOS, SMroot, Ta, Rs, P, SIFspr, SIFsum; SIF=GPP proxy)")
    print("  │   │       ├── partial_p_{var}.tif")
    print("  │   │       ├── vif_retained_{var}.tif")
    print("  │   │       └── R_squared.tif")
    print("  │   ├── Moving_Window/")
    print("  │   │   └── TRshape/")
    print("  │   │       └── partial_r_{var}_{start}-{end}.tif")
    print("  │   └── Sensitivity_Trends/")
    print("  │       └── TRshape/")
    print("  │           ├── {var}_trend_slope.tif")
    print("  │           └── {var}_trend_pvalue.tif")
    print("="*70)


if __name__ == "__main__":
    main()
