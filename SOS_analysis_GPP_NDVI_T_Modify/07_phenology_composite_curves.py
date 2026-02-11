#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module 07: 物候异常分组的年内生长曲线对比图

方法：
1. 像元级：根据SOS/EOS相对于多年平均的提前/推迟分类年份
2. 四类组合：SOS提前+EOS推迟, SOS提前+EOS提前, SOS推迟+EOS推迟, SOS推迟+EOS提前
3. 计算每类年份的多年平均日数据（像元级），再空间平均
4. 每个类别单独绘制一张图，叠加基准线，标注SOS/POS/EOS位置

支持变量替换：通过_config.py中的MIDDLE_VAR_NAME配置（如NDVI、GPP等）

输出：
- 每个类别一张图：该类别曲线 + 全年份基准线 + 物候标注
"""

import numpy as np
import rasterio
from pathlib import Path
from datetime import datetime, timedelta
from tqdm import tqdm
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings('ignore')

# 全局字体配置（与06_plotting.py一致）
mpl.rcParams.update({
    "font.family": "Arial",
    "axes.unicode_minus": False,
    "mathtext.fontset": "stixsans",
    "mathtext.default": "regular",
    "figure.dpi": 300,
    "savefig.bbox": "tight",
})

# 导入配置
from _config import (
    PHENO_SOURCE_DIR, DAILY_DATA_DIR, DAILY_DATA_FORMAT, OUTPUT_ROOT,
    YEAR_START, YEAR_END, NODATA_OUT, MASK_FILE, TEMPLATE_RASTER,
    PHENO_SOURCE_FORMAT, MIDDLE_VAR_NAME, FIGURES_DIR,
    RUN_CONFIGS_07
)

# ==================== 配置 ====================
# 输出目录（与06_plotting.py相同）
OUTPUT_FIG_DIR = FIGURES_DIR / "Phenology_Composite"

def ensure_output_fig_dir():
    OUTPUT_FIG_DIR.mkdir(parents=True, exist_ok=True)

# 平滑参数
SG_WINDOW = 15  # Savitzky-Golay滤波窗口
SG_POLYORDER = 3

# 物候分类阈值列表（会依次运行每个阈值）
# 0 表示无阈值（原有逻辑），0.02 表示 ±2%
PHENO_THRESHOLD_LIST = [0, 0.02, 0.05]

# 颜色配置
COLORS = {
    'baseline': '#808080',            # 灰色 - 全部年份基准
    'sos_early_eos_late': '#2E86AB',  # 蓝色 - SOS提前+EOS推迟（生长季延长）
    'sos_early_eos_early': '#A23B72', # 紫色 - SOS提前+EOS提前
    'sos_late_eos_late': '#F18F01',   # 橙色 - SOS推迟+EOS推迟
    'sos_late_eos_early': '#C73E1D',  # 红色 - SOS推迟+EOS提前（生长季缩短）
}

# 标签配置
LABELS = {
    'baseline': 'All years (baseline)',
    'sos_early_eos_late': 'SOS early + EOS late',
    'sos_early_eos_early': 'SOS early + EOS early',
    'sos_late_eos_late': 'SOS late + EOS late',
    'sos_late_eos_early': 'SOS late + EOS early',
}

# 输出文件名格式
OUTPUT_FILE_FORMAT = f"phenology_composite_{MIDDLE_VAR_NAME}_{{category}}.png"


# ==================== 辅助函数 ====================

def is_leap_year(year):
    """判断是否为闰年"""
    return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)


def noleap_doy_from_date(date_obj):
    """
    将真实日期映射为无闰日DOY（1-365）。
    闰年跳过2月29日，3月1日及之后的DOY减1。
    返回None表示2月29日（需跳过）。
    """
    doy = date_obj.timetuple().tm_yday
    if is_leap_year(date_obj.year):
        if doy == 60:
            return None
        if doy > 60:
            return doy - 1
    return doy


def load_mask():
    """加载研究区掩膜"""
    if not MASK_FILE.exists():
        raise FileNotFoundError(f"掩膜文件不存在: {MASK_FILE}")
    with rasterio.open(MASK_FILE) as src:
        mask_data = src.read(1)
        mask_nodata = src.nodata
    mask = np.isfinite(mask_data)
    if mask_nodata is not None and np.isfinite(mask_nodata):
        mask &= mask_data != mask_nodata
    mask &= mask_data > 0
    return mask


def load_phenology_all_years(years, pheno_dir=None, pheno_format=None):
    """
    加载所有年份的物候数据

    Parameters:
    -----------
    years : list
        年份列表
    pheno_dir : Path, optional
        物候数据目录，默认使用全局PHENO_SOURCE_DIR
    pheno_format : dict, optional
        物候文件命名格式，默认使用全局PHENO_SOURCE_FORMAT

    Returns:
    --------
    sos_stack : ndarray, shape (n_years, H, W)
    eos_stack : ndarray, shape (n_years, H, W)
    valid_years : list
    """
    if pheno_dir is None:
        pheno_dir = PHENO_SOURCE_DIR
    if pheno_format is None:
        pheno_format = PHENO_SOURCE_FORMAT

    sos_list = []
    eos_list = []
    valid_years = []

    for year in tqdm(years, desc="加载物候数据"):
        sos_file = pheno_dir / pheno_format['SOS'].format(year=year)
        eos_file = pheno_dir / pheno_format['EOS'].format(year=year)

        if not sos_file.exists() or not eos_file.exists():
            continue

        with rasterio.open(sos_file) as src:
            sos = src.read(1).astype(np.float32)
            sos_nodata = src.nodata

        with rasterio.open(eos_file) as src:
            eos = src.read(1).astype(np.float32)
            eos_nodata = src.nodata

        # 处理nodata
        if sos_nodata is not None:
            sos[sos == sos_nodata] = np.nan
        if eos_nodata is not None:
            eos[eos == eos_nodata] = np.nan

        # 基本有效性检查
        sos[(sos <= 0) | (sos > 365)] = np.nan
        eos[(eos <= 0) | (eos > 365)] = np.nan

        sos_list.append(sos)
        eos_list.append(eos)
        valid_years.append(year)

    sos_stack = np.stack(sos_list, axis=0)
    eos_stack = np.stack(eos_list, axis=0)

    return sos_stack, eos_stack, valid_years


def classify_years_by_phenology(sos_stack, eos_stack, valid_years, mask, threshold_pct=0.10):
    """
    按物候异常分类年份（像元级）

    Parameters:
    -----------
    sos_stack : ndarray
        SOS数据堆栈 (n_years, H, W)
    eos_stack : ndarray
        EOS数据堆栈 (n_years, H, W)
    valid_years : list
        有效年份列表
    mask : ndarray
        有效像元掩膜 (H, W)
    threshold_pct : float
        物候偏差阈值（相对于多年平均的百分比）
        例如 0.10 表示筛选偏离多年平均至少10%的年份
        设为 0 则不设阈值（原有逻辑）
    """
    n_years, H, W = sos_stack.shape

    # 计算像元级多年平均
    with np.errstate(invalid='ignore'):
        sos_mean = np.nanmean(sos_stack, axis=0)
        eos_mean = np.nanmean(eos_stack, axis=0)

    # 扩展平均值到年份维度
    sos_mean_3d = sos_mean[np.newaxis, :, :]
    eos_mean_3d = eos_mean[np.newaxis, :, :]

    # 分类掩膜（带阈值）
    # SOS提前：SOS < SOSav * (1 - threshold)，即提前至少threshold%
    # SOS推迟：SOS > SOSav * (1 + threshold)，即推迟至少threshold%
    if threshold_pct > 0:
        sos_early = sos_stack < sos_mean_3d * (1 - threshold_pct)
        sos_late = sos_stack > sos_mean_3d * (1 + threshold_pct)
        eos_early = eos_stack < eos_mean_3d * (1 - threshold_pct)
        eos_late = eos_stack > eos_mean_3d * (1 + threshold_pct)
        print(f"  使用物候偏差阈值: ±{threshold_pct*100:.0f}%")
    else:
        # 原有逻辑：无阈值
        sos_early = sos_stack < sos_mean_3d
        sos_late = sos_stack >= sos_mean_3d
        eos_early = eos_stack < eos_mean_3d
        eos_late = eos_stack >= eos_mean_3d
        print("  无物候偏差阈值（原有逻辑）")

    # 四类组合
    categories = {
        'sos_early_eos_late': sos_early & eos_late,
        'sos_early_eos_early': sos_early & eos_early,
        'sos_late_eos_late': sos_late & eos_late,
        'sos_late_eos_early': sos_late & eos_early,
    }

    # 统计每个像元每个类别的年份数
    category_year_counts = {}
    for cat, cat_mask in categories.items():
        valid_pheno = np.isfinite(sos_stack) & np.isfinite(eos_stack)
        cat_valid = cat_mask & valid_pheno
        category_year_counts[cat] = np.sum(cat_valid, axis=0)

    return categories, category_year_counts, sos_mean, eos_mean


def _read_daily_file(args):
    """并行读取单个日数据文件（用于线程池）"""
    data_file, mask, data_nodata_default = args
    try:
        with rasterio.open(data_file) as src:
            data = src.read(1).astype(np.float32)
            data_nodata = src.nodata if src.nodata is not None else data_nodata_default
        valid = mask & np.isfinite(data)
        if data_nodata is not None:
            valid &= (data != data_nodata)
        return data, valid
    except Exception:
        return None, None


def compute_daily_data_by_category(years, categories, mask, use_cache=True, n_workers=8,
                                    daily_dir=None, daily_format=None, var_label=None,
                                    threshold_pct=0):
    """
    计算每个类别的像元级多年平均日数据，然后空间平均

    Parameters:
    -----------
    daily_dir : Path, optional
        日尺度数据目录，默认使用全局DAILY_DATA_DIR
    daily_format : str, optional
        日尺度数据文件命名格式，默认使用全局DAILY_DATA_FORMAT
    var_label : str, optional
        变量标签（用于缓存文件名和打印），默认使用全局MIDDLE_VAR_NAME
    threshold_pct : float
        物候阈值。>0时基准曲线只用四类并集的像元年份
    """
    from concurrent.futures import ThreadPoolExecutor
    import hashlib

    if daily_dir is None:
        daily_dir = DAILY_DATA_DIR
    if daily_format is None:
        daily_format = DAILY_DATA_FORMAT
    if var_label is None:
        var_label = MIDDLE_VAR_NAME

    H, W = mask.shape
    n_years = len(years)

    # 生成缓存文件名（基于变量标签、年份范围、类别和阈值）
    cache_key = f"{var_label}_{years[0]}_{years[-1]}_{len(categories)}_thr{threshold_pct:.2f}"
    cache_hash = hashlib.md5(cache_key.encode()).hexdigest()[:8]
    cache_file = OUTPUT_FIG_DIR / f"daily_data_cache_{cache_hash}.npz"

    # 尝试加载缓存
    if use_cache and cache_file.exists():
        print(f"\n从缓存加载日数据: {cache_file}")
        try:
            cached = np.load(cache_file, allow_pickle=True)
            daily_data = dict(cached['daily_data'].item())
            print(f"  缓存加载成功，包含 {len(daily_data)} 个类别")
            return daily_data
        except Exception as e:
            print(f"  缓存加载失败: {e}，重新计算...")

    print(f"\n计算各类别的日{var_label}...")

    # 阈值>0时，基准曲线只用四类并集的像元年份
    if threshold_pct > 0:
        any_category_mask = None
        for cat_mask in categories.values():
            if any_category_mask is None:
                any_category_mask = cat_mask.copy()
            else:
                any_category_mask |= cat_mask
        print(f"  基准曲线使用四类并集的像元年份（阈值={threshold_pct*100:.0f}%）")
    else:
        any_category_mask = None  # None表示使用所有像元年份
        print(f"  基准曲线使用所有像元年份（无阈值）")

    # 初始化累加器（使用更小的数据类型减少内存）
    baseline_sum = np.zeros((365, H, W), dtype=np.float32)
    baseline_count = np.zeros((365, H, W), dtype=np.uint8)

    cat_sum = {cat: np.zeros((365, H, W), dtype=np.float32) for cat in categories}
    cat_count = {cat: np.zeros((365, H, W), dtype=np.uint8) for cat in categories}

    # 预先构建所有文件路径列表
    file_tasks = []
    for year_idx, year in enumerate(years):
        days_in_year = 366 if is_leap_year(year) else 365
        for day_offset in range(days_in_year):
            date_obj = datetime(year, 1, 1) + timedelta(days=day_offset)
            doy = noleap_doy_from_date(date_obj)
            if doy is None:
                continue
            data_file = daily_dir / daily_format.format(date=date_obj.strftime("%Y%m%d"))
            if data_file.exists():
                file_tasks.append((year_idx, doy - 1, data_file))

    print(f"  共 {len(file_tasks)} 个文件待处理")

    # 使用多线程并行读取和处理（分批处理以控制内存）
    batch_size = 100  # 每批处理100个文件
    n_batches = (len(file_tasks) + batch_size - 1) // batch_size

    for batch_idx in tqdm(range(n_batches), desc=f"处理{var_label}数据"):
        batch_start = batch_idx * batch_size
        batch_end = min(batch_start + batch_size, len(file_tasks))
        batch_tasks = file_tasks[batch_start:batch_end]

        # 并行读取这批文件
        read_args = [(task[2], mask, NODATA_OUT) for task in batch_tasks]

        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            results = list(executor.map(_read_daily_file, read_args))

        # 处理结果
        for (year_idx, doy_idx, _), (data, valid) in zip(batch_tasks, results):
            if data is None:
                continue

            # 基准曲线累加：阈值>0时只用四类并集的像元年份
            if any_category_mask is not None:
                baseline_valid = valid & any_category_mask[year_idx]
            else:
                baseline_valid = valid
            baseline_sum[doy_idx, baseline_valid] += data[baseline_valid]
            baseline_count[doy_idx, baseline_valid] += 1

            for cat, cat_mask in categories.items():
                cat_year_mask = cat_mask[year_idx]
                cat_valid = valid & cat_year_mask
                if np.any(cat_valid):
                    cat_sum[cat][doy_idx, cat_valid] += data[cat_valid]
                    cat_count[cat][doy_idx, cat_valid] += 1

    # 计算空间平均（向量化）
    print("计算空间平均...")

    daily_data = {}
    mask_3d = mask[np.newaxis, :, :]  # (1, H, W)

    # 基准线 - 时间平均（每个像元在多年上的平均值）
    with np.errstate(divide='ignore', invalid='ignore'):
        baseline_mean = np.where(baseline_count > 0,
                                  baseline_sum / baseline_count,
                                  np.nan)

    # 各类别 - 时间平均
    cat_mean_dict = {}
    for cat in categories:
        with np.errstate(divide='ignore', invalid='ignore'):
            cat_mean_dict[cat] = np.where(cat_count[cat] > 0,
                                           cat_sum[cat] / cat_count[cat],
                                           np.nan)

    # 构建公共有效像元掩膜：对每个DOY，只有基准线和所有4个类别都有有效数据的像元才参与空间平均
    # 这确保基准曲线确实是四类曲线的加权平均，避免因不同像元集合导致的系统性偏差
    common_valid = mask_3d & np.isfinite(baseline_mean)  # (365, H, W)
    for cat in categories:
        common_valid = common_valid & np.isfinite(cat_mean_dict[cat])

    n_common = np.sum(common_valid, axis=(1, 2))
    print(f"  公共有效像元数: 最小={np.min(n_common)}, 最大={np.max(n_common)}, 平均={np.mean(n_common):.0f}")

    # 使用公共掩膜进行空间平均
    baseline_mean_masked = np.where(common_valid, baseline_mean, np.nan)
    baseline_spatial = np.nanmean(baseline_mean_masked, axis=(1, 2))  # (365,)
    daily_data['baseline'] = baseline_spatial

    for cat in categories:
        cat_mean_masked = np.where(common_valid, cat_mean_dict[cat], np.nan)
        cat_spatial = np.nanmean(cat_mean_masked, axis=(1, 2))
        daily_data[cat] = cat_spatial

    # 保存缓存
    if use_cache:
        print(f"保存缓存到: {cache_file}")
        np.savez_compressed(cache_file, daily_data=daily_data)

    # 释放大数组内存
    del baseline_sum, baseline_count, cat_sum, cat_count, cat_mean_dict, common_valid

    return daily_data


def smooth_curve(y, window=SG_WINDOW, polyorder=SG_POLYORDER):
    """Savitzky-Golay平滑"""
    valid = np.isfinite(y)
    if np.sum(valid) < window:
        return y

    y_filled = y.copy()
    if np.any(~valid):
        x = np.arange(len(y))
        y_filled = np.interp(x, x[valid], y[valid])

    return savgol_filter(y_filled, window, polyorder, mode='nearest')


# ==================== 双Logistic拟合 ====================

def stable_sigmoid(x, clip=30.0):
    x = np.clip(x, -clip, clip)
    return 1.0 / (1.0 + np.exp(-x))


def double_logistic(t, c1, c2, c3, r1, r2, t1, t2):
    rise = stable_sigmoid(r1 * (t - t1))
    fall = stable_sigmoid(r2 * (t - t2))
    return c1 + c2 * rise - c3 * fall


def fit_double_logistic(x, y, data_type=None):
    """拟合双Logistic曲线，根据数据类型自动调整边界"""
    if data_type is None:
        data_type = MIDDLE_VAR_NAME.lower()

    valid = np.isfinite(y)
    if np.sum(valid) < 20:
        return None

    x_valid = x[valid]
    y_valid = y[valid]

    data_min, data_max = np.min(y_valid), np.max(y_valid)
    dr = data_max - data_min

    print(f"    [DEBUG拟合] 数据范围: min={data_min:.6f}, max={data_max:.6f}, 振幅dr={dr:.6f}")

    # 对于GPP，振幅阈值可以更严格一些；对于NDVI保持原阈值
    min_amplitude = 0.01 if data_type == 'ndvi' else 0.01
    if dr < min_amplitude:
        print(f"    [DEBUG拟合] 振幅dr={dr:.6f} < {min_amplitude}，拒绝拟合")
        return None

    p0 = [data_min, dr/2, dr/2, 0.1, 0.1, 100, 250]

    if data_type == 'ndvi':
        # NDVI固定边界（值域0-1）
        bounds_lower = [-0.5, 0, 0, 0, 0,   0,  50]
        bounds_upper = [ 1.5, 2, 2, 2, 2, 183, 366]
    else:  # gpp或其他变量
        # 动态边界，根据实际数据范围自适应
        bounds_lower = [data_min-dr, 0,     0,     0, 0,   0,  50]
        bounds_upper = [data_max+dr, dr*3,  dr*3,  2, 2, 183, 366]

    try:
        popt, _ = curve_fit(double_logistic, x_valid, y_valid,
                            p0=p0, bounds=(bounds_lower, bounds_upper),
                            maxfev=10000)
        return popt
    except Exception as e:
        print(f"    [DEBUG拟合] curve_fit异常: {type(e).__name__}")
        return None


def extract_phenology_from_curve(curve, x_daily, percent=0.20):
    """从曲线提取物候参数"""
    y = np.asarray(curve, float)
    smin, smax = np.nanmin(y), np.nanmax(y)

    print(f"      [DEBUG] 曲线范围: min={smin:.6f}, max={smax:.6f}, 振幅={smax-smin:.6f}")

    if not np.isfinite(smin) or smax - smin < 1e-6:
        print(f"      [DEBUG] 振幅检查失败")
        return np.nan, np.nan, np.nan

    yn = (y - smin) / (smax - smin)
    peak_idx = int(np.nanargmax(yn))

    print(f"      [DEBUG] peak_idx={peak_idx}, peak_doy={x_daily[peak_idx]}")

    left = yn[:peak_idx+1]
    right = yn[peak_idx:]

    left_min_idx = int(np.nanargmin(left))
    right_min_idx = peak_idx + int(np.nanargmin(right))
    minL, minR = yn[left_min_idx], yn[right_min_idx]
    peak = yn[peak_idx]

    thrL = (peak - minL) * percent + minL
    thrR = (peak - minR) * percent + minR

    sos_doy = np.nan
    for i in range(max(left_min_idx, 0), peak_idx+1):
        if yn[i] >= thrL:
            sos_doy = x_daily[i]
            break

    eos_doy = np.nan
    for i in range(peak_idx, right_min_idx+1):
        if yn[i] >= thrR:
            eos_doy = x_daily[i]

    pos_doy = x_daily[peak_idx]

    return sos_doy, eos_doy, pos_doy


# ==================== 绘图函数 ====================

def plot_decomposition_schematic(baseline_curve, cat_curve, baseline_pheno, cat_pheno,
                                  category, output_path, var_label=None):
    """
    绘制完整的分解示意图，在一张图上展示整个生长季的分解

    分解公式（参照03c）：
      TRc_y = TRc_av + TR_sos_change + TR_pos_change + TR_eos_change + TR_fixed_window

    填充区域：
    - TR_c_av: 基准曲线在[SOSav, EOSav]窗口内的累积（灰色填充）
    - TR_sos_change: SOS变化带来的累积（绿色）
    - TR_pos_change: POS变化带来的累积（金黄色，注意两侧方向相反会抵消）
    - TR_eos_change: EOS变化带来的累积（青色）
    - TR_fixed_window: 固定窗口内的速率差异（紫色斜线填充）
    """
    if var_label is None:
        var_label = MIDDLE_VAR_NAME

    fig, ax = plt.subplots(figsize=(18, 12))

    x_daily = np.arange(1, 366)

    # 获取物候参数
    sos_av = baseline_pheno.get('SOS', 100)
    pos_av = baseline_pheno.get('POS', 200)
    eos_av = baseline_pheno.get('EOS', 280)

    sos_y = cat_pheno.get('SOS', 90)
    pos_y = cat_pheno.get('POS', 200)
    eos_y = cat_pheno.get('EOS', 290)

    # 确保物候参数有效
    if not all(np.isfinite([sos_av, pos_av, eos_av, sos_y, pos_y, eos_y])):
        print(f"  物候参数无效，跳过分解示意图: {category}")
        return

    # 转换为整数索引
    sos_av_idx = max(0, min(364, int(sos_av) - 1))
    pos_av_idx = max(0, min(364, int(pos_av) - 1))
    eos_av_idx = max(0, min(364, int(eos_av) - 1))
    sos_y_idx = max(0, min(364, int(sos_y) - 1))
    pos_y_idx = max(0, min(364, int(pos_y) - 1))
    eos_y_idx = max(0, min(364, int(eos_y) - 1))

    # 颜色定义
    color_baseline = '#3498db'  # 蓝色 - 基准曲线
    color_year = '#e67e22'      # 橙色 - 某年曲线
    color_trc_av = '#E0E0E0'    # 浅灰色 - TR_c_av 基准累积
    color_sos = '#2ecc71'       # 绿色 - TR_sos_change
    color_pos = '#f39c12'       # 金黄色 - TR_pos_change
    color_eos = '#1abc9c'       # 青色 - TR_eos_change
    # TR_fixed: 使用对比明显的颜色区分正负差异
    color_fixed_pos = '#e74c3c'   # 红色 - TR_fixed正差异（当年>基准）
    color_fixed_neg = '#2980b9'   # 深蓝色 - TR_fixed负差异（当年<基准）

    # 字体大小配置：标题独立，其余统一
    FONTSIZE_TITLE = 38
    FONTSIZE_BODY = 36  # 所有非标题字体统一大小

    # 预先计算y轴下限（用于fill_between的底部）
    all_values = np.concatenate([baseline_curve, cat_curve])
    valid_values = all_values[np.isfinite(all_values)]
    data_min = np.min(valid_values)
    if var_label.upper() == 'NDVI':
        y_bottom = max(0, np.floor(data_min / 0.05) * 0.05)
    elif 'GPP' in var_label.upper() or 'TR' in var_label.upper():
        y_bottom = 0
    else:
        y_bottom = max(0, np.floor(data_min * 0.9 * 10) / 10)

    # ==================== 1. TR_c_av: 基准累积 [SOSav, EOSav] ====================
    x_range_av = x_daily[sos_av_idx:eos_av_idx+1]
    y_baseline_av = baseline_curve[sos_av_idx:eos_av_idx+1]
    ax.fill_between(x_range_av, y_bottom, y_baseline_av,
                    color=color_trc_av, alpha=0.6, label=r'$TR_{c,av}$ (baseline)')

    # ==================== 2. TR_sos_change: SOS变化分量 ====================
    if sos_y < sos_av:  # SOS提前 → 正贡献
        x_sos = x_daily[sos_y_idx:sos_av_idx+1]
        y_sos = cat_curve[sos_y_idx:sos_av_idx+1]
        ax.fill_between(x_sos, y_bottom, y_sos,
                        color=color_sos, alpha=0.5, label=r'$\Delta TR_{SOS}$ (+)')
    elif sos_y > sos_av:  # SOS推迟 → 负贡献
        x_sos = x_daily[sos_av_idx:sos_y_idx+1]
        y_sos = cat_curve[sos_av_idx:sos_y_idx+1]
        ax.fill_between(x_sos, y_bottom, y_sos,
                        color=color_sos, alpha=0.5,
                        label=r'$\Delta TR_{SOS}$ (−)')

    # ==================== 3. TR_pos_change: POS变化分量（两侧都绘制） ====================
    # 注意：POS变化在左右两侧的贡献方向相反
    # 左侧 [SOSav, POSav]: POS延后=正贡献, POS提前=负贡献
    # 右侧 [POSav, EOSav]: POS延后=负贡献, POS提前=正贡献
    # 两者会相互抵消，但仍分别绘制以展示变化区域

    if pos_y > pos_av:  # POS延后
        # 左右两侧贡献相互抵消，用同一色块表示
        x_pos = x_daily[pos_av_idx:pos_y_idx+1]
        y_pos = cat_curve[pos_av_idx:pos_y_idx+1]
        ax.fill_between(x_pos, y_bottom, y_pos,
                        color=color_pos, alpha=0.5, label=r'$\Delta TR_{POS}$ (cancel)')
    elif pos_y < pos_av:  # POS提前
        # 左右两侧贡献相互抵消，用同一色块表示
        x_pos = x_daily[pos_y_idx:pos_av_idx+1]
        y_pos = cat_curve[pos_y_idx:pos_av_idx+1]
        ax.fill_between(x_pos, y_bottom, y_pos,
                        color=color_pos, alpha=0.5, label=r'$\Delta TR_{POS}$ (cancel)')

    # ==================== 4. TR_eos_change: EOS变化分量 ====================
    if eos_y > eos_av:  # EOS延后 → 正贡献
        x_eos = x_daily[eos_av_idx:eos_y_idx+1]
        y_eos = cat_curve[eos_av_idx:eos_y_idx+1]
        ax.fill_between(x_eos, y_bottom, y_eos,
                        color=color_eos, alpha=0.5, label=r'$\Delta TR_{EOS}$ (+)')
    elif eos_y < eos_av:  # EOS提前 → 负贡献
        x_eos = x_daily[eos_y_idx:eos_av_idx+1]
        y_eos = cat_curve[eos_y_idx:eos_av_idx+1]
        ax.fill_between(x_eos, y_bottom, y_eos,
                        color=color_eos, alpha=0.5,
                        label=r'$\Delta TR_{EOS}$ (−)')

    # ==================== 5. TR_fixed_window: 固定窗口内速率差异 ====================
    # 在整个固定窗口[SOSav, EOSav]内绘制
    x_fixed = x_daily[sos_av_idx:eos_av_idx+1]
    y_baseline_fixed = baseline_curve[sos_av_idx:eos_av_idx+1]
    y_year_fixed = cat_curve[sos_av_idx:eos_av_idx+1]

    # 正差异（当年 > 基准）- 红色填充
    ax.fill_between(x_fixed, y_baseline_fixed, y_year_fixed,
                    where=(y_year_fixed > y_baseline_fixed),
                    color=color_fixed_pos, alpha=0.5,
                    label=r'$\Delta TR_{fixed}$ (+)')
    # 负差异（当年 < 基准）- 蓝色填充
    ax.fill_between(x_fixed, y_baseline_fixed, y_year_fixed,
                    where=(y_year_fixed <= y_baseline_fixed),
                    color=color_fixed_neg, alpha=0.5,
                    label=r'$\Delta TR_{fixed}$ (−)')

    # ==================== 6. 绘制曲线 ====================
    ax.plot(x_daily, baseline_curve, '-', color=color_baseline,
            linewidth=4.0, label='Baseline (all years)', zorder=5)
    ax.plot(x_daily, cat_curve, '-', color=color_year,
            linewidth=4.0, label=f'{LABELS[category]}', zorder=6)

    # ==================== 7. 标注关键物候日期 ====================
    # 垂直虚线（统一使用 '--' 样式）
    DASH_STYLE = (8, 5)  # (线段长, 间隔长)
    ax.axvline(sos_av, color=color_baseline, linestyle='--', linewidth=3.0, alpha=0.7, dashes=DASH_STYLE)
    ax.axvline(sos_y, color=color_year, linestyle='--', linewidth=3.0, alpha=0.7, dashes=DASH_STYLE)
    ax.axvline(pos_av, color=color_baseline, linestyle='--', linewidth=3.0, alpha=0.7, dashes=DASH_STYLE)
    ax.axvline(pos_y, color=color_year, linestyle='--', linewidth=3.0, alpha=0.7, dashes=DASH_STYLE)
    ax.axvline(eos_av, color=color_baseline, linestyle='--', linewidth=3.0, alpha=0.7, dashes=DASH_STYLE)
    ax.axvline(eos_y, color=color_year, linestyle='--', linewidth=3.0, alpha=0.7, dashes=DASH_STYLE)

    # 设置x轴范围：完整一年 1-365
    ax.set_xlim(1, 365)

    # 设置y轴范围（y_bottom已在前面计算）
    ax.set_ylim(bottom=y_bottom)

    # 物候标注文字（原始位置，zorder最高以浮于所有图层之上）
    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min
    TXT_Z = 20  # 高于曲线(5-6)和填充区域
    FONTSIZE_ANNOT = 40  # 标注字体稍大，使下标/上标更清晰

    # SOS标注（底部）：SOS'在上，SOSav在下（间距0.05统一）
    ax.text(sos_y, y_min + y_range * 0.08, r"$\mathrm{SOS\!'}$", fontsize=FONTSIZE_ANNOT, color=color_year,
            ha='center', va='bottom', fontweight='bold', zorder=TXT_Z)
    ax.text(sos_av, y_min + y_range * 0.03, r"$\mathrm{SOS}_{\!\mathrm{av}}$", fontsize=FONTSIZE_ANNOT, color=color_baseline,
            ha='center', va='bottom', fontweight='bold', zorder=TXT_Z)

    # POS标注（顶部）：POS'在上，POSav在下（间距0.05统一）
    ax.text(pos_y, y_min + y_range * 0.92, r"$\mathrm{POS\!'}$", fontsize=FONTSIZE_ANNOT, color=color_year,
            ha='center', va='top', fontweight='bold', zorder=TXT_Z)
    ax.text(pos_av, y_min + y_range * 0.87, r"$\mathrm{POS}_{\!\mathrm{av}}$", fontsize=FONTSIZE_ANNOT, color=color_baseline,
            ha='center', va='top', fontweight='bold', zorder=TXT_Z)

    # EOS标注（底部）：EOS'在上，EOSav在下（间距0.05统一）
    ax.text(eos_y, y_min + y_range * 0.08, r"$\mathrm{EOS\!'}$", fontsize=FONTSIZE_ANNOT, color=color_year,
            ha='center', va='bottom', fontweight='bold', zorder=TXT_Z)
    ax.text(eos_av, y_min + y_range * 0.03, r"$\mathrm{EOS}_{\!\mathrm{av}}$", fontsize=FONTSIZE_ANNOT, color=color_baseline,
            ha='center', va='bottom', fontweight='bold', zorder=TXT_Z)

    # ==================== 8. 图例和标签 ====================
    ax.set_ylabel(f'Daily {var_label}', fontsize=FONTSIZE_BODY, fontweight='bold')
    ax.tick_params(axis='both', labelsize=FONTSIZE_BODY)
    ax.grid(True, alpha=0.3, linestyle='--')

    # 图例：无边框，字体略小于正文
    FONTSIZE_LEGEND = 28
    ax.legend(loc='upper left', fontsize=FONTSIZE_LEGEND, frameon=False, ncol=1,
              handletextpad=0.4, labelspacing=0.3, handlelength=1.5)

    # 标题：灰色矩形紧贴数据框顶部，与数据区等宽（使用 ax.transAxes 坐标）
    title_text = f'SOS effects on {var_label}: {LABELS[category]}'
    ax.add_patch(mpl.patches.Rectangle(
        (0, 1.0), 1.0, 0.09,
        transform=ax.transAxes, facecolor='#EDEDED',
        edgecolor='#BBBBBB', linewidth=1.5, clip_on=False, zorder=10))
    ax.text(0.5, 1.045, title_text, fontsize=FONTSIZE_TITLE, fontweight='bold',
            ha='center', va='center', transform=ax.transAxes, zorder=11)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved decomposition schematic: {output_path}")


def plot_single_category(daily_data, category, category_stats, output_path,
                         baseline_curve, baseline_pheno, cat_curve, cat_pheno,
                         var_label=None):
    """
    绘制单个类别的生长曲线图（与基准线叠加）

    Parameters:
    -----------
    daily_data : dict
        各类别的日数据
    category : str
        类别名称
    category_stats : dict
        类别统计信息
    output_path : Path
        输出路径
    baseline_curve : ndarray
        基准线拟合曲线
    baseline_pheno : dict
        基准线物候参数
    cat_curve : ndarray
        该类别拟合曲线
    cat_pheno : dict
        该类别物候参数
    var_label : str, optional
        变量标签，默认使用全局MIDDLE_VAR_NAME
    """
    if var_label is None:
        var_label = MIDDLE_VAR_NAME
    fig, ax = plt.subplots(figsize=(14, 8))

    x_daily = np.arange(1, 366)

    # 绘制基准线（灰色，较粗）
    ax.plot(x_daily, baseline_curve, '-', color=COLORS['baseline'],
            linewidth=3.0, alpha=0.7, label=LABELS['baseline'], zorder=1)

    # 绘制该类别曲线（彩色，较粗）
    cat_color = COLORS[category]
    cat_label = LABELS[category]
    if category in category_stats:
        n_years = category_stats[category]['n_years_mean']
        cat_label += f" (n={n_years:.0f} yrs/pixel)"

    ax.plot(x_daily, cat_curve, '-', color=cat_color,
            linewidth=3.0, alpha=1.0, label=cat_label, zorder=2)

    # 标注基准线的SOS、POS、EOS
    y_min, y_max = ax.get_ylim()
    text_offset_y = (y_max - y_min) * 0.03

    if np.isfinite(baseline_pheno['SOS']):
        sos_doy = baseline_pheno['SOS']
        sos_idx = int(sos_doy) - 1
        sos_y = baseline_curve[sos_idx] if 0 <= sos_idx < 365 else y_min
        ax.axvline(sos_doy, color=COLORS['baseline'], linestyle='--', linewidth=1.5, alpha=0.6)
        ax.plot(sos_doy, sos_y, 'o', color=COLORS['baseline'], markersize=8, zorder=3)
        ax.text(sos_doy, sos_y + text_offset_y, f'SOS={sos_doy:.0f}',
                fontsize=10, color=COLORS['baseline'], ha='center', va='bottom', fontweight='bold')

    if np.isfinite(baseline_pheno['POS']):
        pos_doy = baseline_pheno['POS']
        pos_idx = int(pos_doy) - 1
        pos_y = baseline_curve[pos_idx] if 0 <= pos_idx < 365 else y_max
        ax.axvline(pos_doy, color=COLORS['baseline'], linestyle=':', linewidth=1.5, alpha=0.6)
        ax.plot(pos_doy, pos_y, 's', color=COLORS['baseline'], markersize=8, zorder=3)
        ax.text(pos_doy, pos_y + text_offset_y, f'POS={pos_doy:.0f}',
                fontsize=10, color=COLORS['baseline'], ha='center', va='bottom', fontweight='bold')

    if np.isfinite(baseline_pheno['EOS']):
        eos_doy = baseline_pheno['EOS']
        eos_idx = int(eos_doy) - 1
        eos_y = baseline_curve[eos_idx] if 0 <= eos_idx < 365 else y_min
        ax.axvline(eos_doy, color=COLORS['baseline'], linestyle='--', linewidth=1.5, alpha=0.6)
        ax.plot(eos_doy, eos_y, '^', color=COLORS['baseline'], markersize=8, zorder=3)
        ax.text(eos_doy, eos_y + text_offset_y, f'EOS={eos_doy:.0f}',
                fontsize=10, color=COLORS['baseline'], ha='center', va='bottom', fontweight='bold')

    # 标注该类别的SOS、POS、EOS
    if np.isfinite(cat_pheno['SOS']):
        sos_doy = cat_pheno['SOS']
        sos_idx = int(sos_doy) - 1
        sos_y = cat_curve[sos_idx] if 0 <= sos_idx < 365 else y_min
        ax.plot(sos_doy, sos_y, 'o', color=cat_color, markersize=10, zorder=4,
                markeredgecolor='white', markeredgewidth=1.5)
        ax.text(sos_doy, sos_y - text_offset_y * 2, f'SOS={sos_doy:.0f}',
                fontsize=10, color=cat_color, ha='center', va='top', fontweight='bold')

    if np.isfinite(cat_pheno['POS']):
        pos_doy = cat_pheno['POS']
        pos_idx = int(pos_doy) - 1
        pos_y = cat_curve[pos_idx] if 0 <= pos_idx < 365 else y_max
        ax.plot(pos_doy, pos_y, 's', color=cat_color, markersize=10, zorder=4,
                markeredgecolor='white', markeredgewidth=1.5)
        ax.text(pos_doy, pos_y - text_offset_y * 2, f'POS={pos_doy:.0f}',
                fontsize=10, color=cat_color, ha='center', va='top', fontweight='bold')

    if np.isfinite(cat_pheno['EOS']):
        eos_doy = cat_pheno['EOS']
        eos_idx = int(eos_doy) - 1
        eos_y = cat_curve[eos_idx] if 0 <= eos_idx < 365 else y_min
        ax.plot(eos_doy, eos_y, '^', color=cat_color, markersize=10, zorder=4,
                markeredgecolor='white', markeredgewidth=1.5)
        ax.text(eos_doy, eos_y - text_offset_y * 2, f'EOS={eos_doy:.0f}',
                fontsize=10, color=cat_color, ha='center', va='top', fontweight='bold')

    # 添加物候差异标注
    _add_phenology_difference_annotation(ax, baseline_pheno, cat_pheno,
                                          baseline_curve, cat_curve, cat_color)

    # 图形设置
    ax.set_xlabel('Day of Year (DOY)', fontsize=14, fontweight='bold')
    ax.set_ylabel(var_label, fontsize=14, fontweight='bold')

    title = f'{var_label} Growth Curve: {LABELS[category]}'
    ax.set_title(title, fontsize=16, fontweight='bold', pad=15)

    ax.set_xlim(1, 365)
    ax.set_ylim(bottom=0)

    ax.legend(loc='upper left', fontsize=11, framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle='--')

    # 添加月份刻度
    month_starts = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
    month_names = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                   'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(month_starts)
    ax2.set_xticklabels(month_names, fontsize=10)

    # 添加物候差异统计文本框
    _add_stats_textbox(ax, baseline_pheno, cat_pheno, category)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  Saved: {output_path}")


def _add_phenology_difference_annotation(ax, baseline_pheno, cat_pheno,
                                          baseline_curve, cat_curve, cat_color):
    """添加物候差异的箭头标注"""
    x_daily = np.arange(1, 366)

    # SOS差异箭头
    if np.isfinite(baseline_pheno['SOS']) and np.isfinite(cat_pheno['SOS']):
        base_sos = baseline_pheno['SOS']
        cat_sos = cat_pheno['SOS']
        if abs(cat_sos - base_sos) > 1:  # 差异大于1天才标注
            # 在y轴底部画水平箭头表示SOS差异
            y_arrow = ax.get_ylim()[0] + (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.08
            ax.annotate('', xy=(cat_sos, y_arrow), xytext=(base_sos, y_arrow),
                        arrowprops=dict(arrowstyle='<->', color=cat_color, lw=2))
            mid_x = (base_sos + cat_sos) / 2
            diff = cat_sos - base_sos
            diff_text = f'{diff:+.0f}d' if diff != 0 else '0d'
            ax.text(mid_x, y_arrow * 0.7, f'ΔSOS={diff_text}',
                    fontsize=9, color=cat_color, ha='center', va='top', fontweight='bold')

    # EOS差异箭头
    if np.isfinite(baseline_pheno['EOS']) and np.isfinite(cat_pheno['EOS']):
        base_eos = baseline_pheno['EOS']
        cat_eos = cat_pheno['EOS']
        if abs(cat_eos - base_eos) > 1:
            y_arrow = ax.get_ylim()[0] + (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.15
            ax.annotate('', xy=(cat_eos, y_arrow), xytext=(base_eos, y_arrow),
                        arrowprops=dict(arrowstyle='<->', color=cat_color, lw=2))
            mid_x = (base_eos + cat_eos) / 2
            diff = cat_eos - base_eos
            diff_text = f'{diff:+.0f}d' if diff != 0 else '0d'
            ax.text(mid_x, y_arrow * 0.85, f'ΔEOS={diff_text}',
                    fontsize=9, color=cat_color, ha='center', va='top', fontweight='bold')


def _add_stats_textbox(ax, baseline_pheno, cat_pheno, category):
    """添加物候差异统计文本框"""
    # 计算差异
    stats_lines = []

    if np.isfinite(baseline_pheno['SOS']) and np.isfinite(cat_pheno['SOS']):
        diff_sos = cat_pheno['SOS'] - baseline_pheno['SOS']
        stats_lines.append(f"ΔSOS: {diff_sos:+.1f} days")

    if np.isfinite(baseline_pheno['POS']) and np.isfinite(cat_pheno['POS']):
        diff_pos = cat_pheno['POS'] - baseline_pheno['POS']
        stats_lines.append(f"ΔPOS: {diff_pos:+.1f} days")

    if np.isfinite(baseline_pheno['EOS']) and np.isfinite(cat_pheno['EOS']):
        diff_eos = cat_pheno['EOS'] - baseline_pheno['EOS']
        stats_lines.append(f"ΔEOS: {diff_eos:+.1f} days")

    # 计算LOS (Length of Season)
    if (np.isfinite(baseline_pheno['SOS']) and np.isfinite(baseline_pheno['EOS']) and
        np.isfinite(cat_pheno['SOS']) and np.isfinite(cat_pheno['EOS'])):
        base_los = baseline_pheno['EOS'] - baseline_pheno['SOS']
        cat_los = cat_pheno['EOS'] - cat_pheno['SOS']
        diff_los = cat_los - base_los
        stats_lines.append(f"ΔLOS: {diff_los:+.1f} days")

    if stats_lines:
        stats_text = '\n'.join(stats_lines)
        props = dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9, edgecolor='gray')
        ax.text(0.98, 0.75, stats_text, transform=ax.transAxes, fontsize=11,
                verticalalignment='top', horizontalalignment='right',
                bbox=props, fontfamily='monospace')


def compute_category_statistics(categories, sos_stack, eos_stack, mask):
    """计算各类别的统计信息"""
    stats = {}

    for cat, cat_mask in categories.items():
        valid_pheno = np.isfinite(sos_stack) & np.isfinite(eos_stack)
        cat_valid = cat_mask & valid_pheno
        n_years_per_pixel = np.sum(cat_valid, axis=0)
        n_years_in_mask = n_years_per_pixel[mask]

        stats[cat] = {
            'n_years_mean': np.mean(n_years_in_mask) if n_years_in_mask.size > 0 else 0,
            'n_years_min': np.min(n_years_in_mask) if n_years_in_mask.size > 0 else 0,
            'n_years_max': np.max(n_years_in_mask) if n_years_in_mask.size > 0 else 0,
            'n_pixels_with_data': np.sum(n_years_per_pixel > 0),
        }

    return stats


# ==================== 主函数 ====================

def main():
    ensure_output_fig_dir()
    print("\n" + "="*70)
    print("物候异常分组的年内生长曲线分析（多数据源组合 × 多阈值）")
    print("="*70)

    years = list(range(YEAR_START, YEAR_END + 1))
    print(f"分析年份: {YEAR_START}-{YEAR_END} ({len(years)}年)")
    print(f"数据源组合: {len(RUN_CONFIGS_07)} 种")
    print(f"阈值列表: {PHENO_THRESHOLD_LIST}")

    # 加载掩膜（所有组合共用）
    print("\n[0] 加载掩膜...")
    mask = load_mask()
    print(f"  有效像元数: {np.sum(mask)}")

    total_generated = 0
    total_skipped = 0

    # 外层循环：遍历阈值
    for thr_pct in PHENO_THRESHOLD_LIST:
        thr_tag = f"thr{int(thr_pct*100)}pct" if thr_pct > 0 else "no_threshold"

        print(f"\n{'='*70}")
        print(f"阈值: {thr_pct*100:.0f}% ({thr_tag})")
        print(f"{'='*70}")

        # 内层循环：遍历数据源组合
        for run_idx, cfg in enumerate(RUN_CONFIGS_07, 1):
            run_name = cfg['run_name']
            var_label = cfg['var_label']
            data_type = cfg['data_type']

            # 输出子目录：run_name/thr_tag/
            run_output_dir = OUTPUT_FIG_DIR / run_name / thr_tag
            run_output_dir.mkdir(parents=True, exist_ok=True)

            # 检查是否已有结果（看分解示意图是否齐全，4个类别）
            expected_files = [
                run_output_dir / f"decomposition_schematic_{var_label}_{cat}.png"
                for cat in ['sos_early_eos_late', 'sos_early_eos_early',
                            'sos_late_eos_late', 'sos_late_eos_early']
            ]
            if all(f.exists() for f in expected_files):
                print(f"\n  [跳过] {run_name}/{thr_tag} - 已有完整结果")
                total_skipped += 8
                continue

            print(f"\n{'#'*70}")
            print(f"## 组合 {run_idx}/{len(RUN_CONFIGS_07)}: {run_name} "
                  f"(标签: {var_label}, 阈值: {thr_pct*100:.0f}%)")
            print(f"{'#'*70}")

            # [1/4] 加载物候数据
            print(f"\n  [1/4] 加载物候数据 ({cfg['pheno_dir'].name})...")
            sos_stack, eos_stack, valid_years = load_phenology_all_years(
                years, pheno_dir=cfg['pheno_dir'], pheno_format=cfg['pheno_format'])
            print(f"    有效年份数: {len(valid_years)}")
            print(f"    数据形状: {sos_stack.shape}")

            # [2/4] 按物候异常分类年份
            print(f"\n  [2/4] 按物候异常分类年份（阈值={thr_pct*100:.0f}%）...")
            categories, category_year_counts, sos_mean, eos_mean = \
                classify_years_by_phenology(sos_stack, eos_stack, valid_years, mask,
                                             threshold_pct=thr_pct)

            category_stats = compute_category_statistics(categories, sos_stack, eos_stack, mask)
            for cat, stats in category_stats.items():
                print(f"    {LABELS[cat]}: 平均{stats['n_years_mean']:.1f}年/像元")

            # [3/4] 计算日数据
            print(f"\n  [3/4] 计算各类别的日{var_label}...")
            daily_data = compute_daily_data_by_category(
                valid_years, categories, mask,
                daily_dir=cfg['daily_dir'], daily_format=cfg['daily_format'],
                var_label=var_label, threshold_pct=thr_pct)

            # [4/4] 拟合曲线并绘图
            print(f"\n  [4/4] 拟合曲线并绘图...")

            x_daily = np.arange(1, 366)
            fitted_curves = {}
            pheno_params = {}

            for cat in ['baseline'] + list(categories.keys()):
                y = daily_data[cat]
                y_smooth = smooth_curve(y)

                valid_count = np.sum(np.isfinite(y))
                print(f"    {cat}: {valid_count} valid days")

                popt = fit_double_logistic(x_daily, y_smooth, data_type=data_type)
                if popt is not None:
                    y_fitted = double_logistic(x_daily, *popt)
                    fitted_curves[cat] = y_fitted
                    sos, eos, pos = extract_phenology_from_curve(y_fitted, x_daily)
                    pheno_params[cat] = {'SOS': sos, 'EOS': eos, 'POS': pos}
                    print(f"      拟合成功: SOS={sos:.1f}, POS={pos:.1f}, EOS={eos:.1f}")
                else:
                    fitted_curves[cat] = y_smooth
                    pheno_params[cat] = {'SOS': np.nan, 'EOS': np.nan, 'POS': np.nan}
                    print(f"      拟合失败！")

            baseline_curve = fitted_curves['baseline']
            baseline_pheno = pheno_params['baseline']

            file_fmt = f"phenology_composite_{var_label}_{{category}}.png"

            for category in categories.keys():
                # 曲线图
                output_path = run_output_dir / file_fmt.format(category=category)
                plot_single_category(
                    daily_data=daily_data,
                    category=category,
                    category_stats=category_stats,
                    output_path=output_path,
                    baseline_curve=baseline_curve,
                    baseline_pheno=baseline_pheno,
                    cat_curve=fitted_curves[category],
                    cat_pheno=pheno_params[category],
                    var_label=var_label
                )

                # 分解示意图
                decomp_path = run_output_dir / f"decomposition_schematic_{var_label}_{category}.png"
                plot_decomposition_schematic(
                    baseline_curve=baseline_curve,
                    cat_curve=fitted_curves[category],
                    baseline_pheno=baseline_pheno,
                    cat_pheno=pheno_params[category],
                    category=category,
                    output_path=decomp_path,
                    var_label=var_label
                )
                total_generated += 2

            # 打印物候参数对比
            print(f"\n  [{var_label}, {thr_tag}] 物候参数对比:")
            print(f"  {'类别':<30} {'SOS':>10} {'POS':>10} {'EOS':>10} {'LOS':>10}")
            print("  " + "-"*70)
            for cat in ['baseline'] + list(categories.keys()):
                if cat in pheno_params:
                    p = pheno_params[cat]
                    los = p['EOS'] - p['SOS'] if np.isfinite(p['EOS']) and np.isfinite(p['SOS']) else np.nan
                    print(f"  {LABELS[cat]:<30} {p['SOS']:>10.0f} {p['POS']:>10.0f} "
                          f"{p['EOS']:>10.0f} {los:>10.0f}")

            print(f"\n  输出目录: {run_output_dir}")

    # 汇总
    print(f"\n{'='*70}")
    print(f"全部完成！新生成 {total_generated} 张图，跳过 {total_skipped} 张已有结果")
    print(f"输出根目录: {OUTPUT_FIG_DIR}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
