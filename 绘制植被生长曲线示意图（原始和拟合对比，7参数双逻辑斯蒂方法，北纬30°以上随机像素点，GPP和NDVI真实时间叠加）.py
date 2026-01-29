# -*- coding: utf-8 -*-
"""
NDVI与GPP物候对比图（同像元双Y轴）
—— 使用真实DOY映射，对比NDVI和GPP的物候差异
"""

import os, re, time, numpy as np, rasterio
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit

# 配置中文字体显示
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

# ========== 配置 ==========
NDVI_FOLDER = r'I:\F\Data4\GIMMS_NDVI\GIMMS_NDVI_4'  # NDVI半月数据
GPP_FOLDER  = r'I:\F\Data4\GLASS_GPP\GLASS_GPP_8days_1'  # GPP 8天数据（已重采样到NDVI网格）
YEAR        = 2005
AMP_PERCENT = 0.20        # 动态阈值百分比
USE_RANDOM  = True        # True: 随机挑一个植被像元；False: 用下方 ROW/COL
ROW, COL    = 300, 500    # 指定像元时生效
RAND_SEED   = None        # None=每次运行自动生成新种子（在main中动态设置）
SG_WINDOW   = 15          # 平滑窗口
LAT_MIN     = 30.0        # 纬度下限

# ========== 稳定 sigmoid 与双 Logistic ==========
def stable_sigmoid(x, clip=30.0):
    x = np.clip(x, -clip, clip)
    return 1.0 / (1.0 + np.exp(-x))

def double_logistic(t, c1, c2, c3, r1, r2, t1, t2):
    rise = stable_sigmoid(r1 * (t - t1))
    fall = stable_sigmoid(r2 * (t - t2))
    return c1 + c2 * rise - c3 * fall

# ========== NDVI文件处理（真实DOY映射） ==========
def list_ndvi_files(folder, year):
    files = [f for f in os.listdir(folder) if f.endswith('.tif') and f'_{year}' in f]
    return [os.path.join(folder, f) for f in sorted(files)]

def ndvi_halfmonth_to_doy(fname):
    """NDVI半月期到真实DOY的映射"""
    HALFMONTH_DOY_MAP = {
        (1, 'a'): 8,    (1, 'b'): 24,   (2, 'a'): 39,   (2, 'b'): 53,
        (3, 'a'): 67,   (3, 'b'): 83,   (4, 'a'): 98,   (4, 'b'): 113,
        (5, 'a'): 128,  (5, 'b'): 144,  (6, 'a'): 159,  (6, 'b'): 174,
        (7, 'a'): 189,  (7, 'b'): 205,  (8, 'a'): 220,  (8, 'b'): 236,
        (9, 'a'): 251,  (9, 'b'): 266,  (10, 'a'): 281, (10, 'b'): 297,
        (11, 'a'): 312, (11, 'b'): 327, (12, 'a'): 342, (12, 'b'): 358,
    }
    m = re.search(r'_(\d{4})(\d{2})([ab])', os.path.basename(fname))
    if not m: return None
    month, half = int(m.group(2)), m.group(3)
    return HALFMONTH_DOY_MAP.get((month, half), None)

# ========== GPP文件处理 ==========
def list_gpp_files(folder, year):
    files = []
    for f in os.listdir(folder):
        if f.endswith('.tif') and f.startswith('GLASS_'):
            m = re.search(r'GLASS_(\d{4})\d{3}\.tif', f)
            if m and int(m.group(1)) == year:
                files.append(os.path.join(folder, f))
    return sorted(files)

def gpp_filename_to_doy(fname):
    """GPP文件名到DOY"""
    m = re.search(r'GLASS_\d{4}(\d{3})\.tif', os.path.basename(fname))
    return int(m.group(1)) if m else None

# ========== 通用读取函数 ==========
def read_pixel_series(files, row, col, doy_func, nodata_check):
    """读取指定像元的时间序列"""
    vals, doys = [], []
    for fp in files:
        with rasterio.open(fp) as src:
            v = src.read(1)[row, col].astype(float)
        if nodata_check(v):
            v = np.nan
        vals.append(v)
        doys.append(doy_func(fp))
    return np.array(doys, float), np.array(vals, float)

# ========== 随机选点 ==========
def pick_random_pixel(ndvi_files, gpp_files, seed=None, tries=200, lat_min=LAT_MIN):
    """
    选择同时满足NDVI和GPP条件的像元
    """
    rng = np.random.default_rng(seed)

    with rasterio.open(ndvi_files[0]) as src:
        H, W = src.height, src.width
        transform = src.transform

        if lat_min is not None:
            lat_by_row = np.array([(transform * (0, r))[1] for r in range(H)], dtype=float)
            valid_rows = np.where(lat_by_row >= lat_min)[0]
        else:
            valid_rows = np.arange(H)

    if valid_rows.size == 0:
        return None, None

    for attempt in range(tries):
        r = int(rng.choice(valid_rows))
        c = int(rng.integers(0, W))

        # 检查NDVI
        doys_ndvi, vals_ndvi = read_pixel_series(
            ndvi_files, r, c, ndvi_halfmonth_to_doy,
            lambda v: v == -9999 or v < 0
        )
        ndvi_valid = np.isfinite(vals_ndvi).sum() >= max(3, int(0.7*len(vals_ndvi)))
        ndvi_amp = np.nanmax(vals_ndvi) - np.nanmin(vals_ndvi) if ndvi_valid else 0

        # 检查GPP
        doys_gpp, vals_gpp = read_pixel_series(
            gpp_files, r, c, gpp_filename_to_doy,
            lambda v: v == -9999 or v < 0 or not np.isfinite(v)
        )
        gpp_valid = np.isfinite(vals_gpp).sum() >= max(3, int(0.7*len(vals_gpp)))
        gpp_amp = np.nanmax(vals_gpp) - np.nanmin(vals_gpp) if gpp_valid else 0

        # 同时满足条件
        if ndvi_valid and gpp_valid and ndvi_amp >= 0.10 and gpp_amp >= 0.3:
            print(f"找到合格像元：Row={r}, Col={c} (尝试了{attempt+1}次)")
            print(f"  NDVI振幅={ndvi_amp:.3f}, GPP振幅={gpp_amp:.2f}")
            return r, c

        if (attempt + 1) % 50 == 0:
            print(f"正在搜索... 已尝试 {attempt+1}/{tries} 次")

    return None, None

# ========== 动态阈值提取物候 ==========
def extract_phenology(curve, x_daily, percent=0.20):
    """提取SOS/EOS/POS"""
    y = np.asarray(curve, float)
    smin, smax = np.nanmin(y), np.nanmax(y)
    if not np.isfinite(smin) or smax - smin < 1e-6:
        return np.nan, np.nan, np.nan, np.nan

    yn = (y - smin) / (smax - smin)
    peak_idx = int(np.nanargmax(yn))
    left = yn[:peak_idx+1]
    right = yn[peak_idx:]

    left_min_idx = int(np.nanargmin(left))
    right_min_idx = peak_idx + int(np.nanargmin(right))
    minL, minR = yn[left_min_idx], yn[right_min_idx]
    peak = yn[peak_idx]

    thrL = (peak - minL) * percent + minL
    thrR = (peak - minR) * percent + minR

    sos_day = np.nan
    for i in range(max(left_min_idx, 0), peak_idx+1):
        if yn[i] >= thrL:
            sos_day = x_daily[i]
            break

    eos_day = np.nan
    for i in range(peak_idx, right_min_idx+1):
        if yn[i] >= thrR:
            eos_day = x_daily[i]

    pos_doy = x_daily[peak_idx]
    pos_val = y[peak_idx]
    return sos_day, eos_day, pos_doy, pos_val

# ========== 拟合双Logistic ==========
def fit_double_logistic(x_valid, y_valid, data_type='ndvi'):
    """拟合双Logistic，根据数据类型调整边界"""
    data_min, data_max = float(np.min(y_valid)), float(np.max(y_valid))
    dr = data_max - data_min

    p0 = [data_min, dr/2, dr/2, 0.1, 0.1, 100, 250]

    if data_type == 'ndvi':
        bounds_lower = [-0.5, 0, 0, 0, 0,   0,  50]
        bounds_upper = [ 1.5, 2, 2, 2, 2, 183, 366]
    else:  # gpp
        bounds_lower = [data_min-dr, 0,   0,   0, 0,   0,  50]
        bounds_upper = [data_max+dr, dr*3, dr*3, 2, 2, 183, 366]

    popt = None
    for it in range(3):
        try:
            popt, _ = curve_fit(double_logistic, x_valid, y_valid,
                                p0=p0, bounds=(bounds_lower, bounds_upper),
                                maxfev=10000)
            break
        except Exception:
            p0 = [data_min-0.02*dr, 0.9*dr, 0.9*dr, 0.08, 0.08, 110, 240]
    return popt

# ========== 主流程 ==========
def main():
    ndvi_files = list_ndvi_files(NDVI_FOLDER, YEAR)
    gpp_files = list_gpp_files(GPP_FOLDER, YEAR)

    if not ndvi_files or not gpp_files:
        raise RuntimeError("未找到NDVI或GPP文件")

    print(f"NDVI文件: {len(ndvi_files)}个, GPP文件: {len(gpp_files)}个")

    # 动态生成随机种子（每次运行不同，但同一次运行中NDVI和GPP使用相同种子）
    if USE_RANDOM and RAND_SEED is None:
        import time
        dynamic_seed = int(time.time() * 1000) % (2**32)
        print(f"本次运行使用的随机种子: {dynamic_seed}")
    else:
        dynamic_seed = RAND_SEED

    # 选择像元
    if USE_RANDOM:
        r, c = pick_random_pixel(ndvi_files, gpp_files, seed=dynamic_seed)
        if r is None:
            raise RuntimeError("未找到合格像元")
    else:
        r, c = ROW, COL
        print(f"使用指定像元：Row={r}, Col={c}")

    # 读取NDVI数据
    doys_ndvi, ndvi_raw = read_pixel_series(
        ndvi_files, r, c, ndvi_halfmonth_to_doy,
        lambda v: v == -9999 or v < 0
    )

    # 读取GPP数据
    doys_gpp, gpp_raw = read_pixel_series(
        gpp_files, r, c, gpp_filename_to_doy,
        lambda v: v == -9999 or v < 0 or not np.isfinite(v)
    )

    # 插值到日尺度
    x_daily = np.arange(1, 366, dtype=float)

    # NDVI处理
    mask_ndvi = np.isfinite(ndvi_raw)
    ndvi_daily = np.interp(x_daily, doys_ndvi[mask_ndvi], ndvi_raw[mask_ndvi])
    if SG_WINDOW and SG_WINDOW >= 3:
        ndvi_daily_s = savgol_filter(ndvi_daily, SG_WINDOW, 3, mode='nearest')
    else:
        ndvi_daily_s = ndvi_daily

    # GPP处理
    mask_gpp = np.isfinite(gpp_raw)
    gpp_daily = np.interp(x_daily, doys_gpp[mask_gpp], gpp_raw[mask_gpp])
    if SG_WINDOW and SG_WINDOW >= 3:
        gpp_daily_s = savgol_filter(gpp_daily, SG_WINDOW, 3, mode='nearest')
    else:
        gpp_daily_s = gpp_daily

    # 拟合
    popt_ndvi = fit_double_logistic(x_daily, ndvi_daily_s, 'ndvi')
    popt_gpp = fit_double_logistic(x_daily, gpp_daily_s, 'gpp')

    if popt_ndvi is None or popt_gpp is None:
        raise RuntimeError("拟合失败")

    fitted_ndvi = double_logistic(x_daily, *popt_ndvi)
    fitted_gpp = double_logistic(x_daily, *popt_gpp)

    # 提取物候参数
    sos_ndvi, eos_ndvi, pos_ndvi, val_ndvi = extract_phenology(fitted_ndvi, x_daily, AMP_PERCENT)
    sos_gpp, eos_gpp, pos_gpp, val_gpp = extract_phenology(fitted_gpp, x_daily, AMP_PERCENT)

    los_ndvi = eos_ndvi - sos_ndvi if np.isfinite(sos_ndvi) and np.isfinite(eos_ndvi) else np.nan
    los_gpp = eos_gpp - sos_gpp if np.isfinite(sos_gpp) and np.isfinite(eos_gpp) else np.nan

    # ========== 绘图（双Y轴） ==========
    fig, ax1 = plt.subplots(figsize=(12, 6))

    # 左Y轴：NDVI
    color_ndvi = 'tab:green'
    ax1.set_xlabel('Time (Julian Day)', fontsize=12)
    ax1.set_ylabel('NDVI', color=color_ndvi, fontsize=12)
    ax1.plot(doys_ndvi, ndvi_raw, 'o', color=color_ndvi, alpha=0.5,
             markersize=5, label='NDVI (half-month)')
    ax1.plot(x_daily, fitted_ndvi, '-', color=color_ndvi, linewidth=2.5,
             label='NDVI fitted curve')
    ax1.tick_params(axis='y', labelcolor=color_ndvi)
    ax1.set_ylim(bottom=-0.05)

    # 右Y轴：GPP
    ax2 = ax1.twinx()
    color_gpp = 'tab:orange'
    ax2.set_ylabel('GPP (gC m$^{-2}$ day$^{-1}$)', color=color_gpp, fontsize=12)
    ax2.plot(doys_gpp, gpp_raw, 's', color=color_gpp, alpha=0.5,
             markersize=4, label='GPP (8-day)')
    ax2.plot(x_daily, fitted_gpp, '-', color=color_gpp, linewidth=2.5,
             label='GPP fitted curve')
    ax2.tick_params(axis='y', labelcolor=color_gpp)
    ax2.set_ylim(bottom=-0.1)

    # 标记物候点（NDVI - 左Y轴）
    if np.isfinite(sos_ndvi):
        ax1.axvline(sos_ndvi, color=color_ndvi, linestyle='--', alpha=0.6, linewidth=1)
        ax1.text(sos_ndvi, ax1.get_ylim()[1]*0.95, f'NDVI SOS\n{int(sos_ndvi)}',
                ha='center', fontsize=9, color=color_ndvi, fontweight='bold')

    if np.isfinite(eos_ndvi):
        ax1.axvline(eos_ndvi, color=color_ndvi, linestyle='--', alpha=0.6, linewidth=1)
        ax1.text(eos_ndvi, ax1.get_ylim()[1]*0.95, f'NDVI EOS\n{int(eos_ndvi)}',
                ha='center', fontsize=9, color=color_ndvi, fontweight='bold')

    # 标记物候点（GPP - 右Y轴）
    if np.isfinite(sos_gpp):
        ax2.axvline(sos_gpp, color=color_gpp, linestyle=':', alpha=0.6, linewidth=1)
        ax2.text(sos_gpp, ax2.get_ylim()[1]*0.85, f'GPP SOS\n{int(sos_gpp)}',
                ha='center', fontsize=9, color=color_gpp, fontweight='bold')

    if np.isfinite(eos_gpp):
        ax2.axvline(eos_gpp, color=color_gpp, linestyle=':', alpha=0.6, linewidth=1)
        ax2.text(eos_gpp, ax2.get_ylim()[1]*0.85, f'GPP EOS\n{int(eos_gpp)}',
                ha='center', fontsize=9, color=color_gpp, fontweight='bold')

    # 图例
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=9)

    plt.title(f'NDVI vs GPP Phenology Comparison — Row {r}, Col {c} — {YEAR}',
              fontsize=13, fontweight='bold', pad=15)
    plt.grid(alpha=0.3, linestyle='--')
    fig.tight_layout()
    plt.show()

    # ========== 打印对比结果 ==========
    print("\n" + "="*70)
    print(f"NDVI与GPP物候对比 ({YEAR}年) — Row {r}, Col {c}")
    print("="*70)
    print(f"{'参数':<15} {'NDVI':>12} {'GPP':>12} {'差异(GPP-NDVI)':>18}")
    print("-"*70)
    print(f"{'SOS (开始日)':<15} {f'第{sos_ndvi:.0f}天':>12} {f'第{sos_gpp:.0f}天':>12} "
          f"{f'{sos_gpp-sos_ndvi:+.0f}天':>18}")
    print(f"{'EOS (结束日)':<15} {f'第{eos_ndvi:.0f}天':>12} {f'第{eos_gpp:.0f}天':>12} "
          f"{f'{eos_gpp-eos_ndvi:+.0f}天':>18}")
    print(f"{'POS (峰值日)':<15} {f'第{pos_ndvi:.0f}天':>12} {f'第{pos_gpp:.0f}天':>12} "
          f"{f'{pos_gpp-pos_ndvi:+.0f}天':>18}")
    print(f"{'LOS (生长季长)':<15} {f'{los_ndvi:.0f}天':>12} {f'{los_gpp:.0f}天':>12} "
          f"{f'{los_gpp-los_ndvi:+.0f}天':>18}")
    print("-"*70)
    print(f"{'峰值':<15} {f'{val_ndvi:.4f}':>12} {f'{val_gpp:.2f} gC/m²/d':>12}")
    print("="*70)

    # 物候差异分析
    print("\n物候差异分析:")
    if np.isfinite(sos_ndvi) and np.isfinite(sos_gpp):
        diff = sos_gpp - sos_ndvi
        if abs(diff) < 5:
            print(f"  ✓ SOS基本一致 (差异{diff:+.0f}天)")
        elif diff > 0:
            print(f"  ⚠ GPP的生长季开始晚于NDVI {diff:.0f}天")
        else:
            print(f"  ⚠ GPP的生长季开始早于NDVI {abs(diff):.0f}天")

    if np.isfinite(eos_ndvi) and np.isfinite(eos_gpp):
        diff = eos_gpp - eos_ndvi
        if abs(diff) < 5:
            print(f"  ✓ EOS基本一致 (差异{diff:+.0f}天)")
        elif diff > 0:
            print(f"  ⚠ GPP的生长季结束晚于NDVI {diff:.0f}天")
        else:
            print(f"  ⚠ GPP的生长季结束早于NDVI {abs(diff):.0f}天")

if __name__ == "__main__":
    main()
