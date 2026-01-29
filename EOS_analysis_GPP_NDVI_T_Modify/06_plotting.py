#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module 06: 统一绘图入口
04c 结果：偏相关/回归/趋势图使用“偏相关系数与显著性图”模板样式
05b/05c/05d 结果：SEM 路径图使用"EOS-NDVI全路径图像"模板样式
"""

from dataclasses import dataclass
import argparse
from pathlib import Path
import os
import colorsys

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.path as mpath
import matplotlib.patches as patches
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from osgeo import gdal


# ==================== 运行模式配置 ====================
# 支持环境变量控制是否覆盖已存在的图片
# 方式1: 设置 WANG_OVERWRITE=true 强制重新绘制所有图片
# 方式2: 设置 WANG_RUN_MODE=overwrite (与pipeline保持一致)
_overwrite_flag = os.getenv("WANG_OVERWRITE", "false").lower() in ("true", "1", "yes")
_run_mode = os.getenv("WANG_RUN_MODE", "skip").lower()
OVERWRITE = _overwrite_flag or (_run_mode == "overwrite")


def should_plot(out_path):
    """检查是否应该绘图（如果文件已存在且非覆盖模式则跳过）"""
    out_path = Path(out_path)
    if OVERWRITE:
        return True
    if out_path.exists():
        print(f"  [skip] 已存在: {out_path.name}")
        return False
    return True


# 导入配置中的变量名和输出根目录
from _config import MIDDLE_VAR_NAME, OUTPUT_ROOT, FIGURES_DIR

# ==================== 路径与全局配置 ====================
ANALYSIS_DIR = OUTPUT_ROOT
OUTPUT_FIG_DIR = FIGURES_DIR  # 使用_config.py中的统一配置
OUTPUT_FIG_DIR.mkdir(parents=True, exist_ok=True)

SEM_FALLBACK_DIRS = [
    Path(r"D:\claude-project\结果分析"),
    Path("/mnt/d/claude-project/结果分析"),
]

FONT_SIZE = 9
NODATA_IN = -9999.0
PIXEL_SIG_FRAC_THRESHOLD = 0.5
USE_NEIGHBORHOOD_FILTER = True
NEIGHBORHOOD_THRESHOLD = 0.5
PLOT_ONLY_DEFAULT = "all"  # "all", "04c", "05b", "05c", "05d"


# ==================== 04c 模板样式函数（偏相关图模板） ====================
def _set_corr_style():
    mpl.rcParams.update({
        "font.family": "Arial",
        "font.size": FONT_SIZE,
        "axes.titlesize": FONT_SIZE,
        "axes.labelsize": FONT_SIZE,
        "xtick.labelsize": FONT_SIZE,
        "ytick.labelsize": FONT_SIZE,
        "axes.linewidth": 0.3,
        "figure.dpi": 300,
        "savefig.bbox": "tight",
        "mathtext.fontset": "stixsans",
        "mathtext.default": "regular",
    })
    plt.rcParams["agg.path.chunksize"] = 10000


def _load_pvalue_data(pvalue_file, factor=1):
    if pvalue_file is None or (not Path(pvalue_file).exists()):
        return None
    try:
        ds = gdal.Open(str(pvalue_file), gdal.GA_ReadOnly)
        if ds is None:
            print(f"警告: 无法打开 p 值文件 {pvalue_file}")
            return None
        p_arr = ds.ReadAsArray().astype(float)
        nodata = ds.GetRasterBand(1).GetNoDataValue()
        if nodata is not None:
            p_arr[p_arr == nodata] = np.nan
        if factor > 1:
            p_arr = p_arr[::factor, ::factor]
        return p_arr
    except Exception as exc:
        print(f"警告: 处理 p 值文件时出错 {pvalue_file}: {exc}")
        return None


def _apply_neighborhood_filter(sig_mask, threshold_ratio=0.5):
    try:
        from scipy import ndimage
    except ImportError:
        print("警告：scipy 未安装，跳过邻域过滤")
        return sig_mask

    structure = np.ones((3, 3), dtype=int)
    sig_count = ndimage.generic_filter(
        sig_mask.astype(int),
        np.sum,
        footprint=structure,
        mode="constant",
        cval=0,
    )
    valid_mask = np.ones_like(sig_mask, dtype=int)
    valid_count = ndimage.generic_filter(
        valid_mask,
        np.sum,
        footprint=structure,
        mode="constant",
        cval=0,
    )

    with np.errstate(divide="ignore", invalid="ignore"):
        sig_ratio = sig_count / valid_count
        sig_ratio[valid_count == 0] = 0

    return sig_mask & (sig_ratio >= threshold_ratio)


def _add_significance_overlay(ax, pvalue_arr, extent, proj_pc):
    if pvalue_arr is None:
        return

    nrows, ncols = pvalue_arr.shape
    lon_min, lon_max, lat_min, lat_max = extent
    lon_res = (lon_max - lon_min) / ncols
    lat_res = (lat_max - lat_min) / nrows

    lon_centers = lon_min + (np.arange(ncols) + 0.5) * lon_res
    lat_centers = lat_max - (np.arange(nrows) + 0.5) * lat_res
    xx, yy = np.meshgrid(lon_centers, lat_centers)

    sig_mask = (pvalue_arr < 0.05) & (~np.isnan(pvalue_arr))
    if USE_NEIGHBORHOOD_FILTER and np.any(sig_mask):
        sig_mask = _apply_neighborhood_filter(sig_mask, threshold_ratio=NEIGHBORHOOD_THRESHOLD)

    if np.any(sig_mask):
        sig_x = xx[sig_mask]
        sig_y = yy[sig_mask]
        marker_w = lon_res
        marker_h = lat_res
        for i in range(len(sig_x)):
            rect = patches.Rectangle(
                (sig_x[i] - marker_w / 2, sig_y[i] - marker_h / 2),
                marker_w,
                marker_h,
                linewidth=0,
                edgecolor="none",
                facecolor="black",
                alpha=0.4,
                transform=proj_pc,
                zorder=10,
            )
            ax.add_patch(rect)


def _compute_axis_scaling(ticks):
    ticks = np.asarray(ticks)
    if ticks.size == 0:
        return 1.0, 0
    max_abs = np.nanmax(np.abs(ticks))
    if max_abs == 0 or np.isnan(max_abs):
        return 1.0, 0
    exp = int(np.floor(np.log10(max_abs)))
    if abs(exp) <= 1:
        return 1.0, 0
    return 10.0**exp, exp


def _format_tick_label(value, decimals=1):
    if np.isnan(value):
        return ""
    threshold = 0.5 * (10 ** (-decimals))
    if abs(value) < threshold:
        return "0"
    fmt = f"{{:.{decimals}f}}"
    s = fmt.format(value)
    return s.rstrip("0").rstrip(".")


def _add_polar_lon_labels(ax, lon_values, lat_label=35, proj_src=ccrs.PlateCarree()):
    """外圈经度标注，沿圆环切线方向（与模板一致）"""
    proj = ax.projection
    data_to_axes = ax.transAxes.inverted().transform
    cx, cy = 0.5, 0.5
    r_ring = 0.5
    margin = 0.02
    r_text = r_ring + margin

    for lon in lon_values:
        # 下方不标，避免和色条/直方图冲突（与模板一致）
        if -60 <= lon <= 60:
            continue

        if lon == 0:
            label = "0°"
        elif lon == 180 or lon == -180:
            label = "180°"
        else:
            hemi = 'E' if lon > 0 else 'W'
            label = f"{abs(int(lon))}°{hemi}"

        x_data, y_data = proj.transform_point(lon, lat_label, proj_src)
        x_disp, y_disp = ax.transData.transform((x_data, y_data))
        x_axes, y_axes = data_to_axes((x_disp, y_disp))

        vx = x_axes - cx
        vy = y_axes - cy
        r = np.hypot(vx, vy)
        if r == 0:
            continue

        scale = r_text / r
        x_text = cx + vx * scale
        y_text = cy + vy * scale

        # 去掉最顶部的 180°（与模板一致）
        if label == "180°" and y_text > 0.95:
            continue

        theta = np.arctan2(vy, vx)
        angle_deg = np.degrees(theta) + 90
        if angle_deg > 90:
            angle_deg -= 180
        if angle_deg < -90:
            angle_deg += 180

        ax.text(
            x_text,
            y_text,
            label,
            transform=ax.transAxes,
            ha="center",
            va="center",
            rotation=angle_deg,
            rotation_mode="anchor",
            fontsize=FONT_SIZE,
        )


def _add_lat_labels(ax, lat_values=(45, 60, 75), lon_for_labels=-130, proj_src=ccrs.PlateCarree()):
    for lat in lat_values:
        ax.text(
            lon_for_labels,
            lat,
            f"{int(lat)}°N",
            transform=proj_src,
            ha="left",
            va="center",
            fontsize=FONT_SIZE,
        )


def _get_pair_short_labels(var_label):
    """
    从文件名中解析出 X ~ Y（与模板一致）
    例如：EOS_vs_SpringET_controls_... → EOS ~ ET
    """
    base = var_label
    if base.endswith("_pcorr_r"):
        base = base[:-8]

    if "_vs_" in base:
        left, right = base.split("_vs_", 1)
        x_raw = left
        if "_controls_" in right:
            y_raw = right.split("_controls_", 1)[0]
        else:
            y_raw = right
    else:
        # 兜底：无法解析就用原始文件名
        return var_label, ""

    # 解析 Y：去掉前面的季节 / Annual 前缀
    def map_y(raw):
        seasons = ["Spring", "Summer", "Autumn", "Winter", "Annual"]
        raw2 = raw
        for s in seasons:
            if raw2.startswith(s):
                raw2 = raw2[len(s):]
                break
        # 现在 raw2 应该类似 "ET", "T", "E", "GPP"
        return raw2 or raw

    # 解析 X
    def map_x(raw):
        # 物候：EOS / POS / deltaEOS
        if raw in ["EOS", "POS"]:
            return raw
        if raw in ["deltaEOS", "DeltaEOS"]:
            return "ΔEOS"
        # 气象：Spring_Tem / Summer_Pre / Autumn_DSW
        parts = raw.split("_")
        if len(parts) == 2 and parts[1] in ["Tem", "Pre", "DSW"]:
            return parts[1]
        seasons = ["Spring", "Summer", "Autumn", "Winter", "Annual"]
        for s in seasons:
            if raw.startswith(s):
                r2 = raw[len(s):]
                if r2 in ["Tem", "Pre", "DSW"]:
                    return r2
        return raw

    x_short = map_x(x_raw)
    y_short = map_y(y_raw)
    return x_short, y_short


def _make_bright_diverging_cmap():
    colors = [
        (95 / 255, 175 / 255, 250 / 255),
        (175 / 255, 215 / 255, 255 / 255),
        (1, 1, 1),
        (1, 220 / 255, 185 / 255),
        (1, 160 / 255, 115 / 255),
    ]
    brighten = 0.9
    for i, (r, g, b) in enumerate(colors):
        h, l, s = colorsys.rgb_to_hls(r, g, b)
        l = max(0, min(1, l * brighten))
        colors[i] = colorsys.hls_to_rgb(h, l, s)
    return mpl.colors.LinearSegmentedColormap.from_list("BrightDiv", colors, N=256)


def plot_corr_style_map(r_file, p_file, out_figure, var_label,
                        cbar_title=None, value_range=None, auto_range=None):
    if not should_plot(out_figure):
        return
    if not Path(r_file).exists():
        print(f"相关系数文件不存在，跳过: {r_file}")
        return

    _set_corr_style()
    ds = gdal.Open(str(r_file), gdal.GA_ReadOnly)
    if ds is None:
        print(f"GDAL 无法打开文件 {r_file}，跳过。")
        return

    arr = ds.ReadAsArray().astype(float)
    gt = ds.GetGeoTransform()
    nodata = ds.GetRasterBand(1).GetNoDataValue()
    if nodata is not None:
        arr[arr == nodata] = np.nan

    nrows, ncols = arr.shape
    extent = [gt[0], gt[0] + gt[1] * ncols, gt[3] + gt[5] * nrows, gt[3]]

    max_pixels = 10_000_000
    factor = int(np.ceil(np.sqrt((nrows * ncols) / max_pixels)))
    if factor > 1:
        arr = arr[::factor, ::factor]

    p_arr = _load_pvalue_data(p_file, factor=factor)
    valid_data = arr[~np.isnan(arr)]
    if valid_data.size == 0:
        print("全部为 NaN，跳过此图。")
        return

    if value_range is None:
        if auto_range == "symmetric":
            abs_max = np.nanpercentile(np.abs(valid_data), 98)
            if not np.isfinite(abs_max) or abs_max == 0:
                abs_max = 1.0
            value_range = (-abs_max, abs_max)
        else:
            value_range = (-1.0, 1.0)

    boundaries = np.linspace(value_range[0], value_range[1], 9)
    cmap = _make_bright_diverging_cmap()
    norm = mpl.colors.BoundaryNorm(boundaries, cmap.N)

    proj_map = ccrs.NorthPolarStereo()
    proj_pc = ccrs.PlateCarree()

    fig = plt.figure(figsize=(3.2, 3.8))
    ax = plt.axes(projection=proj_map)
    ax.set_extent([-180, 180, 30, 90], crs=proj_pc)

    ax.add_feature(cfeature.LAND.with_scale("110m"), facecolor="none", edgecolor="k", linewidth=0.3, zorder=1)
    ax.add_feature(cfeature.COASTLINE.with_scale("110m"), linewidth=0.3, zorder=1)

    theta = np.linspace(0, 2 * np.pi, 360)
    verts = np.column_stack([0.5 + 0.5 * np.sin(theta), 0.5 + 0.5 * np.cos(theta)])
    ax.set_boundary(mpath.Path(verts, closed=True), transform=ax.transAxes)

    mesh = ax.imshow(arr, origin="upper", extent=extent, transform=proj_pc, cmap=cmap, norm=norm,
                     interpolation="nearest", zorder=2)

    _add_significance_overlay(ax, p_arr, extent, proj_pc)

    gl = ax.gridlines(crs=proj_pc, draw_labels=False, linewidth=0.4, color="grey",
                      linestyle=(0, (5, 2)), alpha=0.8)
    gl.xlocator = mticker.FixedLocator(np.arange(-180, 181, 30))
    gl.ylocator = mticker.FixedLocator([45, 60, 75])

    _add_polar_lon_labels(ax, np.arange(-180, 181, 30), lat_label=35)
    _add_lat_labels(ax, lat_values=(45, 60, 75), lon_for_labels=-130)

    ticks = boundaries[::2]
    scale_factor, _ = _compute_axis_scaling(ticks)
    cax = inset_axes(ax, width="52%", height="3%", loc="lower right",
                     bbox_to_anchor=(0.0, -0.08, 1, 1),
                     bbox_transform=ax.transAxes, borderpad=0)
    cbar = plt.colorbar(mesh, cax=cax, orientation="horizontal", ticks=ticks)
    scaled_values = [t / scale_factor for t in ticks]
    scaled_ticklabels = [_format_tick_label(v, decimals=1) for v in scaled_values]
    cbar.ax.set_xticklabels(scaled_ticklabels)
    cax.tick_params(length=0, width=0, pad=1)
    for spine in cax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.3)

    if cbar_title is None:
        x_short, y_short = _get_pair_short_labels(var_label)
        cbar_title = f"Partial r, {x_short} ~ {y_short}"
    cbar.ax.text(0.5, 1.03, cbar_title, transform=cbar.ax.transAxes,
                 ha="center", va="bottom", fontsize=FONT_SIZE)

    hist_ax = inset_axes(ax, width="34%", height="34%", loc="lower left",
                         bbox_to_anchor=(0.05, -0.06, 1, 1),
                         bbox_transform=ax.transAxes, borderpad=0)
    bins = boundaries
    hist, _ = np.histogram(valid_data, bins=bins)
    percent = hist / hist.sum() * 100
    centers = (bins[:-1] + bins[1:]) / 2
    bar_width = centers[1] - centers[0]

    # 使用色带颜色（与模板一致）
    hist_ax.bar(centers / scale_factor, percent,
                width=(bar_width / scale_factor) * 0.9,
                color=cmap(norm(centers)), edgecolor='none', align='center')

    # 设置直方图刻度
    hist_ticks = ticks
    hist_tick_values = [t / scale_factor for t in hist_ticks]
    hist_ax.set_xticks(hist_tick_values)
    hist_ax.set_xticklabels(
        [_format_tick_label(v, decimals=1) for v in hist_tick_values]
    )

    # 设置Y轴范围和刻度（与模板一致）
    hist_ax.set_ylim(0, max(60, np.nanmax(percent) * 1.2))
    hist_ax.set_yticks([15, 30, 45, 60])
    hist_ax.tick_params(direction='out', length=2, width=.3, pad=1)
    for spine in hist_ax.spines.values():
        spine.set_linewidth(.3)

    hist_ax.set_ylabel('Frequency (%)', labelpad=3, fontsize=FONT_SIZE)
    hist_ax.axhline(0, color='black', linewidth=.3)

    # P / N 比例（正/负相关 + 显著比例）（与模板一致）
    total_pixels = valid_data.size
    pos_corr = np.sum(valid_data >= 0)
    neg_corr = np.sum(valid_data < 0)
    pos_percent = pos_corr / total_pixels * 100
    neg_percent = neg_corr / total_pixels * 100

    sig_pos_percent = 0.0
    sig_neg_percent = 0.0

    if p_arr is not None:
        valid_mask = ~np.isnan(p_arr) & ~np.isnan(arr)
        valid_p = p_arr[valid_mask]
        valid_r = arr[valid_mask]
        if valid_p.size > 0:
            sig_mask = valid_p < 0.05
            sig_pos_pixels = np.sum(sig_mask & (valid_r > 0))
            sig_neg_pixels = np.sum(sig_mask & (valid_r < 0))
            sig_pos_percent = sig_pos_pixels / total_pixels * 100
            sig_neg_percent = sig_neg_pixels / total_pixels * 100

    pos_text = f'P: {pos_percent:.1f}% ({sig_pos_percent:.1f}%)'
    neg_text = f'N: {neg_percent:.1f}% ({sig_neg_percent:.1f}%)'

    ax.text(0.5, 0.93, pos_text, transform=ax.transAxes,
            ha='center', va='center', fontsize=FONT_SIZE, color='black')
    ax.text(0.5, 0.88, neg_text, transform=ax.transAxes,
            ha='center', va='center', fontsize=FONT_SIZE, color='black')

    # 布局调整（与模板一致）
    plt.subplots_adjust(bottom=0.12, top=0.98, left=0.03, right=0.97)

    out_figure.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_figure, format="tiff", dpi=600, bbox_inches="tight")
    plt.close()


# ==================== 05 模板样式函数（SEM 全路径图） ====================
def _set_sem_style():
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["axes.unicode_minus"] = False


def _load_sem_params(param_file):
    df = pd.read_csv(param_file)
    coef_map = {}

    if "op" in df.columns and "lhs" in df.columns:
        for _, row in df.iterrows():
            if row.get("op") != "~":
                continue
            lhs = row.get("lhs")
            rhs = row.get("rhs")
            coef = row.get("std.all")
            if not np.isfinite(coef):
                coef = row.get("est")
            pval = row.get("pvalue")
            try:
                coef = float(coef)
                pval = float(pval) if np.isfinite(pval) else np.nan
            except (TypeError, ValueError):
                continue
            sig = bool(np.isfinite(pval) and pval < 0.05)
            coef_map[(lhs, rhs)] = {"coef": coef, "sig": sig}
        return coef_map, "lavaan"

    if "path" in df.columns:
        for _, row in df.iterrows():
            path = row.get("path")
            if not isinstance(path, str) or "~" not in path:
                continue
            lhs, rhs = path.split("~", 1)
            coef = row.get("mean")
            sig_frac = row.get("sig_frac")
            try:
                coef = float(coef)
            except (TypeError, ValueError):
                continue
            sig = False
            boot_low = row.get("boot_ci_low", np.nan)
            boot_high = row.get("boot_ci_high", np.nan)
            try:
                boot_low = float(boot_low)
            except (TypeError, ValueError):
                boot_low = np.nan
            try:
                boot_high = float(boot_high)
            except (TypeError, ValueError):
                boot_high = np.nan
            if np.isfinite(boot_low) and np.isfinite(boot_high):
                sig = not (boot_low <= 0 <= boot_high)
            elif np.isfinite(sig_frac):
                sig = sig_frac >= PIXEL_SIG_FRAC_THRESHOLD
            coef_map[(lhs, rhs)] = {"coef": coef, "sig": sig}
        return coef_map, "pixelwise"

    raise ValueError(f"无法识别参数文件格式: {param_file}")


def _load_r2(r2_file):
    df = pd.read_csv(r2_file)
    r2_map = {}
    if "variable" in df.columns:
        for _, row in df.iterrows():
            try:
                r2_map[str(row["variable"])] = float(row["R2"])
            except (TypeError, ValueError):
                continue
    else:
        for _, row in df.iterrows():
            if len(row) < 2:
                continue
            try:
                r2_map[str(row[0])] = float(row[1])
            except (TypeError, ValueError):
                continue
    return r2_map


def _build_path(coef_map, lhs, rhs):
    item = coef_map.get((lhs, rhs))
    if not item:
        return None
    coef = item.get("coef")
    sig = item.get("sig", False)
    if not np.isfinite(coef):
        return None
    return {"coef": float(coef), "sig": bool(sig)}


def plot_sem_full_path(param_file, r2_file, out_png, title_text,
                       middle_var, outcome_var, middle_label, outcome_label,
                       bar_direct_label, bar_indirect_label, bar_ylabel):
    if not should_plot(out_png):
        return
    if not param_file.exists() or not r2_file.exists():
        print(f"SEM文件缺失，跳过: {param_file}")
        return

    _set_sem_style()
    coef_map, mode = _load_sem_params(param_file)
    r2_map = _load_r2(r2_file)

    r2_values = {
        "EOS": r2_map.get("EOS", np.nan),
        "NDVI": r2_map.get(middle_var, np.nan),
        "TRATE": r2_map.get(outcome_var, np.nan),
    }

    display_labels = {
        "NDVI": middle_label,
        "TRATE": outcome_label,
    }

    direct_paths = {}
    for key, (lhs, rhs) in {
        "EOS~P_pre": ("EOS", "P_pre"),
        "EOS~T_pre": ("EOS", "T_pre"),
        "EOS~SW_pre": ("EOS", "SW_pre"),
        "EOS~P_season": ("EOS", "P_season"),
        "EOS~T_season": ("EOS", "T_season"),
        "EOS~SW_season": ("EOS", "SW_season"),
        f"{middle_var}~EOS": (middle_var, "EOS"),
        f"{middle_var}~P_pre": (middle_var, "P_pre"),
        f"{middle_var}~T_pre": (middle_var, "T_pre"),
        f"{middle_var}~SW_pre": (middle_var, "SW_pre"),
        f"{middle_var}~P_season": (middle_var, "P_season"),
        f"{middle_var}~T_season": (middle_var, "T_season"),
        f"{middle_var}~SW_season": (middle_var, "SW_season"),
        f"{outcome_var}~EOS": (outcome_var, "EOS"),
        f"{outcome_var}~{middle_var}": (outcome_var, middle_var),
        f"{outcome_var}~P_pre": (outcome_var, "P_pre"),
        f"{outcome_var}~T_pre": (outcome_var, "T_pre"),
        f"{outcome_var}~SW_pre": (outcome_var, "SW_pre"),
        f"{outcome_var}~P_season": (outcome_var, "P_season"),
        f"{outcome_var}~T_season": (outcome_var, "T_season"),
        f"{outcome_var}~SW_season": (outcome_var, "SW_season"),
    }.items():
        path = _build_path(coef_map, lhs, rhs)
        if path:
            direct_paths[key] = path

    fig = plt.figure(figsize=(11, 12))
    fig.patch.set_facecolor("white")
    ax1 = plt.subplot2grid((4, 1), (0, 0), rowspan=3)
    ax2 = plt.subplot2grid((4, 1), (3, 0))
    ax1.set_facecolor("white")
    ax2.set_facecolor("#F8F8F8")
    ax1.set_xlim(0, 11)
    ax1.set_ylim(0, 12)

    positions = {
        "P_PRE": (1.5, 9.5),
        "T_PRE": (5.5, 9.5),
        "SW_PRE": (9.5, 9.5),
        "EOS": (1.5, 5.5),
        "NDVI": (5.5, 5.5),
        "TRATE": (9.5, 5.5),
        "P_SEASON": (1.5, 1.5),
        "T_SEASON": (5.5, 1.5),
        "SW_SEASON": (9.5, 1.5),
    }

    var_to_pos = {
        "P_pre": "P_PRE",
        "T_pre": "T_PRE",
        "SW_pre": "SW_PRE",
        "EOS": "EOS",
        middle_var: "NDVI",
        outcome_var: "TRATE",
        "P_season": "P_SEASON",
        "T_season": "T_SEASON",
        "SW_season": "SW_SEASON",
    }

    def draw_box(ax, name, pos, is_outcome=False):
        x, y = pos
        if is_outcome:
            width, height = 2.2, 1.4
            fontsize, r2_fontsize = 24, 18
        else:
            width, height = 2.0, 1.2
            fontsize, r2_fontsize = 20, 16

        shadow = patches.Rectangle(
            (x - width / 2 + 0.03, y - height / 2 - 0.03),
            width, height,
            facecolor="lightgray", alpha=0.3, zorder=1
        )
        ax.add_patch(shadow)

        box = patches.Rectangle(
            (x - width / 2, y - height / 2),
            width, height,
            facecolor="#E8F5E9" if is_outcome else "#F5F5F5",
            edgecolor="#4CAF50" if is_outcome else "#757575",
            linewidth=2.5 if is_outcome else 2, zorder=2
        )
        ax.add_patch(box)

        label = display_labels.get(name, name)
        if name in r2_values and np.isfinite(r2_values[name]):
            ax.text(x, y + 0.15, label, ha="center", va="center",
                    fontsize=fontsize, fontweight="bold", color="black")
            ax.text(x, y - 0.25, f"R² = {r2_values[name]:.3f}",
                    ha="center", va="center",
                    fontsize=r2_fontsize, color="black", fontweight="bold")
        else:
            ax.text(x, y, label, ha="center", va="center",
                    fontsize=fontsize, fontweight="bold", color="black")

    def draw_arrow(ax, start_pos, end_pos, coef, significant, start_name, end_name, is_eos_to_out=False):
        x1, y1 = start_pos
        x2, y2 = end_pos

        if start_name in ["EOS", "NDVI", "TRATE"]:
            start_width, start_height = 2.2, 1.4
        else:
            start_width, start_height = 2.0, 1.2

        if end_name in ["EOS", "NDVI", "TRATE"]:
            end_width, end_height = 2.2, 1.4
        else:
            end_width, end_height = 2.0, 1.2

        dx = x2 - x1
        dy = y2 - y1

        if start_name in ["P_PRE", "T_PRE", "SW_PRE"]:
            start_x, start_y = x1, y1 - start_height / 2
        elif start_name in ["P_SEASON", "T_SEASON", "SW_SEASON"]:
            start_x, start_y = x1, y1 + start_height / 2
        else:
            if end_name in ["P_SEASON", "T_SEASON", "SW_SEASON"]:
                start_x, start_y = x1, y1 - start_height / 2
            elif end_name in ["P_PRE", "T_PRE", "SW_PRE"]:
                start_x, start_y = x1, y1 + start_height / 2
            else:
                if dx > 0:
                    start_x, start_y = x1 + start_width / 2, y1
                else:
                    start_x, start_y = x1 - start_width / 2, y1

        if end_name in ["P_PRE", "T_PRE", "SW_PRE"]:
            end_x, end_y = x2, y2 - end_height / 2
        elif end_name in ["P_SEASON", "T_SEASON", "SW_SEASON"]:
            end_x, end_y = x2, y2 + end_height / 2
        else:
            if start_name in ["P_PRE", "T_PRE", "SW_PRE"]:
                end_x, end_y = x2, y2 + end_height / 2
            elif start_name in ["P_SEASON", "T_SEASON", "SW_SEASON"]:
                end_x, end_y = x2, y2 - end_height / 2
            else:
                if dx > 0:
                    end_x, end_y = x2 - end_width / 2, y2
                else:
                    end_x, end_y = x2 + end_width / 2, y2

        if significant:
            color = "#D32F2F" if coef > 0 else "#1976D2"
            linestyle = "-"
            abs_coef = abs(coef)
            if abs_coef >= 0.4:
                linewidth = 6
            elif abs_coef >= 0.3:
                linewidth = 5
            elif abs_coef >= 0.2:
                linewidth = 4
            elif abs_coef >= 0.1:
                linewidth = 2.5
            else:
                linewidth = 1.5
            alpha = 1.0
        else:
            color = "#9E9E9E"
            linestyle = "--"
            linewidth = 2
            alpha = 0.6

        is_pre_to_mid = (end_name == "NDVI" and start_name in ["P_PRE", "T_PRE", "SW_PRE"])
        is_season_to_mid = (end_name == "NDVI" and start_name in ["P_SEASON", "T_SEASON", "SW_SEASON"])
        is_eos_to_mid = (start_name == "EOS" and end_name == "NDVI")
        is_mid_to_out = (start_name == "NDVI" and end_name == "TRATE")

        if is_eos_to_out:
            arrow = FancyArrowPatch(
                (start_x, start_y), (end_x, end_y),
                connectionstyle="arc3,rad=0.3",
                arrowstyle="->",
                mutation_scale=25,
                color=color,
                linewidth=linewidth,
                linestyle=linestyle,
                alpha=alpha,
                zorder=2
            )
            ax.add_patch(arrow)
            t = 0.2
            straight_x = start_x + (end_x - start_x) * t
            straight_y = start_y + (end_y - start_y) * t
            chord_length = ((end_x - start_x) ** 2 + (end_y - start_y) ** 2) ** 0.5
            arc_offset_y = -0.3 * np.sin(np.pi * t) * chord_length * 0.6
            label_x = straight_x
            label_y = straight_y + arc_offset_y
            dx_dt = (end_x - start_x)
            dy_dt = (end_y - start_y) - 0.3 * np.pi * np.cos(np.pi * t) * chord_length * 0.6
            angle = np.arctan2(dy_dt, dx_dt) * 180 / np.pi
            if angle > 90:
                angle -= 180
            if angle < -90:
                angle += 180
        else:
            ax.annotate("", xy=(end_x, end_y), xytext=(start_x, start_y),
                        arrowprops=dict(arrowstyle="->", color=color,
                                        linewidth=linewidth, mutation_scale=25,
                                        linestyle=linestyle, alpha=alpha, zorder=2))

            if is_pre_to_mid or is_season_to_mid:
                label_x = start_x + (end_x - start_x) * 0.33
                label_y = start_y + (end_y - start_y) * 0.33
            elif is_eos_to_mid or is_mid_to_out:
                label_x = start_x + (end_x - start_x) * 0.5
                label_y = start_y + (end_y - start_y) * 0.5
            else:
                label_x = start_x + (end_x - start_x) * 0.67
                label_y = start_y + (end_y - start_y) * 0.67

            angle = np.arctan2(end_y - start_y, end_x - start_x) * 180 / np.pi
            if angle > 90:
                angle -= 180
            if angle < -90:
                angle += 180

        ax.text(label_x, label_y, f"{coef:.2f}" if abs(coef) >= 0.1 else f"{coef:.3f}",
                ha="center", va="center", fontsize=18, fontweight="bold", color="black",
                rotation=angle, rotation_mode="anchor",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor="none", alpha=0.95),
                zorder=100)

    for name, pos in positions.items():
        is_outcome = name in ["EOS", "NDVI", "TRATE"]
        draw_box(ax1, name, pos, is_outcome)

    for path_key, path_data in direct_paths.items():
        parts = path_key.split("~")
        if len(parts) == 2:
            end_var, start_var = parts
            if start_var in var_to_pos and end_var in var_to_pos:
                start_name = var_to_pos[start_var]
                end_name = var_to_pos[end_var]
                is_eos_out = (start_name == "EOS" and end_name == "TRATE")
                draw_arrow(ax1, positions[start_name], positions[end_name],
                           path_data["coef"], path_data["sig"],
                           start_name, end_name, is_eos_out)

    ax1.text(5.5, 11.2, title_text,
             fontsize=22, fontweight="bold", ha="center",
             bbox=dict(boxstyle="round,pad=0.8", facecolor="white",
                       edgecolor="#666666", linewidth=3))

    legend_elements = [
        plt.Line2D([0], [0], color="#D32F2F", linewidth=6, label="Positive effect"),
        plt.Line2D([0], [0], color="#1976D2", linewidth=6, label="Negative effect"),
        plt.Line2D([0], [0], color="#9E9E9E", linewidth=4, linestyle="--", label="Non-significant"),
    ]
    ax1.legend(handles=legend_elements, loc="upper left", fontsize=14, frameon=True,
               fancybox=True, shadow=True, bbox_to_anchor=(0.01, 0.99))

    ax1.set_xticks([])
    ax1.set_yticks([])
    for spine in ax1.spines.values():
        spine.set_visible(False)

    def get_coef(lhs, rhs):
        item = coef_map.get((lhs, rhs))
        return item["coef"] if item else np.nan

    e1 = get_coef(outcome_var, "P_season")
    e2 = get_coef(outcome_var, "T_season")
    e3 = get_coef(outcome_var, "SW_season")
    c1 = get_coef(middle_var, "P_season")
    c2 = get_coef(middle_var, "T_season")
    c3 = get_coef(middle_var, "SW_season")
    d_coef = get_coef(outcome_var, middle_var)

    direct_effects = [e1, e2, e3]
    indirect_via_mid = [c1 * d_coef, c2 * d_coef, c3 * d_coef]
    total_effects = [d + i for d, i in zip(direct_effects, indirect_via_mid)]
    mediation_ratios = [
        (i / t) if np.isfinite(t) and t != 0 else np.nan
        for i, t in zip(indirect_via_mid, total_effects)
    ]

    x_pos = np.arange(3)
    width = 0.5
    ax2.bar(x_pos, direct_effects, width,
            label=bar_direct_label,
            color="#9C27B0", alpha=0.85, edgecolor="#7B1FA2", linewidth=2)
    ax2.bar(x_pos, indirect_via_mid, width, bottom=direct_effects,
            label=bar_indirect_label,
            color="#FF9800", alpha=0.85, edgecolor="#F57C00", linewidth=2)

    for i, (total, mediation_pct) in enumerate(zip(total_effects, mediation_ratios)):
        pct_text = f"{mediation_pct*100:.1f}%" if np.isfinite(mediation_pct) else "NA"
        ax2.text(x_pos[i], total + 0.02, f"{total:.3f}\n({pct_text})",
                 ha="center", va="bottom", fontsize=16, fontweight="bold")

    ymin = min(0, np.nanmin(total_effects + direct_effects + indirect_via_mid)) - 0.05
    ymax = max(0, np.nanmax(total_effects + direct_effects + indirect_via_mid)) + 0.05
    ax2.set_xlim(-0.5, 2.5)
    ax2.set_ylim(ymin, ymax)
    ax2.axhline(y=0, color="black", linewidth=2.5)
    ax2.set_ylabel(bar_ylabel, fontsize=18, fontweight="bold")
    ax2.set_xlabel("Season Climate Variables (Fixed Window)", fontsize=18, fontweight="bold")
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(["Precipitation", "Temperature", "Solar Radiation"],
                        fontsize=16, fontweight="bold")
    ax2.tick_params(axis="both", which="major", labelsize=16)
    legend = ax2.legend(loc="upper left", fontsize=16, frameon=True, fancybox=True, shadow=True)
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_alpha(0.9)
    ax2.grid(True, alpha=0.3, axis="y", linestyle="-", linewidth=0.5)
    ax2.set_axisbelow(True)
    for spine in ["top", "right"]:
        ax2.spines[spine].set_visible(False)
    ax2.spines["left"].set_linewidth(2.5)
    ax2.spines["bottom"].set_linewidth(2.5)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()

    print(f"已保存: {out_png} (mode={mode})")


# ==================== 04c 绘图任务 ====================
@dataclass
class PlotTask:
    r_file: Path
    p_file: Path | None
    out_file: Path
    var_label: str
    cbar_title: str | None = None
    value_range: tuple | None = None
    auto_range: str | None = None


def _label_to_pair(label):
    if "_vs_" in label:
        left, right = label.split("_vs_", 1)
        return f"{left} ~ {right}"
    return label


def _collect_04c_tasks(base_dir, out_dir):
    tasks = []
    sec32 = base_dir / "Section_3.2_Phenology_Impact"
    if sec32.exists():
        for slope_file in sorted(sec32.glob("*_slope.tif")):
            p_file = slope_file.with_name(slope_file.name.replace("_slope.tif", "_pvalue.tif"))
            r2_file = slope_file.with_name(slope_file.name.replace("_slope.tif", "_R2.tif"))
            label = slope_file.stem.replace("_slope", "")
            cbar_title = f"Slope, {_label_to_pair(label)}"
            tasks.append(PlotTask(
                r_file=slope_file,
                p_file=p_file if p_file.exists() else None,
                out_file=out_dir / "Section_3.2_Phenology_Impact" / f"{slope_file.stem}.tif",
                var_label=label,
                cbar_title=cbar_title,
                auto_range="symmetric",
            ))
            if r2_file.exists():
                tasks.append(PlotTask(
                    r_file=r2_file,
                    p_file=None,
                    out_file=out_dir / "Section_3.2_Phenology_Impact" / f"{r2_file.stem}.tif",
                    var_label=f"R2_{label}",
                    cbar_title=f"R², {_label_to_pair(label)}",
                    value_range=(0.0, 1.0),
                ))

    full_period = base_dir / "Section_3.3_Drivers" / "Full_Period"
    if full_period.exists():
        for resp_dir in sorted(full_period.iterdir()):
            if not resp_dir.is_dir():
                continue
            for r_file in sorted(resp_dir.glob("partial_r_*.tif")):
                var = r_file.stem.replace("partial_r_", "")
                p_file = resp_dir / f"partial_p_{var}.tif"
                label = f"{var}_vs_{resp_dir.name}"
                tasks.append(PlotTask(
                    r_file=r_file,
                    p_file=p_file if p_file.exists() else None,
                    out_file=out_dir / "Section_3.3_Drivers" / "Full_Period" / resp_dir.name / f"{r_file.stem}.tif",
                    var_label=label,
                ))

            r2_file = resp_dir / "R_squared.tif"
            if r2_file.exists():
                tasks.append(PlotTask(
                    r_file=r2_file,
                    p_file=None,
                    out_file=out_dir / "Section_3.3_Drivers" / "Full_Period" / resp_dir.name / "R_squared.tif",
                    var_label=f"R2_vs_{resp_dir.name}",
                    cbar_title=f"R², {resp_dir.name}",
                    value_range=(0.0, 1.0),
                ))

            for v_file in sorted(resp_dir.glob("vif_retained_*.tif")):
                var = v_file.stem.replace("vif_retained_", "")
                tasks.append(PlotTask(
                    r_file=v_file,
                    p_file=None,
                    out_file=out_dir / "Section_3.3_Drivers" / "Full_Period" / resp_dir.name / f"{v_file.stem}.tif",
                    var_label=f"VIF_{var}_vs_{resp_dir.name}",
                    cbar_title=f"VIF retained, {var}",
                    value_range=(0.0, 1.0),
                ))

    moving_dir = base_dir / "Section_3.3_Drivers" / "Moving_Window" / "TR_fixed_window"
    if moving_dir.exists():
        for r_file in sorted(moving_dir.glob("partial_r_*.tif")):
            label = r_file.stem.replace("partial_r_", "")
            tasks.append(PlotTask(
                r_file=r_file,
                p_file=None,
                out_file=out_dir / "Section_3.3_Drivers" / "Moving_Window" / "TR_fixed_window" / f"{r_file.stem}.tif",
                var_label=label,
            ))

    trend_dir = base_dir / "Section_3.3_Drivers" / "Sensitivity_Trends" / "TR_fixed_window"
    if trend_dir.exists():
        for slope_file in sorted(trend_dir.glob("*_trend_slope.tif")):
            base = slope_file.stem.replace("_trend_slope", "")
            p_file = trend_dir / f"{base}_trend_pvalue.tif"
            cbar_title = f"Trend slope, {base}"
            tasks.append(PlotTask(
                r_file=slope_file,
                p_file=p_file if p_file.exists() else None,
                out_file=out_dir / "Section_3.3_Drivers" / "Sensitivity_Trends" / "TR_fixed_window" / f"{slope_file.stem}.tif",
                var_label=base,
                cbar_title=cbar_title,
                auto_range="symmetric",
            ))
    return tasks


def plot_04c_all():
    base_root = ANALYSIS_DIR / "Statistical_Analysis_FixedWindow"
    for run_label in ["Raw", "Detrended"]:
        base_dir = base_root / run_label
        if not base_dir.exists():
            print(f"04c 输出目录不存在，跳过: {base_dir}")
            continue
        out_dir = OUTPUT_FIG_DIR / "04c" / run_label
        tasks = _collect_04c_tasks(base_dir, out_dir)
        print(f"04c({run_label}) 绘图任务数: {len(tasks)}")
        for task in tasks:
            plot_corr_style_map(
                r_file=task.r_file,
                p_file=task.p_file,
                out_figure=task.out_file,
                var_label=task.var_label,
                cbar_title=task.cbar_title,
                value_range=task.value_range,
                auto_range=task.auto_range,
            )


# ==================== 05 SEM 绘图任务 ====================
def _resolve_sem_dir(dir_name, fallback_file):
    primary = ANALYSIS_DIR / dir_name
    if primary.exists():
        return primary
    for fb in SEM_FALLBACK_DIRS:
        fb_dir = fb / dir_name
        if fb_dir.exists():
            return fb_dir
        if (fb / fallback_file).exists():
            return fb
    return primary


def plot_sem_all(only_group=None):
    sem_out_root = OUTPUT_FIG_DIR / "SEM"
    sem_out_global = sem_out_root / "Global"
    sem_out_pixelwise = sem_out_root / "Pixelwise"
    sem_out_global.mkdir(parents=True, exist_ok=True)
    sem_out_pixelwise.mkdir(parents=True, exist_ok=True)

    run_variants = [
        {"label": "Raw", "suffix": "", "title_suffix": "Raw"},
        {"label": "Detrended", "suffix": "_detrended", "title_suffix": "Detrended"},
    ]
    for run in run_variants:
        (sem_out_global / run["label"]).mkdir(parents=True, exist_ok=True)
        (sem_out_pixelwise / run["label"]).mkdir(parents=True, exist_ok=True)

    sem_tasks = []

    sem_tasks.append({
        "name": f"05b_EOS_{MIDDLE_VAR_NAME}_Trate",
        "group": "05b",
        "dir_name": "SEM_Results_Dual_Fixed",
        "param": "SEM_dual_timescale_parameters.csv",
        "r2": "SEM_dual_timescale_R2.csv",
        "middle_var": "Fixed_GPPrate",
        "outcome_var": "Fixed_Trate",
        "middle_label": MIDDLE_VAR_NAME,
        "outcome_label": r"T$_{rate}$",
        "title": "EOS - Dual-Timescale SEM (Fixed Window, Full Paths)",
        "bar_direct": "Direct (Season Climate → TRate)",
        "bar_indirect": f"Indirect (Season Climate → {MIDDLE_VAR_NAME} → TRate)",
        "bar_ylabel": "Standardized Effect on TRate",
    })

    sem_tasks.append({
        "name": "05c_Pooled_EOS",
        "group": "05c",
        "dir_name": "SEM_Results_Dual_Fixed_Robust_Pooled_EOS",
        "param": "SEM_dual_timescale_parameters.csv",
        "r2": "SEM_dual_timescale_R2.csv",
        "middle_var": "Fixed_GPPrate",
        "outcome_var": "Fixed_Trate",
        "middle_label": MIDDLE_VAR_NAME,
        "outcome_label": r"T$_{rate}$",
        "title": "EOS - Pooled SEM (Pixel-Year, Full Paths)",
        "bar_direct": "Direct (Season Climate → TRate)",
        "bar_indirect": f"Indirect (Season Climate → {MIDDLE_VAR_NAME} → TRate)",
        "bar_ylabel": "Standardized Effect on TRate",
    })

    sem_tasks.append({
        "name": "05d_Lavaan_Ours",
        "group": "05d",
        "dir_name": "SEM_Results_Dual_Fixed_Lavaan_Compare_Ours",
        "param": "SEM_dual_timescale_parameters.csv",
        "r2": "SEM_dual_timescale_R2.csv",
        "middle_var": "Fixed_GPPrate",
        "outcome_var": "Fixed_Trate",
        "middle_label": MIDDLE_VAR_NAME,
        "outcome_label": r"T$_{rate}$",
        "title": "EOS - Lavaan Pixelwise SEM (Ours)",
        "bar_direct": "Direct (Season Climate → TRate)",
        "bar_indirect": f"Indirect (Season Climate → {MIDDLE_VAR_NAME} → TRate)",
        "bar_ylabel": "Standardized Effect on TRate",
        "pixel_only": True,
    })

    sem_tasks.append({
        "name": "05d_Lavaan_Other",
        "group": "05d",
        "dir_name": "SEM_Results_Dual_Fixed_Lavaan_Compare_Other",
        "param": "SEM_dual_timescale_parameters.csv",
        "r2": "SEM_dual_timescale_R2.csv",
        "middle_var": "Fixed_GPPrate",
        "outcome_var": "Fixed_Trate",
        "middle_label": MIDDLE_VAR_NAME,
        "outcome_label": r"T$_{rate}$",
        "title": "EOS - Lavaan Pixelwise SEM (Other)",
        "bar_direct": "Direct (Season Climate → TRate)",
        "bar_indirect": f"Indirect (Season Climate → {MIDDLE_VAR_NAME} → TRate)",
        "bar_ylabel": "Standardized Effect on TRate",
        "pixel_only": True,
    })

    for task in sem_tasks:
        if only_group and task.get("group") != only_group:
            continue
        for run in run_variants:
            run_dir = _resolve_sem_dir(f"{task['dir_name']}{run['suffix']}", task["param"])
            if not task.get("pixel_only", False):
                param_file = run_dir / task["param"]
                r2_file = run_dir / task["r2"]
                if not (param_file.exists() and r2_file.exists()):
                    print(f"05 SEM 输出缺失，跳过: {run_dir}")
                    continue
                out_png = (sem_out_global / run["label"]) / f"{task['name']}.png"
                plot_sem_full_path(
                    param_file=param_file,
                    r2_file=r2_file,
                    out_png=out_png,
                    title_text=f"{task['title']} ({run['title_suffix']})",
                    middle_var=task["middle_var"],
                    outcome_var=task["outcome_var"],
                    middle_label=task["middle_label"],
                    outcome_label=task["outcome_label"],
                    bar_direct_label=task["bar_direct"],
                    bar_indirect_label=task["bar_indirect"],
                    bar_ylabel=task["bar_ylabel"],
                )

            pix_dir = run_dir / "Pixelwise"
            if pix_dir.exists():
                pix_param = pix_dir / task["param"]
                pix_r2 = pix_dir / task["r2"]
                if pix_param.exists() and pix_r2.exists():
                    out_png = (sem_out_pixelwise / run["label"]) / f"{task['name']}_pixelwise.png"
                    plot_sem_full_path(
                        param_file=pix_param,
                        r2_file=pix_r2,
                        out_png=out_png,
                        title_text=f"{task['title']} ({run['title_suffix']} · Pixelwise summary)",
                        middle_var=task["middle_var"],
                        outcome_var=task["outcome_var"],
                        middle_label=task["middle_label"],
                        outcome_label=task["outcome_label"],
                        bar_direct_label=task["bar_direct"],
                        bar_indirect_label=task["bar_indirect"],
                        bar_ylabel=task["bar_ylabel"],
                    )
                elif task.get("pixel_only", False):
                    print(f"05 SEM 像元级输出缺失，跳过: {pix_dir}")


def _parse_plot_args():
    parser = argparse.ArgumentParser(description="Wang2025 plotting entrypoint")
    parser.add_argument(
        "--only",
        choices=["all", "04c", "05b", "05c", "05d"],
        default=PLOT_ONLY_DEFAULT,
        help="only plot a specific module (default: all)"
    )
    return parser.parse_args()


def main():
    args = _parse_plot_args()
    only = args.only
    print("=" * 70)
    print("绘图入口：04c + 05b/05c/05d")
    print("=" * 70)
    if only in ("all", "04c"):
        plot_04c_all()
    if only in ("all", "05b", "05c", "05d"):
        plot_sem_all(only_group=None if only == "all" else only)
    print("\n✓ 绘图完成")
    print(f"输出目录: {OUTPUT_FIG_DIR}")


if __name__ == "__main__":
    main()
