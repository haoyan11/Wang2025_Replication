# -*- coding: utf-8 -*-
"""
Fast north-polar maps of Sen trend (1982-2022) with significance overlay
=======================================================================
· 绘制蒸腾量 T 和非蒸腾量 E_veg（canopy + bare soil）
  —— 年尺度 & 春 / 夏 / 秋 季节尺度趋势图（Sen + MK 显著性叠加）
· 额外绘制 SOS / POS / EOS 的趋势图（Sen + MK 显著性叠加）
· 所有图片保存在同一个输出文件夹 OUT_FIG_DIR 下
"""

import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.path as mpath
import matplotlib.patches as patches
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from osgeo import gdal
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap
import colorsys

# ------------------------------------------------------------------
# A. 全局样式
# ------------------------------------------------------------------
FONT_SIZE = 9  # 字号

mpl.rcParams.update({
    'font.family': 'Arial',
    'font.size': FONT_SIZE,
    'axes.titlesize': FONT_SIZE,
    'axes.labelsize': FONT_SIZE,
    'xtick.labelsize': FONT_SIZE,
    'ytick.labelsize': FONT_SIZE,
    'axes.linewidth': .3,
    'figure.dpi': 300,
    'savefig.bbox': 'tight',
    'mathtext.fontset': 'stixsans',  # 使用无衬线数学字体，与 Arial 更匹配
    'mathtext.default': 'regular',
})
plt.rcParams['agg.path.chunksize'] = 10000

# ------------------------------------------------------------------
# 基础路径 & 输出图片统一文件夹
# ------------------------------------------------------------------
ERA5_BASE = r"I:\F\Data4\Meteorological Data\ERA5_Land"

# 物候 Sen & MK 结果目录
PHENO_SEN_DIR = r"I:\F\Data4\GIMMS_NDVI_phenology_correct\GIMMS_NDVI_phenology_1_northern\sen_trend"

OUT_FIG_DIR = r"I:\F\Data4\Sen_trend_maps_T_E_Pheno"
os.makedirs(OUT_FIG_DIR, exist_ok=True)

# ---------------- 各变量单位 ----------------
# var_root = var_label.split('_')[0]
# T/E 来自 ET（mm），物候是天
VAR_UNITS = {
    "ET": "mm",       # T / E → mm
    "sos": "days",    # SOS/POS/EOS → days/y
    "pos": "days",
    "eos": "days",
}

# ---------------- 显著性配置 ----------------
significance_levels = {
    'p05': {
        'edgecolor': 'none',
        'facecolor': 'black',
        'linewidth': 0,
        'alpha': 0.4
    }
}

USE_NEIGHBORHOOD_FILTER = True
NEIGHBORHOOD_THRESHOLD = 0.5

NODATA_IN = -9999.0
NODATA_OUT = -9999.0


def load_zvalue_data(zvalue_file, factor=1):
    """读取 MK 检验 z 值文件并返回处理后的数组"""
    if zvalue_file is None or (not os.path.exists(zvalue_file)):
        return None

    try:
        ds = gdal.Open(zvalue_file, gdal.GA_ReadOnly)
        if ds is None:
            print(f'警告: 无法打开 z 值文件 {zvalue_file}')
            return None

        zvalue_arr = ds.ReadAsArray().astype(float)
        nodata = ds.GetRasterBand(1).GetNoDataValue()
        if nodata is not None:
            zvalue_arr[zvalue_arr == nodata] = np.nan

        if factor > 1:
            zvalue_arr = zvalue_arr[::factor, ::factor]

        return zvalue_arr

    except Exception as e:
        print(f'警告: 处理 z 值文件时出错 {zvalue_file}: {e}')
        return None


def apply_neighborhood_filter(sig_mask, threshold_ratio=0.5):
    """空间邻域过滤"""
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
        mode='constant',
        cval=0
    )

    valid_mask = np.ones_like(sig_mask, dtype=int)
    valid_count = ndimage.generic_filter(
        valid_mask,
        np.sum,
        footprint=structure,
        mode='constant',
        cval=0
    )

    with np.errstate(divide='ignore', invalid='ignore'):
        sig_ratio = sig_count / valid_count
        sig_ratio[valid_count == 0] = 0

    filtered_mask = sig_mask & (sig_ratio >= threshold_ratio)

    return filtered_mask


def add_significance_overlay(ax, zvalue_arr, extent, proj_pc, significance_levels,
                             use_neighborhood_filter=True):
    """在地图上添加 MK 显著性方块"""
    if zvalue_arr is None:
        return

    nrows, ncols = zvalue_arr.shape

    lon_min, lon_max, lat_min, lat_max = extent
    lon_res = (lon_max - lon_min) / ncols
    lat_res = (lat_max - lat_min) / nrows

    print(f"数据分辨率: 经度 {lon_res:.4f}°/像元, 纬度 {lat_res:.4f}°/像元")
    print("显著性方块大小与像元一致。")

    lon_centers = lon_min + (np.arange(ncols) + 0.5) * lon_res
    lat_centers = lat_max - (np.arange(nrows) + 0.5) * lat_res
    xx, yy = np.meshgrid(lon_centers, lat_centers)

    z_threshold = 1.96
    sig_mask = (np.abs(zvalue_arr) > z_threshold) & (~np.isnan(zvalue_arr))

    if use_neighborhood_filter and np.any(sig_mask):
        try:
            sig_mask = apply_neighborhood_filter(sig_mask, threshold_ratio=NEIGHBORHOOD_THRESHOLD)
            print(f"邻域过滤完成：只保留邻域内显著性比例≥{NEIGHBORHOOD_THRESHOLD}的聚集区域")
        except Exception as e:
            print(f"警告：邻域过滤失败 ({e})，使用原始显著性标记")

    if np.any(sig_mask):
        sig_x = xx[sig_mask]
        sig_y = yy[sig_mask]

        cfg = significance_levels['p05']
        marker_w = lon_res
        marker_h = lat_res

        for i in range(len(sig_x)):
            x_center = sig_x[i]
            y_center = sig_y[i]

            rect = patches.Rectangle(
                (x_center - marker_w / 2, y_center - marker_h / 2),
                marker_w,
                marker_h,
                linewidth=cfg['linewidth'],
                edgecolor=cfg['edgecolor'],
                facecolor=cfg['facecolor'],
                alpha=cfg['alpha'],
                transform=proj_pc,
                zorder=10
            )
            ax.add_patch(rect)

        print(f"已添加 {len(sig_x)} 个显著性标记")


def compute_axis_scaling(ticks):
    """决定是否做 10^n 缩放"""
    ticks = np.asarray(ticks)
    if ticks.size == 0:
        return 1.0, 0

    max_abs = np.nanmax(np.abs(ticks))
    if max_abs == 0 or np.isnan(max_abs):
        return 1.0, 0

    exp = int(np.floor(np.log10(max_abs)))
    if abs(exp) <= 1:
        return 1.0, 0
    else:
        scale_factor = 10.0 ** exp
        return scale_factor, exp


def format_tick_label(value: float, decimals: int = 1) -> str:
    """把 0.0 变成 0，其他去掉多余 0"""
    if np.isnan(value):
        return ""
    threshold = 0.5 * (10 ** (-decimals))
    if abs(value) < threshold:
        return "0"
    fmt = f"{{:.{decimals}f}}"
    s = fmt.format(value)
    s = s.rstrip('0').rstrip('.')
    return s


def add_polar_lon_labels(ax,
                         lon_values,
                         lat_label=35,
                         proj_src=ccrs.PlateCarree(),
                         font_size=FONT_SIZE):
    """
    外圈经度标注，沿圆环切线方向
    """
    proj = ax.projection
    data_to_axes = ax.transAxes.inverted().transform

    cx, cy = 0.5, 0.5
    r_ring = 0.5
    margin = 0.02
    r_text = r_ring + margin

    for lon in lon_values:
        # 下方不标，避免和色条/直方图冲突
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

        # 强制把标注锚点缩放到统一半径 r_text，上下左右到圆环的“距离”一致
        scale = r_text / r
        x_text = cx + vx * scale
        y_text = cy + vy * scale

        # 去掉最顶部的 180°
        if label == "180°" and y_text > 0.95:
            continue

        theta = np.arctan2(vy, vx)
        angle_deg = np.degrees(theta) + 90
        if angle_deg > 90:
            angle_deg -= 180
        if angle_deg < -90:
            angle_deg += 180

        ax.text(
            x_text, y_text, label,
            transform=ax.transAxes,
            ha='center', va='center',
            rotation=angle_deg,
            rotation_mode='anchor',
            fontsize=font_size
        )


def add_lat_labels(ax,
                   lat_values=(45, 60, 75),
                   lon_for_labels=-130,
                   proj_src=ccrs.PlateCarree(),
                   font_size=FONT_SIZE):
    """
    在同一条经线（半径方向）上添加 45°N, 60°N, 75°N 标注
    · lon_for_labels = -130°，整体略向右
    """
    for lat in lat_values:
        label = f"{int(lat)}°N"
        ax.text(
            lon_for_labels, lat, label,
            transform=proj_src,
            ha='left', va='center',
            fontsize=font_size
        )


def plot_trend_map(trend_file, mk_file, out_figure, var_label):
    """
    核心绘图函数
    """
    if not os.path.exists(trend_file):
        print(f"趋势文件不存在，跳过: {trend_file}")
        return

    print(f"\n绘图: {var_label}")
    print(f"  trend: {trend_file}")
    print(f"  mk   : {mk_file if mk_file and os.path.exists(mk_file) else 'None or missing'}")

    ds = gdal.Open(trend_file, gdal.GA_ReadOnly)
    if ds is None:
        print(f"GDAL 无法打开文件 {trend_file}，跳过。")
        return

    arr = ds.ReadAsArray().astype(float)
    gt = ds.GetGeoTransform()
    nodata = ds.GetRasterBand(1).GetNoDataValue()
    if nodata is not None:
        arr[arr == nodata] = np.nan

    nrows, ncols = arr.shape
    extent = [gt[0], gt[0] + gt[1] * ncols,
              gt[3] + gt[5] * nrows, gt[3]]

    # 下采样（大图时提速）
    max_pixels = 10_000_000
    factor = int(np.ceil(np.sqrt((nrows * ncols) / max_pixels)))
    if factor > 1:
        arr = arr[::factor, ::factor]

    zvalue_arr = load_zvalue_data(mk_file, factor=factor)

    valid_data = arr[~np.isnan(arr)]
    if valid_data.size == 0:
        print("全部为 NaN，跳过此图。")
        return

    vmax = np.nanpercentile(np.abs(valid_data), 99)
    if vmax <= 0:
        vmax = np.nanmax(np.abs(valid_data))
    if vmax == 0 or np.isnan(vmax):
        vmax = 1.0

    boundaries = np.linspace(-vmax, vmax, 9)

    colors = [
        (95 / 255, 175 / 255, 250 / 255),
        (175 / 255, 215 / 255, 255 / 255),
        (1, 1, 1),
        (1, 220 / 255, 185 / 255),
        (1, 160 / 255, 115 / 255)
    ]

    brighten = 0.9
    for i, (r, g, b) in enumerate(colors):
        h, l, s = colorsys.rgb_to_hls(r, g, b)
        l = max(0, min(1, l * brighten))
        colors[i] = colorsys.hls_to_rgb(h, l, s)

    cmap = LinearSegmentedColormap.from_list('BrightDiv', colors, N=256)
    norm = BoundaryNorm(boundaries, cmap.N)

    proj_map = ccrs.NorthPolarStereo()
    proj_pc = ccrs.PlateCarree()

    fig = plt.figure(figsize=(3.2, 3.8))
    ax = plt.axes(projection=proj_map)
    ax.set_extent([-180, 180, 30, 90], crs=proj_pc)

    ax.add_feature(cfeature.LAND.with_scale('110m'),
                   facecolor='none', edgecolor='k',
                   linewidth=.3, zorder=1)
    ax.add_feature(cfeature.COASTLINE.with_scale('110m'),
                   linewidth=.3, zorder=1)

    theta = np.linspace(0, 2 * np.pi, 360)
    verts = np.column_stack([0.5 + 0.5 * np.sin(theta),
                             0.5 + 0.5 * np.cos(theta)])
    ax.set_boundary(mpath.Path(verts, closed=True), transform=ax.transAxes)

    mesh = ax.imshow(arr, origin='upper', extent=extent,
                     transform=proj_pc, cmap=cmap, norm=norm,
                     interpolation='nearest', zorder=2)

    add_significance_overlay(ax, zvalue_arr, extent, proj_pc,
                             significance_levels,
                             use_neighborhood_filter=USE_NEIGHBORHOOD_FILTER)

    # 网格线：只画线，不画文字
    gl = ax.gridlines(crs=proj_pc,
                      draw_labels=False,
                      linewidth=.4, color='grey',
                      linestyle=(0, (5, 2)), alpha=.8)
    gl.xlocator = mticker.FixedLocator(np.arange(-180, 181, 30))
    gl.ylocator = mticker.FixedLocator([45, 60, 75])

    # 经度：沿圆环切线方向的一圈数字
    lon_values = np.arange(-180, 181, 30)
    add_polar_lon_labels(ax, lon_values, lat_label=35,
                         proj_src=proj_pc, font_size=FONT_SIZE)

    # 纬度：沿同一条经线（半径方向）45°N, 60°N, 75°N
    add_lat_labels(ax,
                   lat_values=(45, 60, 75),
                   lon_for_labels=-130,
                   proj_src=proj_pc,
                   font_size=FONT_SIZE)

    # 颜色条 & 缩放
    ticks = boundaries[::2]
    scale_factor, exp = compute_axis_scaling(ticks)

    cax = inset_axes(ax, width="52%", height="3%",
                     loc='lower right',
                     bbox_to_anchor=(0.0, -0.08, 1, 1),
                     bbox_transform=ax.transAxes,
                     borderpad=0)

    cbar = plt.colorbar(mesh, cax=cax, orientation='horizontal',
                        ticks=ticks)

    # ---- 变量名和单位 ---- #
    parts = var_label.split('_')
    var_root = parts[0]          # 如 ET_transp_Annual → "ET"; sos → "sos"
    unit = VAR_UNITS.get(var_root, "")

    # 显示名映射：T/E/SOS/POS/EOS
    if var_root in ("sos", "pos", "eos"):
        display_code = var_root.upper()   # "SOS" / "POS" / "EOS"
    else:
        # ET 组件：ET_transp_* → T；ET_Eveg_* → E
        if var_label.startswith("ET_transp"):
            display_code = "T"
        elif var_label.startswith("ET_Eveg"):
            display_code = "E"
        else:
            display_code = var_label  # fallback（理论上不会用到）

    scaled_values = [t / scale_factor for t in ticks]
    scaled_ticklabels = [format_tick_label(v, decimals=1) for v in scaled_values]
    cbar.ax.set_xticklabels(scaled_ticklabels)

    cax.tick_params(length=0, width=0, pad=1)
    for spine in cax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(.3)

    if exp != 0:
        unit_text = rf"$\times 10^{{{exp}}}$ {unit}"
    else:
        unit_text = unit

    # 例：Trend, T (mm/y)；Trend, E (mm/y)；Trend, SOS (days/y)
    cbar_title = f"Trend, {display_code} ({unit_text}/y)"

    cbar.ax.text(0.5, 1.03, cbar_title,
                 transform=cbar.ax.transAxes,
                 ha='center', va='bottom', fontsize=FONT_SIZE)

    # 直方图
    hist_ax = inset_axes(ax, width="34%", height="34%",
                         loc='lower left',
                         bbox_to_anchor=(0.05, -0.06, 1, 1),
                         bbox_transform=ax.transAxes,
                         borderpad=0)

    bins = boundaries
    hist, _ = np.histogram(valid_data, bins=bins)
    percent = hist / hist.sum() * 100
    centers = (bins[:-1] + bins[1:]) / 2

    bar_width = centers[1] - centers[0]
    hist_ax.bar(centers / scale_factor, percent,
                width=(bar_width / scale_factor) * 0.9,
                color=cmap(norm(centers)), edgecolor='none', align='center')

    hist_ticks = ticks
    hist_tick_values = [t / scale_factor for t in hist_ticks]
    hist_ax.set_xticks(hist_tick_values)
    hist_ax.set_xticklabels(
        [format_tick_label(v, decimals=1) for v in hist_tick_values]
    )

    hist_ax.set_ylim(0, max(60, np.nanmax(percent) * 1.2))
    hist_ax.set_yticks([15, 30, 45, 60])
    hist_ax.tick_params(direction='out', length=2, width=.3, pad=1)
    for spine in hist_ax.spines.values():
        spine.set_linewidth(.3)

    hist_ax.set_ylabel('Frequency (%)', labelpad=3, fontsize=FONT_SIZE)
    hist_ax.axhline(0, color='black', linewidth=.3)

    # P / N 比例
    total_pixels = valid_data.size
    pos_trend = np.sum(valid_data >= 0)
    neg_trend = np.sum(valid_data < 0)
    pos_percent = pos_trend / total_pixels * 100
    neg_percent = neg_trend / total_pixels * 100

    sig_pos_percent = 0.0
    sig_neg_percent = 0.0

    if zvalue_arr is not None:
        valid_mask = ~np.isnan(zvalue_arr) & ~np.isnan(arr)
        valid_z = zvalue_arr[valid_mask]
        valid_trend = arr[valid_mask]
        if valid_z.size > 0:
            sig_mask = np.abs(valid_z) > 1.96
            sig_pos_pixels = np.sum(sig_mask & (valid_trend > 0))
            sig_neg_pixels = np.sum(sig_mask & (valid_trend < 0))
            sig_pos_percent = sig_pos_pixels / total_pixels * 100
            sig_neg_percent = sig_neg_pixels / total_pixels * 100

    pos_text = f'P: {pos_percent:.1f}% ({sig_pos_percent:.1f}%)'
    neg_text = f'N: {neg_percent:.1f}% ({sig_neg_percent:.1f}%)'

    ax.text(0.5, 0.93, pos_text, transform=ax.transAxes,
            ha='center', va='center', fontsize=FONT_SIZE, color='black')
    ax.text(0.5, 0.88, neg_text, transform=ax.transAxes,
            ha='center', va='center', fontsize=FONT_SIZE, color='black')

    plt.subplots_adjust(bottom=0.12, top=0.98, left=0.03, right=0.97)
    plt.savefig(out_figure, format='tiff', dpi=600)
    plt.close(fig)
    print(f"已保存: {out_figure}")


# ====================== 任务列表：T & E_veg + SOS/POS/EOS ======================= #

TASKS = [
    # ---------- 年尺度 T ---------- #
    dict(
        name="T_Annual",
        trend_dir=os.path.join(ERA5_BASE, r"ET\ET_transp_Annual\Sen_trend"),
        var_label="ET_transp_Annual",
        start_year=1982,
        end_year=2022,
        mk_pattern="{var_label}_mk_{y1}_{y2}.tif",
    ),

    # ---------- 年尺度 E_veg ---------- #
    dict(
        name="E_Annual",
        trend_dir=os.path.join(ERA5_BASE, r"ET\ET_Eveg_Annual\Sen_trend"),
        var_label="ET_Eveg_Annual",
        start_year=1982,
        end_year=2022,
        mk_pattern="{var_label}_mk_{y1}_{y2}.tif",
    ),

    # ---------- 春季 T ---------- #
    dict(
        name="T_Spring",
        trend_dir=os.path.join(ERA5_BASE, r"ET\ET_transp_Seasonal_1\Spring\Sen_trend"),
        var_label="ET_transp_Spring",
        start_year=1982,
        end_year=2022,
        mk_pattern="{var_label}_mk_{y1}_{y2}.tif",
    ),

    # ---------- 春季 E_veg ---------- #
    dict(
        name="E_Spring",
        trend_dir=os.path.join(ERA5_BASE, r"ET\ET_Eveg_Seasonal_1\Spring\Sen_trend"),
        var_label="ET_Eveg_Spring",
        start_year=1982,
        end_year=2022,
        mk_pattern="{var_label}_mk_{y1}_{y2}.tif",
    ),

    # ---------- 夏季 T ---------- #
    dict(
        name="T_Summer",
        trend_dir=os.path.join(ERA5_BASE, r"ET\ET_transp_Seasonal_1\Summer\Sen_trend"),
        var_label="ET_transp_Summer",
        start_year=1982,
        end_year=2022,
        mk_pattern="{var_label}_mk_{y1}_{y2}.tif",
    ),

    # ---------- 夏季 E_veg ---------- #
    dict(
        name="E_Summer",
        trend_dir=os.path.join(ERA5_BASE, r"ET\ET_Eveg_Seasonal_1\Summer\Sen_trend"),
        var_label="ET_Eveg_Summer",
        start_year=1982,
        end_year=2022,
        mk_pattern="{var_label}_mk_{y1}_{y2}.tif",
    ),

    # ---------- 秋季 T ---------- #
    dict(
        name="T_Autumn",
        trend_dir=os.path.join(ERA5_BASE, r"ET\ET_transp_Seasonal_1\Autumn\Sen_trend"),
        var_label="ET_transp_Autumn",
        start_year=1982,
        end_year=2022,
        mk_pattern="{var_label}_mk_{y1}_{y2}.tif",
    ),

    # ---------- 秋季 E_veg ---------- #
    dict(
        name="E_Autumn",
        trend_dir=os.path.join(ERA5_BASE, r"ET\ET_Eveg_Seasonal_1\Autumn\Sen_trend"),
        var_label="ET_Eveg_Autumn",
        start_year=1982,
        end_year=2022,
        mk_pattern="{var_label}_mk_{y1}_{y2}.tif",
    ),

    # ---------- 物候：SOS ---------- #
    dict(
        name="SOS",
        trend_dir=PHENO_SEN_DIR,
        var_label="sos",                     # sos_sen_1982_2022.tif
        start_year=1982,
        end_year=2022,
        mk_pattern="{var_label}_mk_{y1}_{y2}.tif",  # sos_mk_1982_2022.tif
    ),

    # ---------- 物候：POS ---------- #
    dict(
        name="POS",
        trend_dir=PHENO_SEN_DIR,
        var_label="pos",
        start_year=1982,
        end_year=2022,
        mk_pattern="{var_label}_mk_{y1}_{y2}.tif",
    ),

    # ---------- 物候：EOS ---------- #
    dict(
        name="EOS",
        trend_dir=PHENO_SEN_DIR,
        var_label="eos",
        start_year=1982,
        end_year=2022,
        mk_pattern="{var_label}_mk_{y1}_{y2}.tif",
    ),
]


def main():
    for task in TASKS:
        trend_dir = task["trend_dir"]
        var_label = task["var_label"]
        y1, y2 = task["start_year"], task["end_year"]

        # Sen trend 文件
        sen_file = os.path.join(trend_dir, f"{var_label}_sen_{y1}_{y2}.tif")

        # MK z 文件
        mk_pattern = task.get("mk_pattern", "{var_label}_mk_{y1}_{y2}.tif")
        mk_file = os.path.join(
            trend_dir,
            mk_pattern.format(var_label=var_label, y1=y1, y2=y2)
        )

        if not os.path.exists(sen_file):
            print(f"跳过 {task['name']}：找不到趋势栅格 {sen_file}")
            continue

        fig_name = f"{task['name']}_{y1}_{y2}_trend_with_significance.tif"
        out_figure = os.path.join(OUT_FIG_DIR, fig_name)

        plot_trend_map(sen_file, mk_file, out_figure, var_label)


if __name__ == "__main__":
    main()
