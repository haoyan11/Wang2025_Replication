#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module 07: 绘图函数
生成Wang (2025)论文风格的所有图表
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import rasterio
from pathlib import Path
from scipy import stats
import seaborn as sns
import pandas as pd

# ==================== 全局配置 ====================
ROOT = Path(r"I:\F\Data4")
ANALYSIS_DIR = ROOT / "Wang2025_Analysis"
OUTPUT_FIG_DIR = ANALYSIS_DIR / "Figures"
OUTPUT_FIG_DIR.mkdir(parents=True, exist_ok=True)

# 设置matplotlib中文字体
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 10

# ==================== 辅助函数 ====================
def read_geotiff(file_path):
    """读取GeoTIFF并返回数据和坐标"""
    with rasterio.open(file_path) as src:
        data = src.read(1)
        bounds = src.bounds
        transform = src.transform

        # 生成经纬度网格
        height, width = data.shape
        cols, rows = np.meshgrid(np.arange(width), np.arange(height))
        lon, lat = rasterio.transform.xy(transform, rows, cols)
        lon = np.array(lon)
        lat = np.array(lat)

    return data, lon, lat, bounds

def mask_nodata(data, nodata=-9999.0):
    """将NODATA转换为masked array"""
    return np.ma.masked_equal(data, nodata)

# ==================== 配色方案 ====================
def get_diverging_colormap(name='RdYlBu_r'):
    """获取发散色标"""
    if name == 'RdYlBu_r':
        return plt.cm.RdYlBu_r
    elif name == 'BrBG':
        return plt.cm.BrBG
    else:
        return plt.cm.RdBu_r

def get_sequential_colormap(name='YlGn'):
    """获取连续色标"""
    if name == 'YlGn':
        return plt.cm.YlGn
    elif name == 'Blues':
        return plt.cm.Blues
    else:
        return plt.cm.viridis

# ==================== 地图绘制基础函数 ====================
def create_map_axes(ax, extent=None):
    """创建地图坐标系"""
    if extent is None:
        extent = [-180, 180, 30, 90]  # 默认≥30°N

    ax.set_extent(extent, crs=ccrs.PlateCarree())

    # 添加地理要素
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle='--', alpha=0.5)

    # 添加经纬度网格
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray',
                      alpha=0.3, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}

    return ax

def plot_spatial_map(data, lon, lat, ax, cmap, vmin, vmax, title,
                     cbar_label, extend='both', mask_ocean=True):
    """
    绘制空间分布图

    Parameters:
    -----------
    data : ndarray
        数据数组
    lon, lat : ndarray
        经纬度坐标
    ax : matplotlib.axes
        子图对象
    cmap : colormap
        色标
    vmin, vmax : float
        色标范围
    title : str
        标题
    cbar_label : str
        色标标签
    extend : str
        色标延伸方向
    mask_ocean : bool
        是否遮罩海洋
    """
    # 遮罩NODATA
    data_masked = mask_nodata(data)

    # 绘制
    im = ax.pcolormesh(lon, lat, data_masked,
                       cmap=cmap, vmin=vmin, vmax=vmax,
                       transform=ccrs.PlateCarree(),
                       shading='auto')

    # 添加色标
    cbar = plt.colorbar(im, ax=ax, orientation='horizontal',
                        pad=0.05, shrink=0.8, extend=extend)
    cbar.set_label(cbar_label, fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    # 添加标题
    ax.set_title(title, fontsize=11, fontweight='bold')

    return im

# ==================== Figure 1: 研究区域与数据概览 ====================
def plot_figure1_study_area():
    """绘制研究区域与森林掩膜"""
    print("\n绘制 Figure 1: 研究区域...")

    fig = plt.figure(figsize=(12, 8))

    # (a) 森林掩膜
    ax1 = fig.add_subplot(2, 2, 1, projection=ccrs.Robinson())
    mask_file = ANALYSIS_DIR / "masks" / "forest_mask.tif"
    data, lon, lat, _ = read_geotiff(mask_file)
    create_map_axes(ax1, extent=[-180, 180, 30, 90])

    data_masked = np.ma.masked_equal(data, 0)
    ax1.pcolormesh(lon, lat, data_masked, cmap='Greens',
                   transform=ccrs.PlateCarree(), shading='auto')
    ax1.set_title("(a) Forest mask (≥30°N)", fontsize=11, fontweight='bold')

    # (b) 多年平均TRc
    ax2 = fig.add_subplot(2, 2, 2, projection=ccrs.Robinson())
    TRc_av_file = ANALYSIS_DIR / "Decomposition" / "TRc_av.tif"
    data, lon, lat, _ = read_geotiff(TRc_av_file)
    create_map_axes(ax2, extent=[-180, 180, 30, 90])

    plot_spatial_map(data, lon, lat, ax2,
                     cmap=get_sequential_colormap('Blues'),
                     vmin=0, vmax=500,
                     title="(b) Mean TRc (1982-2018)",
                     cbar_label="TRc (mm)")

    # (c) 多年平均LSP
    ax3 = fig.add_subplot(2, 2, 3, projection=ccrs.Robinson())
    LSP_av_file = ANALYSIS_DIR / "Decomposition" / "LSP_av.tif"
    data, lon, lat, _ = read_geotiff(LSP_av_file)
    create_map_axes(ax3, extent=[-180, 180, 30, 90])

    plot_spatial_map(data, lon, lat, ax3,
                     cmap=get_sequential_colormap('YlGn'),
                     vmin=60, vmax=150,
                     title="(c) Mean LSP (1982-2018)",
                     cbar_label="LSP (days)")

    # (d) 多年平均SOS
    ax4 = fig.add_subplot(2, 2, 4, projection=ccrs.Robinson())
    SOS_file = ROOT / "Phenology_Output_1" / "SIF_phenology" / "SOS_2000.tif"  # 示例年份
    data, lon, lat, _ = read_geotiff(SOS_file)
    create_map_axes(ax4, extent=[-180, 180, 30, 90])

    plot_spatial_map(data, lon, lat, ax4,
                     cmap=get_sequential_colormap('YlOrRd'),
                     vmin=60, vmax=150,
                     title="(d) Mean SOS (1982-2018)",
                     cbar_label="SOS (DOY)")

    plt.tight_layout()
    plt.savefig(OUTPUT_FIG_DIR / "Figure1_StudyArea.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("  ✓ Figure 1 已保存")

# ==================== Figure 2: 趋势分析 ====================
def plot_figure2_trends():
    """绘制TRpheno和TRproduct的趋势"""
    print("\n绘制 Figure 2: 趋势分析...")

    fig = plt.figure(figsize=(14, 6))

    trend_dir = ANALYSIS_DIR / "Statistical_Analysis" / "Trends"

    # (a) TRpheno趋势
    ax1 = fig.add_subplot(1, 2, 1, projection=ccrs.Robinson())
    slope_file = trend_dir / "TRpheno_slope.tif"
    data, lon, lat, _ = read_geotiff(slope_file)
    create_map_axes(ax1, extent=[-180, 180, 30, 90])

    # 转换为mm/decade
    data_decade = data * 10

    plot_spatial_map(data_decade, lon, lat, ax1,
                     cmap=get_diverging_colormap('RdBu_r'),
                     vmin=-20, vmax=20,
                     title="(a) TRpheno trend (1982-2018)",
                     cbar_label="Trend (mm/decade)",
                     extend='both')

    # (b) TRproduct趋势
    ax2 = fig.add_subplot(1, 2, 2, projection=ccrs.Robinson())
    slope_file = trend_dir / "TRproduct_slope.tif"
    data, lon, lat, _ = read_geotiff(slope_file)
    create_map_axes(ax2, extent=[-180, 180, 30, 90])

    data_decade = data * 10

    plot_spatial_map(data_decade, lon, lat, ax2,
                     cmap=get_diverging_colormap('RdBu_r'),
                     vmin=-20, vmax=20,
                     title="(b) TRproduct trend (1982-2018)",
                     cbar_label="Trend (mm/decade)",
                     extend='both')

    plt.tight_layout()
    plt.savefig(OUTPUT_FIG_DIR / "Figure2_Trends.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("  ✓ Figure 2 已保存")

# ==================== Figure 3: 归因分析 ====================
def plot_figure3_attribution():
    """绘制归因分析结果"""
    print("\n绘制 Figure 3: 归因分析...")

    fig = plt.figure(figsize=(16, 10))

    attr_dir = ANALYSIS_DIR / "Statistical_Analysis" / "Attribution" / "TRproduct"

    variables = ['SOS', 'LSP', 'SIF', 'SM']
    titles = ['(a) SOS contribution', '(b) LSP contribution',
              '(c) SIF contribution', '(d) SM contribution']

    for idx, (var, title) in enumerate(zip(variables, titles), 1):
        ax = fig.add_subplot(2, 2, idx, projection=ccrs.Robinson())
        attr_file = attr_dir / f"attribution_{var}.tif"

        if not attr_file.exists():
            continue

        data, lon, lat, _ = read_geotiff(attr_file)
        create_map_axes(ax, extent=[-180, 180, 30, 90])

        plot_spatial_map(data, lon, lat, ax,
                         cmap=get_diverging_colormap('RdBu_r'),
                         vmin=-0.5, vmax=0.5,
                         title=title,
                         cbar_label="Standardized coefficient",
                         extend='both')

    plt.tight_layout()
    plt.savefig(OUTPUT_FIG_DIR / "Figure3_Attribution.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("  ✓ Figure 3 已保存")

# ==================== Figure 4: 滑动窗口敏感性 ====================
def plot_figure4_moving_window():
    """绘制15年滑动窗口敏感性"""
    print("\n绘制 Figure 4: 滑动窗口敏感性...")

    sens_dir = ANALYSIS_DIR / "Statistical_Analysis" / "Moving_Window_Sensitivity"

    # 获取所有窗口文件
    sens_files = sorted(sens_dir.glob("sensitivity_TRc_SOS_*.tif"))

    if len(sens_files) == 0:
        print("  ⚠ 警告：未找到滑动窗口结果文件")
        return

    # 读取数据计算空间平均
    years_mid = []
    mean_sensitivity = []

    for file in sens_files:
        # 从文件名提取年份
        filename = file.stem
        year_range = filename.split('_')[-1]
        start_year, end_year = map(int, year_range.split('-'))
        mid_year = (start_year + end_year) / 2

        data, _, _, _ = read_geotiff(file)
        data_masked = mask_nodata(data)
        mean_val = np.ma.mean(data_masked)

        years_mid.append(mid_year)
        mean_sensitivity.append(mean_val)

    # 绘制时间序列
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.plot(years_mid, mean_sensitivity, 'o-', linewidth=2, markersize=6,
            color='#2E86AB', label='15-year moving window')

    # 添加零线
    ax.axhline(0, color='gray', linestyle='--', linewidth=1, alpha=0.5)

    ax.set_xlabel('Year (window center)', fontsize=12)
    ax.set_ylabel('Sensitivity (TRc ~ SOS correlation)', fontsize=12)
    ax.set_title('Temporal variation of TRc sensitivity to SOS', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(OUTPUT_FIG_DIR / "Figure4_MovingWindow.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("  ✓ Figure 4 已保存")

# ==================== Figure 5: 分解方法比较 ====================
def plot_figure5_decomposition_comparison():
    """比较原版vs改进版分解方法"""
    print("\n绘制 Figure 5: 分解方法比较...")

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    # 原版方法
    decomp_original = ANALYSIS_DIR / "Decomposition"
    # 改进方法
    decomp_improved = ANALYSIS_DIR / "Decomposition_Improved" / "Method2_Rate"

    test_year = 2000

    # Row 1: 原版方法
    for idx, (var, title) in enumerate([
        ('TRpheno', '(a) TRpheno - Original'),
        ('TRproduct', '(b) TRproduct - Original'),
        ('TRc', '(c) TRc reconstruction')
    ]):
        ax = axes[0, idx]

        if var == 'TRc':
            file_path = ANALYSIS_DIR / "TRc_annual" / f"TRc_{test_year}.tif"
        else:
            file_path = decomp_original / f"{var}_{test_year}.tif"

        data, _, _, _ = read_geotiff(file_path)
        data_masked = mask_nodata(data)

        im = ax.imshow(data_masked, cmap='RdBu_r', vmin=-100, vmax=100)
        ax.set_title(title, fontsize=11, fontweight='bold')
        ax.axis('off')
        plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, shrink=0.8)

    # Row 2: 改进方法
    for idx, (var, title) in enumerate([
        ('TRpheno_rate', '(d) TRpheno - Rate-based'),
        ('TRproduct_rate', '(e) TRproduct - Rate-based'),
        ('diff', '(f) Difference (Original - Rate)')
    ]):
        ax = axes[1, idx]

        if var == 'diff':
            # 计算差异
            orig, _, _, _ = read_geotiff(decomp_original / f"TRpheno_{test_year}.tif")
            rate, _, _, _ = read_geotiff(decomp_improved / f"TRpheno_rate_{test_year}.tif")
            data = orig - rate
        else:
            file_path = decomp_improved / f"{var}_{test_year}.tif"
            data, _, _, _ = read_geotiff(file_path)

        data_masked = mask_nodata(data)

        im = ax.imshow(data_masked, cmap='RdBu_r', vmin=-100, vmax=100)
        ax.set_title(title, fontsize=11, fontweight='bold')
        ax.axis('off')
        plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.05, shrink=0.8)

    plt.tight_layout()
    plt.savefig(OUTPUT_FIG_DIR / "Figure5_DecompositionComparison.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("  ✓ Figure 5 已保存")

# ==================== Figure 6: SEM路径图（从R导出后展示） ====================
def plot_figure6_SEM_summary():
    """展示SEM结果摘要（基于R输出的CSV）"""
    print("\n绘制 Figure 6: SEM结果摘要...")

    sem_dir = ANALYSIS_DIR / "SEM_Results"

    # 读取参数估计
    params_original = pd.read_csv(sem_dir / "SEM_original_parameters.csv")
    params_improved = pd.read_csv(sem_dir / "SEM_improved_parameters.csv")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # (a) 原版SEM路径系数
    ax1 = axes[0]
    paths_original = params_original[params_original['op'] == '~']
    y_pos = np.arange(len(paths_original))

    ax1.barh(y_pos, paths_original['std.all'], color='steelblue', alpha=0.7)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels([f"{row['lhs']} ~ {row['rhs']}" for _, row in paths_original.iterrows()])
    ax1.set_xlabel('Standardized coefficient', fontsize=11)
    ax1.set_title('(a) Original SEM pathways', fontsize=12, fontweight='bold')
    ax1.axvline(0, color='gray', linestyle='--', linewidth=1)
    ax1.grid(axis='x', alpha=0.3)

    # (b) 改进版SEM路径系数
    ax2 = axes[1]
    paths_improved = params_improved[params_improved['op'] == '~']
    y_pos = np.arange(len(paths_improved))

    ax2.barh(y_pos, paths_improved['std.all'], color='darkorange', alpha=0.7)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels([f"{row['lhs']} ~ {row['rhs']}" for _, row in paths_improved.iterrows()])
    ax2.set_xlabel('Standardized coefficient', fontsize=11)
    ax2.set_title('(b) Improved SEM pathways', fontsize=12, fontweight='bold')
    ax2.axvline(0, color='gray', linestyle='--', linewidth=1)
    ax2.grid(axis='x', alpha=0.3)

    plt.tight_layout()
    plt.savefig(OUTPUT_FIG_DIR / "Figure6_SEM_Summary.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("  ✓ Figure 6 已保存")

# ==================== 主程序 ====================
def main():
    print("\n" + "="*70)
    print("绘制所有图表")
    print("="*70)

    plot_figure1_study_area()
    plot_figure2_trends()
    plot_figure3_attribution()
    plot_figure4_moving_window()
    plot_figure5_decomposition_comparison()
    plot_figure6_SEM_summary()

    print("\n" + "="*70)
    print("✓ 所有图表绘制完成！")
    print(f"输出目录: {OUTPUT_FIG_DIR}")
    print("="*70)

if __name__ == "__main__":
    main()
