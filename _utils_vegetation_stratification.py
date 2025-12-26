#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
工具脚本: 按植被类型分层分析结果
用于在全局分析完成后，提取各植被类型的统计结果
"""

import numpy as np
import rasterio
from pathlib import Path
import pandas as pd
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# 导入配置
try:
    from _config import ROOT, LANDCOVER_FILE, NODATA_OUT
except ImportError:
    ROOT = Path(r"I:\F\Data4")
    LANDCOVER_FILE = ROOT / "Landcover" / "MCD12Q1" / "MCD12Q1_IGBP_2018.tif"
    NODATA_OUT = -9999.0

# 输出目录
OUTPUT_DIR = ROOT / "Wang2025_Analysis" / "Vegetation_Stratification"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# IGBP植被类型定义
VEG_TYPES = {
    1: 'ENF (常绿针叶林)',
    2: 'EBF (常绿阔叶林)',
    3: 'DNF (落叶针叶林)',
    4: 'DBF (落叶阔叶林)',
    5: 'MF (混交林)',
    6: 'CSH (郁闭灌木)',
    7: 'OSH (开阔灌木)',
    8: 'WSA (木本稀树草原)',
    9: 'SAV (稀树草原)',
    10: 'GRA (草地)',
    11: 'WET (永久湿地)',
    12: 'CRO (农田)',
    13: 'URB (城市)',
    14: 'CVM (农田/自然植被镶嵌)',
    15: 'SNI (雪/冰)',
    16: 'BSV (贫瘠)',
    17: 'WAT (水体)'
}

# 森林类型组
FOREST_TYPES = {
    'All_Forests': [1, 2, 3, 4, 5],
    'ENF': [1],
    'EBF': [2],
    'DNF': [3],
    'DBF': [4],
    'MF': [5],
    'Needleleaf_Forests': [1, 3],
    'Broadleaf_Forests': [2, 4],
    'Evergreen_Forests': [1, 2],
    'Deciduous_Forests': [3, 4]
}

# ==================== 核心函数 ====================
def extract_statistics_by_vegetation(data_file, landcover_file, veg_classes, mask_file=None):
    """
    按植被类型提取统计量

    Parameters:
    -----------
    data_file : Path
        数据文件路径（如TRc、趋势等）
    landcover_file : Path
        土地覆盖文件路径
    veg_classes : list
        植被类型代码列表
    mask_file : Path, optional
        额外的掩膜文件（如纬度掩膜）

    Returns:
    --------
    stats : dict
        统计结果字典
    """
    # 读取数据
    with rasterio.open(data_file) as src:
        data = src.read(1)

    # 读取土地覆盖
    with rasterio.open(landcover_file) as src:
        lc = src.read(1)

    # 读取掩膜（如果有）
    if mask_file and mask_file.exists():
        with rasterio.open(mask_file) as src:
            mask = src.read(1).astype(bool)
    else:
        mask = np.ones_like(data, dtype=bool)

    # 创建植被掩膜
    veg_mask = np.isin(lc, veg_classes)

    # 组合掩膜
    final_mask = mask & veg_mask & (data != NODATA_OUT) & ~np.isnan(data)

    # 提取有效值
    valid_data = data[final_mask]

    if len(valid_data) == 0:
        return {
            'count': 0,
            'mean': np.nan,
            'std': np.nan,
            'min': np.nan,
            'q25': np.nan,
            'median': np.nan,
            'q75': np.nan,
            'max': np.nan
        }

    # 计算统计量
    stats = {
        'count': len(valid_data),
        'mean': np.mean(valid_data),
        'std': np.std(valid_data),
        'min': np.min(valid_data),
        'q25': np.percentile(valid_data, 25),
        'median': np.median(valid_data),
        'q75': np.percentile(valid_data, 75),
        'max': np.max(valid_data)
    }

    return stats

# ==================== 批量处理 ====================
def analyze_variable_by_vegetation(var_name, file_pattern, veg_types_dict=FOREST_TYPES):
    """
    按植被类型分析某个变量

    Parameters:
    -----------
    var_name : str
        变量名（如'TRc', 'TRpheno', 'trend'）
    file_pattern : str or Path
        文件路径模式（如'Decomposition/TRc_av.tif'）
    veg_types_dict : dict
        植被类型分组字典

    Returns:
    --------
    results_df : DataFrame
        统计结果表
    """
    print(f"\n分析变量: {var_name}")

    # 构建文件路径
    if isinstance(file_pattern, str):
        data_file = ROOT / "Wang2025_Analysis" / file_pattern
    else:
        data_file = file_pattern

    if not data_file.exists():
        print(f"  ✗ 文件不存在: {data_file}")
        return None

    # 纬度掩膜
    mask_file = ROOT / "Wang2025_Analysis" / "masks" / "lat_mask.tif"

    # 逐植被类型提取统计
    results = []

    for veg_name, veg_classes in tqdm(veg_types_dict.items(), desc="植被类型"):
        stats = extract_statistics_by_vegetation(
            data_file, LANDCOVER_FILE, veg_classes, mask_file
        )

        results.append({
            'vegetation_type': veg_name,
            'veg_classes': str(veg_classes),
            **stats
        })

    # 转为DataFrame
    results_df = pd.DataFrame(results)

    # 保存结果
    output_file = OUTPUT_DIR / f"{var_name}_by_vegetation.csv"
    results_df.to_csv(output_file, index=False, float_format='%.4f')
    print(f"  ✓ 已保存: {output_file.name}")

    return results_df

# ==================== 时间序列分析 ====================
def analyze_timeseries_by_vegetation(var_name, year_start, year_end, file_template):
    """
    分析时间序列数据的植被分层统计

    Parameters:
    -----------
    var_name : str
        变量名（如'TRc', 'TRpheno'）
    year_start, year_end : int
        年份范围
    file_template : str
        文件名模板（如'TRc_annual/TRc_{year}.tif'）

    Returns:
    --------
    results_df : DataFrame
        时间序列统计结果
    """
    print(f"\n分析时间序列: {var_name} ({year_start}-{year_end})")

    years = list(range(year_start, year_end + 1))
    mask_file = ROOT / "Wang2025_Analysis" / "masks" / "lat_mask.tif"

    all_results = []

    for year in tqdm(years, desc="年份"):
        # 构建文件路径
        file_path = ROOT / "Wang2025_Analysis" / file_template.format(year=year)

        if not file_path.exists():
            continue

        # 逐植被类型提取
        for veg_name, veg_classes in FOREST_TYPES.items():
            stats = extract_statistics_by_vegetation(
                file_path, LANDCOVER_FILE, veg_classes, mask_file
            )

            all_results.append({
                'year': year,
                'vegetation_type': veg_name,
                **stats
            })

    # 转为DataFrame
    results_df = pd.DataFrame(all_results)

    # 保存结果
    output_file = OUTPUT_DIR / f"{var_name}_timeseries_by_vegetation.csv"
    results_df.to_csv(output_file, index=False, float_format='%.4f')
    print(f"  ✓ 已保存: {output_file.name}")

    return results_df

# ==================== 创建分层掩膜 ====================
def create_vegetation_masks():
    """
    为每种植被类型创建单独的掩膜文件
    """
    print("\n创建植被分层掩膜...")

    # 读取土地覆盖
    with rasterio.open(LANDCOVER_FILE) as src:
        lc = src.read(1)
        profile = src.profile.copy()

    # 读取纬度掩膜
    mask_file = ROOT / "Wang2025_Analysis" / "masks" / "lat_mask.tif"
    with rasterio.open(mask_file) as src:
        lat_mask = src.read(1).astype(bool)

    # 为每个植被类型组创建掩膜
    veg_mask_dir = OUTPUT_DIR / "vegetation_masks"
    veg_mask_dir.mkdir(exist_ok=True)

    profile.update(dtype=rasterio.uint8, count=1, compress='lzw', nodata=0)

    for veg_name, veg_classes in tqdm(FOREST_TYPES.items(), desc="创建掩膜"):
        # 创建植被掩膜
        veg_mask = np.isin(lc, veg_classes)
        combined_mask = (lat_mask & veg_mask).astype(np.uint8)

        # 保存
        output_file = veg_mask_dir / f"mask_{veg_name}.tif"
        with rasterio.open(output_file, 'w', **profile) as dst:
            dst.write(combined_mask, 1)

        pixel_count = np.sum(combined_mask)
        print(f"  {veg_name}: {pixel_count} 像元")

    print(f"  ✓ 掩膜已保存至: {veg_mask_dir}")

# ==================== 主程序 ====================
def main():
    print("\n" + "="*70)
    print("按植被类型分层分析")
    print("="*70)

    # 检查土地覆盖文件
    if not LANDCOVER_FILE.exists():
        print(f"\n✗ 错误：土地覆盖文件不存在")
        print(f"  路径: {LANDCOVER_FILE}")
        print("\n请提供MODIS IGBP土地覆盖数据")
        return

    print(f"\n土地覆盖文件: {LANDCOVER_FILE.name}")

    # 1. 创建植被分层掩膜
    create_vegetation_masks()

    # 2. 分析多年平均值
    print("\n" + "="*70)
    print("分析多年平均值")
    print("="*70)

    variables_avg = {
        'TRc_av': 'Decomposition/TRc_av.tif',
        'LSP_av': 'Decomposition/LSP_av.tif',
    }

    for var_name, file_path in variables_avg.items():
        analyze_variable_by_vegetation(var_name, file_path)

    # 3. 分析趋势
    print("\n" + "="*70)
    print("分析趋势")
    print("="*70)

    variables_trend = {
        'SOS_trend': 'Statistical_Analysis/Trends/SOS_slope.tif',
        'TRpheno_trend': 'Statistical_Analysis/Trends/TRpheno_slope.tif',
        'TRproduct_trend': 'Statistical_Analysis/Trends/TRproduct_slope.tif',
    }

    for var_name, file_path in variables_trend.items():
        analyze_variable_by_vegetation(var_name, file_path)

    # 4. 分析时间序列（可选）
    print("\n" + "="*70)
    print("分析时间序列（可选，较耗时）")
    print("="*70)

    user_input = input("是否分析时间序列数据？这可能需要较长时间。[y/N]: ")
    if user_input.lower() == 'y':
        analyze_timeseries_by_vegetation('TRc', 1982, 2018, 'TRc_annual/TRc_{year}.tif')
        analyze_timeseries_by_vegetation('TRpheno', 1982, 2018, 'Decomposition/TRpheno_{year}.tif')
        analyze_timeseries_by_vegetation('TRproduct', 1982, 2018, 'Decomposition/TRproduct_{year}.tif')

    print("\n" + "="*70)
    print("✓ 植被分层分析完成！")
    print(f"输出目录: {OUTPUT_DIR}")
    print("\n文件说明:")
    print("  - *_by_vegetation.csv        各植被类型统计摘要")
    print("  - *_timeseries_by_vegetation.csv  时间序列统计")
    print("  - vegetation_masks/          各植被类型掩膜")
    print("="*70)

# ==================== 运行 ====================
if __name__ == "__main__":
    main()
