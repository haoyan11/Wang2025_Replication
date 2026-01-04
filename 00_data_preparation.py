#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module 00: 数据准备
创建分析所需的掩膜文件和检查数据完整性
1. 创建纬度掩膜（≥30°N）
2. [可选] 创建森林类型掩膜（IGBP分类1-5）
3. 合并掩膜
4. 检查数据完整性
"""

import numpy as np
import rasterio
from pathlib import Path
from datetime import datetime, timedelta
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# 导入配置
try:
    # 优先尝试导入_config.py
    from _config import (
        ROOT, LAT_MIN, FOREST_CLASSES, NODATA_OUT,
        USE_FOREST_MASK, LANDCOVER_FILE, TR_DAILY_DIR, PHENO_DIR
    )
except ImportError:
    try:
        # 如果_config不存在，尝试config.py
        from config import (
            ROOT, LAT_MIN, FOREST_CLASSES, NODATA_OUT,
            USE_FOREST_MASK, LANDCOVER_FILE, TR_DAILY_DIR, PHENO_DIR
        )
    except ImportError:
        # 如果都无法导入，使用默认配置
        print("⚠ 警告：无法导入_config.py或config.py，使用默认配置")
        ROOT = Path(r"I:\F\Data4")
        LAT_MIN = 30.0
        FOREST_CLASSES = [1, 2, 3, 4, 5]
        NODATA_OUT = -9999.0
        USE_FOREST_MASK = False  # 默认不使用森林掩膜
        LANDCOVER_FILE = ROOT / "Landcover" / "MCD12Q1" / "MCD12Q1_IGBP_2018.tif"
        TR_DAILY_DIR = ROOT / "Meteorological Data" / "ERA5_Land" / "ET_components" / "ET_transp" / "ET_transp_Daily" / "ET_transp_Daily_2"
        PHENO_DIR = ROOT / "Phenology_Output_1" / "GPP_phenology_EPSG4326"

# 输出目录
OUTPUT_DIR = ROOT / "Wang2025_Analysis" / "masks"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ==================== 辅助函数 ====================
def get_template_file():
    """
    获取模板文件（用于生成掩膜的空间参考）
    优先级：TR数据 > 物候（确保使用EPSG:4326 CRS）
    """
    # 优先尝试从TR数据获取（CRS: EPSG:4326）
    # ERA5-Land格式: ERA5L_ET_transp_Daily_mm_YYYYMMDD.tif
    test_date = datetime(2000, 1, 15)
    yyyymmdd = test_date.strftime('%Y%m%d')
    tr_files = [
        TR_DAILY_DIR / f"ERA5L_ET_transp_Daily_mm_{yyyymmdd}.tif",  # ERA5-Land格式
        TR_DAILY_DIR / f"Et_{yyyymmdd}.tif",                         # GLEAM格式（备选）
        TR_DAILY_DIR / f"Et_{test_date.year}_{test_date.timetuple().tm_yday:03d}.tif",
        TR_DAILY_DIR / str(test_date.year) / f"Et_{yyyymmdd}.tif"
    ]

    for f in tr_files:
        if f.exists():
            return f

    # 如果找不到TR数据，尝试从物候数据获取
    pheno_files = list(PHENO_DIR.glob("sos_gpp_*.tif"))
    if pheno_files:
        return pheno_files[0]

    return None

# ==================== 1. 创建纬度掩膜 ====================
def create_lat_mask(template_file, lat_min=30.0):
    """
    创建纬度掩膜

    Parameters:
    -----------
    template_file : Path
        模板文件路径
    lat_min : float
        最小纬度（度）

    Returns:
    --------
    mask : ndarray
        纬度掩膜（布尔数组）
    profile : dict
        栅格配置
    """
    print(f"\n[1] 创建纬度掩膜（≥{lat_min}°N）...")

    with rasterio.open(template_file) as src:
        profile = src.profile.copy()
        height, width = src.height, src.width
        transform = src.transform
        crs = src.crs

    # CRS检查：确保是地理坐标系（EPSG:4326或类似）
    if crs is not None:
        if not crs.is_geographic:
            print(f"  ⚠ 警告：模板CRS不是地理坐标系: {crs}")
            print(f"  纬度掩膜可能不准确，建议使用EPSG:4326数据")
        elif str(crs) != 'EPSG:4326':
            print(f"  ⚠ 注意：模板CRS为 {crs}，非标准EPSG:4326")
            print(f"  纬度计算可能略有偏差")

    # 生成纬度数组
    # 对于每一行，计算中心点纬度
    lat_array = np.zeros((height, width), dtype=np.float32)

    for row in range(height):
        # 计算该行的地理纬度
        lon, lat = rasterio.transform.xy(transform, row, 0, offset='center')
        lat_array[row, :] = lat

    # 创建掩膜
    mask = lat_array >= lat_min

    print(f"  栅格尺寸: {height} × {width}")
    print(f"  纬度范围: {np.min(lat_array):.2f}° - {np.max(lat_array):.2f}°")
    print(f"  有效像元数: {np.sum(mask)} ({np.sum(mask) / (height * width) * 100:.1f}%)")

    return mask.astype(np.uint8), profile

# ==================== 2. 创建森林掩膜 ====================
def create_forest_mask(landcover_file, forest_classes=[1, 2, 3, 4, 5]):
    """
    创建森林类型掩膜

    Parameters:
    -----------
    landcover_file : Path
        MODIS土地覆盖文件路径
    forest_classes : list
        森林类型代码列表
        1: Evergreen Needleleaf Forest (ENF)
        2: Evergreen Broadleaf Forest (EBF)
        3: Deciduous Needleleaf Forest (DNF)
        4: Deciduous Broadleaf Forest (DBF)
        5: Mixed Forest (MF)

    Returns:
    --------
    mask : ndarray
        森林掩膜（布尔数组）
    profile : dict
        栅格配置
    """
    print(f"\n[2] 创建森林类型掩膜...")

    if not landcover_file.exists():
        print(f"  ⚠ 警告：土地覆盖文件不存在: {landcover_file}")
        print("  将创建全1掩膜（跳过森林类型筛选）")

        # 使用模板创建全1掩膜
        template = get_template_file()
        with rasterio.open(template) as src:
            profile = src.profile.copy()
            height, width = src.height, src.width

        mask = np.ones((height, width), dtype=np.uint8)
        return mask, profile

    with rasterio.open(landcover_file) as src:
        lc_data = src.read(1)
        profile = src.profile.copy()

    # 创建森林掩膜
    mask = np.isin(lc_data, forest_classes)

    print(f"  森林类型: {forest_classes}")
    print(f"  有效像元数: {np.sum(mask)} ({np.sum(mask) / mask.size * 100:.1f}%)")

    # 统计各森林类型
    for fc in forest_classes:
        count = np.sum(lc_data == fc)
        fc_names = {
            1: 'ENF (常绿针叶林)',
            2: 'EBF (常绿阔叶林)',
            3: 'DNF (落叶针叶林)',
            4: 'DBF (落叶阔叶林)',
            5: 'MF (混交林)'
        }
        print(f"    {fc_names.get(fc, f'类型{fc}')}: {count} 像元")

    return mask.astype(np.uint8), profile

# ==================== 3. 合并掩膜 ====================
def combine_masks(lat_mask, forest_mask):
    """
    合并纬度掩膜和森林掩膜

    Parameters:
    -----------
    lat_mask : ndarray
        纬度掩膜
    forest_mask : ndarray
        森林掩膜

    Returns:
    --------
    combined_mask : ndarray
        组合掩膜
    """
    print(f"\n[3] 合并掩膜...")

    # 确保尺寸一致
    if lat_mask.shape != forest_mask.shape:
        print(f"  ⚠ 警告：掩膜尺寸不一致")
        print(f"    纬度掩膜: {lat_mask.shape}")
        print(f"    森林掩膜: {forest_mask.shape}")

        # 使用最小公共尺寸
        min_h = min(lat_mask.shape[0], forest_mask.shape[0])
        min_w = min(lat_mask.shape[1], forest_mask.shape[1])

        lat_mask = lat_mask[:min_h, :min_w]
        forest_mask = forest_mask[:min_h, :min_w]

        print(f"    裁剪后: {lat_mask.shape}")

    # 合并（逻辑与）
    combined_mask = (lat_mask.astype(bool) & forest_mask.astype(bool)).astype(np.uint8)

    print(f"  纬度掩膜有效像元: {np.sum(lat_mask)}")
    print(f"  森林掩膜有效像元: {np.sum(forest_mask)}")
    print(f"  组合掩膜有效像元: {np.sum(combined_mask)}")
    print(f"  保留比例: {np.sum(combined_mask) / np.sum(lat_mask) * 100:.1f}%")

    return combined_mask

# ==================== 4. 数据完整性检查 ====================
def check_data_availability(years=[2000, 2010, 2018]):
    """
    检查关键数据的可用性

    Parameters:
    -----------
    years : list
        检查的年份列表
    """
    print(f"\n[4] 数据完整性检查...")

    results = {
        'pheno': {'available': 0, 'missing': 0},
        'tr': {'available': 0, 'missing': 0}
    }

    # 检查物候数据
    print("\n  物候数据 (SOS/POS/EOS):")
    for year in years:
        # 根据配置的格式检查物候文件
        sos_file = PHENO_DIR / f"sos_gpp_{year}.tif"  # 物候代码输出小写
        pos_file = PHENO_DIR / f"pos_doy_gpp_{year}.tif"
        eos_file = PHENO_DIR / f"eos_gpp_{year}.tif"

        files_exist = [
            sos_file.exists(),
            pos_file.exists(),
            eos_file.exists()
        ]

        if all(files_exist):
            print(f"    {year}: ✓")
            results['pheno']['available'] += 1
        else:
            # 修复警告：先计算缺失列表再插入f-string
            missing_types = [name for i, name in enumerate(['SOS', 'POS', 'EOS']) if not files_exist[i]]
            print(f"    {year}: ✗ (缺失: {missing_types})")
            results['pheno']['missing'] += 1

    # 检查TR数据（抽样检查每年1月15日）
    print("\n  TR数据 (每年1月15日抽样 - ERA5-Land):")
    for year in years:
        test_date = datetime(year, 1, 15)
        # ERA5-Land格式: ERA5L_ET_transp_Daily_mm_YYYYMMDD.tif
        tr_file = TR_DAILY_DIR / f"ERA5L_ET_transp_Daily_mm_{test_date.strftime('%Y%m%d')}.tif"

        if tr_file.exists():
            print(f"    {year}: ✓ ({tr_file.name})")
            results['tr']['available'] += 1
        else:
            print(f"    {year}: ✗")
            results['tr']['missing'] += 1

    # 总结
    print("\n  总结:")
    print(f"    物候数据: {results['pheno']['available']}/{len(years)} 年份可用")
    print(f"    TR数据:   {results['tr']['available']}/{len(years)} 年份可用")

    if results['pheno']['missing'] > 0:
        print(f"\n  ⚠ 提示：缺少物候数据，需运行 Module 01: 01_phenology_extraction.py")

    if results['tr']['missing'] > 0:
        print(f"\n  ⚠ 提示：缺少TR数据，请检查路径: {TR_DAILY_DIR}")

    return results

# ==================== 保存掩膜 ====================
def save_mask(file_path, mask, profile):
    """保存掩膜为GeoTIFF"""
    profile.update(
        dtype=rasterio.uint8,
        count=1,
        compress='lzw',
        nodata=0
    )

    with rasterio.open(file_path, 'w', **profile) as dst:
        dst.write(mask.astype(np.uint8), 1)

    print(f"  ✓ 已保存: {file_path.name}")

# ==================== 主程序 ====================
def main():
    print("\n" + "="*70)
    print("Module 00: 数据准备")
    print("="*70)

    # 显示掩膜策略
    print(f"\n掩膜策略: {'仅森林类型' if USE_FOREST_MASK else '所有陆地（≥30°N）'}")
    if not USE_FOREST_MASK:
        print("  → 将仅使用纬度掩膜，不筛选森林类型")
        print("  → 后续可按植被类型分层分析")

    # 获取模板文件
    print("\n获取模板文件...")
    template_file = get_template_file()

    if template_file is None:
        print("  ✗ 错误：找不到模板文件")
        print("  请确保以下路径之一存在数据:")
        print(f"    物候数据: {PHENO_DIR}")
        print(f"    TR数据: {TR_DAILY_DIR}")
        return

    print(f"  ✓ 使用模板: {template_file.name}")

    # 1. 创建纬度掩膜
    lat_mask, profile = create_lat_mask(template_file, lat_min=LAT_MIN)
    save_mask(OUTPUT_DIR / "lat_mask.tif", lat_mask, profile)

    # 2. 创建森林掩膜（可选）
    if USE_FOREST_MASK:
        print("\n  使用森林掩膜模式...")
        forest_mask, forest_profile = create_forest_mask(LANDCOVER_FILE, FOREST_CLASSES)

        # 确保使用相同的profile
        if forest_profile['width'] != profile['width'] or forest_profile['height'] != profile['height']:
            print("\n  ⚠ 警告：森林掩膜与模板尺寸不同，将重采样...")
            # 这里简化处理，实际应用中需要使用rasterio.warp.reproject
            # 暂时使用模板的profile
        else:
            profile = forest_profile

        save_mask(OUTPUT_DIR / "forest_mask.tif", forest_mask, profile)

        # 3. 合并掩膜
        combined_mask = combine_masks(lat_mask, forest_mask)
    else:
        print("\n  跳过森林掩膜（USE_FOREST_MASK=False）")
        print("  组合掩膜将等于纬度掩膜")

        # 创建一个占位的森林掩膜（全1）用于记录
        forest_mask = np.ones_like(lat_mask)
        save_mask(OUTPUT_DIR / "forest_mask.tif", forest_mask, profile)

        # 组合掩膜直接等于纬度掩膜
        combined_mask = lat_mask.copy()

    save_mask(OUTPUT_DIR / "combined_mask.tif", combined_mask, profile)

    # 4. 数据完整性检查
    check_data_availability(years=[1982, 2000, 2010, 2018])

    print("\n" + "="*70)
    print("✓ 数据准备完成！")
    print(f"输出目录: {OUTPUT_DIR}")
    print("\n掩膜文件:")
    print(f"  - lat_mask.tif       (纬度 ≥{LAT_MIN}°N)")
    if USE_FOREST_MASK:
        print(f"  - forest_mask.tif    (森林类型: {FOREST_CLASSES})")
        print(f"  - combined_mask.tif  (纬度 AND 森林)")
    else:
        print(f"  - forest_mask.tif    (占位文件，全1)")
        print(f"  - combined_mask.tif  (仅纬度掩膜，用于后续分析)")
    print("="*70)

# ==================== 运行 ====================
if __name__ == "__main__":
    main()
