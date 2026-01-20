#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
配置文件: Wang (2025) EOS分析流程
请根据您的数据路径修改以下配置
"""

from pathlib import Path
import os

# ==================== 数据路径配置 ====================

# 根目录（请修改为您的数据根目录）
ROOT = Path(r"I:\F\Data4")

# GLEAM数据路径
GLEAM_ROOT = ROOT / "Meteorological Data" / "GLEAM"
# TR (蒸腾) - 来自ERA5-Land（根据物候代码）
TR_DAILY_DIR = ROOT / "Meteorological Data" / "ERA5_Land" / "ET_components" / "ET_transp" / "ET_transp_Daily" / "ET_transp_Daily_2"

# 土壤水分路径（GLEAM） - 推荐使用深层SMrz
SM_DAILY_DIR = GLEAM_ROOT / "SMrz" / "SMrz_Daily_1"    # 深层土壤水分（日尺度）- 推荐
# SM_DAILY_DIR = GLEAM_ROOT / "SMs" / "SMs_Daily"      # 表层土壤水分（日尺度）- 备选

# GPP数据路径
GPP_DAILY_DIR = ROOT / "GLASS_GPP" / "GLASS_GPP_daily_interpolated"  # GPP日数据（已插值）
GPP_8DAY_DIR = ROOT / "GLASS_GPP" / "GLASS_GPP_8days_1"              # GPP 8天数据（原始）

# SIF数据路径（如需要可配置）
SIF_DAILY_DIR = ROOT / "SIF_Data" / "CSIF_daily"       # 日尺度SIF（暂未使用）
SIF_ANNUAL_DIR = ROOT / "SIF_Data" / "CSIF_annual"     # 年度SIF总量（暂未使用）

# 物候数据路径（从物候代码输出）
# 选项1: 使用GPP物候（原始目录；若已重投影则改为 *_EPSG4326）
PHENO_DIR = ROOT / "Phenology_Output_1" / "GPP_phenology"
# 已重投影到EPSG:4326的目录（如需切换，启用这一行）
# PHENO_DIR = ROOT / "Phenology_Output_1" / "GPP_phenology_EPSG4326"
# 选项2: 使用T物候（对比分析）
# PHENO_DIR = ROOT / "Phenology_Output_1" / "T_phenology"
# 选项3: 使用SIF物候（需单独提取）
# PHENO_DIR = ROOT / "Phenology_Output_1" / "SIF_phenology"

# 土地覆盖数据
LANDCOVER_FILE = ROOT / "Landcover" / "MCD12Q1" / "MCD12Q1_IGBP_2018.tif"

# 输出目录（EOS分析+GPP修改版本，避免覆盖原结果）
OUTPUT_ROOT = ROOT / "Wang2025_Analysis_EOS_GPP_Modify"
TRC_ANNUAL_DIR = OUTPUT_ROOT / "TRc_annual"              # TRc年度输出（02代码）
CLIMATOLOGY_DIR = OUTPUT_ROOT / "Climatology"            # 气候态数据（02代码）
DECOMPOSITION_DIR = OUTPUT_ROOT / "Decomposition"        # 分解结果（03a代码）
DECOMPOSITION_TIMING_DIR = OUTPUT_ROOT / "Decomposition_TimingShape"  # 时序/形状分解（03b代码）
DECOMPOSITION_FIXED_DIR = OUTPUT_ROOT / "Decomposition_FixedWindow"   # 固定窗口分解（03c代码）
STATISTICAL_DIR = OUTPUT_ROOT / "Statistical_Analysis"   # 统计分析（04a代码）
STATISTICAL_TIMING_DIR = OUTPUT_ROOT / "Statistical_Analysis_TimingShape"   # 统计分析（04b代码）
STATISTICAL_FIXED_DIR = OUTPUT_ROOT / "Statistical_Analysis_FixedWindow"     # 统计分析（04c代码）

# 统一模板与掩膜（全流程只认 TEMPLATE_RASTER）
TEMPLATE_RASTER = OUTPUT_ROOT / "masks" / "template_grid.tif"
MASK_FILE = OUTPUT_ROOT / "masks" / "combined_mask.tif"

# ==================== 分析参数配置 ====================

# 时间范围
YEAR_START = 1982
YEAR_END = 2018

# 纬度范围（≥30°N）
LAT_MIN = 30.0

# ==================== 掩膜策略配置 ====================
# 重要：控制是否使用森林掩膜
USE_FOREST_MASK = False  # True: 仅森林, False: 所有陆地像元（≥30°N）

# 森林类型（IGBP分类）- 仅当 USE_FOREST_MASK=True 时生效
FOREST_CLASSES = [1, 2, 3, 4, 5]
# 1: Evergreen Needleleaf Forest (ENF)
# 2: Evergreen Broadleaf Forest (EBF)
# 3: Deciduous Needleleaf Forest (DNF)
# 4: Deciduous Broadleaf Forest (DBF)
# 5: Mixed Forest (MF)

# ==================== 处理参数配置 ====================

# 块处理参数
BLOCK_SIZE = 128          # 块大小（像素）
MAX_WORKERS = 2           # 最大并行进程数

# 数据质量控制
NODATA_OUT = -9999.0     # 输出NODATA值
MIN_VALID_FRAC = 0.60    # 最小有效数据比例（趋势分析）

# 物候提取参数
PHENO_THRESHOLD = 0.20   # 物候提取阈值（20%）
EOS_MIN_DOY = 200        # EOS最小DOY（约束条件，经验值）

# 统计分析参数
TREND_PVALUE_THRESHOLD = 0.05      # 趋势显著性阈值
MOVING_WINDOW_SIZE = 15            # 滑动窗口大小（年）

# ==================== 文件命名格式配置 ====================

# TR (蒸腾) 文件命名格式 - ERA5-Land格式
# 格式: ERA5L_ET_transp_Daily_mm_YYYYMMDD.tif
# 示例: ERA5L_ET_transp_Daily_mm_19820101.tif
TR_FILE_FORMAT = "ERA5L_ET_transp_Daily_mm_{date}.tif"  # {date} = YYYYMMDD

# 土壤水分文件命名格式 - GLEAM格式
# 格式: SMrz_YYYYMMDD.tif (日尺度)
# 根据用户提供: SMrz_19801218.tif (YYYYMMDD格式，日尺度)
SM_FILE_FORMAT = "SMrz_{date}.tif"  # 日尺度，{date} = YYYYMMDD

# GPP文件命名格式
# 日尺度（插值后）: GPP_YYYYMMDD.tif
# 8天数据: GLASS_GPP_YYYYDDD.tif
GPP_DAILY_FORMAT = "GPP_{date}.tif"  # {date} = YYYYMMDD, 如: GPP_19820101.tif
GPP_8DAY_FORMAT = "GLASS_GPP_{year}{doy:03d}.tif"  # 如: GLASS_GPP_2000001.tif

# 物候文件命名格式（物候代码输出）
PHENO_FILE_FORMAT = {
    'SOS': 'sos_gpp_{year}.tif',      # 注意：物候代码输出小写
    'POS': 'pos_doy_gpp_{year}.tif',  # POS使用pos_doy
    'EOS': 'eos_gpp_{year}.tif',
    'POS_VALUE': 'pos_value_gpp_{year}.tif',  # 峰值数值
    'QUALITY': 'quality_flags_gpp_{year}.tif'  # 质量标记
}

# 年度数据命名格式（暂未使用）
ANNUAL_FILE_FORMAT = {
    'SIF': 'SIF_annual_{year}.tif',
    'SM': 'SMrz_{year}.tif',
    'TRc': 'TRc_{year}.tif'
}

# ==================== 绘图参数配置 ====================

# 图形尺寸
FIG_DPI = 300
FIG_FORMAT = 'png'  # 或 'pdf', 'svg'

# 地图范围（经度，纬度）
MAP_EXTENT = [-180, 180, 30, 90]  # [lon_min, lon_max, lat_min, lat_max]

# 色标范围（可根据数据调整）
CBAR_RANGES = {
    'TRc': [0, 500],              # mm
    'LSP': [60, 150],             # days
    'EOS': [200, 300],            # DOY
    'trend': [-20, 20],           # mm/decade
    'attribution': [-0.5, 0.5]    # standardized coefficient
}

# ==================== 日志配置 ====================

# 日志级别: DEBUG, INFO, WARNING, ERROR
LOG_LEVEL = "INFO"

# 日志格式
LOG_FORMAT = '%(asctime)s [%(levelname)s] %(message)s'

# ==================== 高级配置（一般无需修改） ====================

# 内存优化
USE_MEMORY_MAPPING = True         # 大文件使用内存映射
CHUNK_SIZE = 1000                 # 逐像元处理时的批次大小

# 数值精度
FLOAT_PRECISION = 'float32'       # 输出数据精度

# 并行处理
USE_MULTIPROCESSING = True        # 是否使用多进程
PREFETCH_SIZE = 2                 # 预读取队列大小

# ==================== 函数：获取完整路径 ====================

def get_TR_file_path(date_obj):
    """
    根据日期获取TR文件路径（ERA5-Land格式）

    Parameters:
    -----------
    date_obj : datetime
        日期对象

    Returns:
    --------
    file_path : Path
        文件完整路径，如果文件不存在返回None
    """
    yyyymmdd = date_obj.strftime("%Y%m%d")

    # 使用TR_FILE_FORMAT
    file_path = TR_DAILY_DIR / TR_FILE_FORMAT.format(date=yyyymmdd)

    if file_path.exists():
        return file_path

    return None

def get_SM_file_path(date_obj):
    """
    根据日期获取土壤水分文件路径（日尺度）

    Parameters:
    -----------
    date_obj : datetime
        日期对象

    Returns:
    --------
    file_path : Path
        文件完整路径，如果文件不存在返回None
    """
    yyyymmdd = date_obj.strftime("%Y%m%d")
    file_path = SM_DAILY_DIR / SM_FILE_FORMAT.format(date=yyyymmdd)
    if file_path.exists():
        return file_path
    return None

def _normalize_pheno_key(var_name):
    key = str(var_name).strip().upper()
    aliases = {
        "SOS": "SOS",
        "POS": "POS",
        "POS_DOY": "POS",
        "POSDOY": "POS",
        "EOS": "EOS",
        "POS_VALUE": "POS_VALUE",
        "POSVALUE": "POS_VALUE",
        "QUALITY": "QUALITY",
        "QUALITY_FLAGS": "QUALITY",
        "QUALITYFLAGS": "QUALITY",
    }
    return aliases.get(key, key)


def get_pheno_file_path(var_name, year, check_exists=True):
    """获取物候文件路径"""
    key = _normalize_pheno_key(var_name)
    if key not in PHENO_FILE_FORMAT:
        raise KeyError(f"Unknown PHENO_FILE_FORMAT key: {var_name}")
    filename = PHENO_FILE_FORMAT[key].format(year=year)
    file_path = PHENO_DIR / filename
    if check_exists and not file_path.exists():
        print(f"  ⚠ Missing pheno file: {file_path}")
        return None
    return file_path


def get_GPP_file_path(date_obj, daily=True):
    """
    获取GPP文件路径

    Parameters:
    -----------
    date_obj : datetime
        日期对象
    daily : bool
        True为日尺度，False为8天尺度
    """
    if daily:
        yyyymmdd = date_obj.strftime("%Y%m%d")
        file_path = GPP_DAILY_DIR / GPP_DAILY_FORMAT.format(date=yyyymmdd)
    else:
        year = date_obj.year
        doy = date_obj.timetuple().tm_yday
        file_path = GPP_8DAY_DIR / GPP_8DAY_FORMAT.format(year=year, doy=doy)
    if file_path.exists():
        return file_path
    return None

def get_annual_file_path(var_name, year):
    """获取年度数据文件路径"""
    if var_name == 'SIF':
        return SIF_ANNUAL_DIR / f"SIF_annual_{year}.tif"
    elif var_name == 'SM':
        # 年度SM文件位于SM_DAILY_DIR的父目录
        return GLEAM_ROOT / "SMrz" / f"SMrz_{year}.tif"
    elif var_name == 'TRc':
        return OUTPUT_ROOT / "TRc_annual" / f"TRc_{year}.tif"
    else:
        raise ValueError(f"Unknown variable: {var_name}")

# ==================== 配置验证 ====================

def validate_config():
    """验证配置的有效性"""
    errors = []

    # 检查关键路径是否存在
    if not ROOT.exists():
        errors.append(f"根目录不存在: {ROOT}")

    if not TR_DAILY_DIR.exists():
        errors.append(f"TR数据目录不存在: {TR_DAILY_DIR}")

    if not GPP_DAILY_DIR.exists():
        errors.append(f"GPP数据目录不存在: {GPP_DAILY_DIR}")

    if not SM_DAILY_DIR.exists():
        errors.append(f"土壤水分目录不存在: {SM_DAILY_DIR}")

    if not PHENO_DIR.exists():
        errors.append(f"物候数据目录不存在: {PHENO_DIR}")

    if USE_FOREST_MASK and not LANDCOVER_FILE.exists():
        errors.append(f"土地覆盖文件不存在: {LANDCOVER_FILE}")

    output_parent = OUTPUT_ROOT.parent
    if OUTPUT_ROOT.exists():
        if not os.access(OUTPUT_ROOT, os.W_OK):
            errors.append(f"输出目录不可写: {OUTPUT_ROOT}")
    else:
        if not os.access(output_parent, os.W_OK):
            errors.append(f"输出目录父路径不可写: {output_parent}")

    # 检查参数合理性
    if YEAR_START >= YEAR_END:
        errors.append(f"年份范围错误: {YEAR_START} >= {YEAR_END}")

    if not 0 < MIN_VALID_FRAC <= 1:
        errors.append(f"有效数据比例错误: {MIN_VALID_FRAC}")

    if BLOCK_SIZE <= 0 or BLOCK_SIZE > 2048:
        errors.append(f"块大小错误: {BLOCK_SIZE}")

    if errors:
        print("\n配置错误:")
        for error in errors:
            print(f"  ✗ {error}")
        return False
    else:
        print("\n✓ 配置验证通过")
        return True

# ==================== 配置信息打印 ====================

def print_config():
    """打印当前配置"""
    print("\n" + "="*70)
    print("当前配置信息")
    print("="*70)
    print(f"\n数据路径:")
    print(f"  根目录: {ROOT}")
    print(f"  TR数据: {TR_DAILY_DIR}")
    print(f"  物候数据: {PHENO_DIR}")
    print(f"  输出目录: {OUTPUT_ROOT}")

    print(f"\n时间范围:")
    print(f"  起始年: {YEAR_START}")
    print(f"  结束年: {YEAR_END}")
    print(f"  总年数: {YEAR_END - YEAR_START + 1}")

    print(f"\n处理参数:")
    print(f"  块大小: {BLOCK_SIZE}")
    print(f"  并行进程: {MAX_WORKERS}")
    print(f"  最小有效率: {MIN_VALID_FRAC*100:.0f}%")

    print("\n" + "="*70)

if __name__ == "__main__":
    print_config()
    validate_config()
