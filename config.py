#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
配置文件: Wang (2025) 分析流程
请根据您的数据路径修改以下配置
"""

from pathlib import Path

# ==================== 数据路径配置 ====================

# 根目录（请修改为您的数据根目录）
ROOT = Path(r"I:\F\Data4")

# GLEAM数据路径
GLEAM_ROOT = ROOT / "Meteorological Data" / "GLEAM"
# TR (蒸腾) - 来自ERA5-Land（根据物候代码）
TR_DAILY_DIR = ROOT / "Meteorological Data" / "ERA5_Land" / "ET_components" / "ET_transp" / "ET_transp_Daily" / "ET_transp_Daily_2"

# 土壤水分路径（GLEAM） - 推荐使用深层SMrz
SM_DAILY_DIR = GLEAM_ROOT / "SMrz" / "SMrz_Daily"      # 深层土壤水分（日尺度）- 推荐
# SM_DAILY_DIR = GLEAM_ROOT / "SMs" / "SMs_Daily"      # 表层土壤水分（日尺度）- 备选

# GPP数据路径
GPP_DAILY_DIR = ROOT / "GLASS_GPP" / "GLASS_GPP_daily_interpolated"  # GPP日数据（已插值）
GPP_8DAY_DIR = ROOT / "GLASS_GPP" / "GLASS_GPP_8days_1"              # GPP 8天数据（原始）

# SIF数据路径（如需要可配置）
SIF_DAILY_DIR = ROOT / "SIF_Data" / "CSIF_daily"       # 日尺度SIF（暂未使用）
SIF_ANNUAL_DIR = ROOT / "SIF_Data" / "CSIF_annual"     # 年度SIF总量（暂未使用）

# 物候数据路径（从物候代码输出）
# 选项1: 使用GPP物候（推荐，已计算）
PHENO_DIR = ROOT / "Phenology_Output_1" / "GPP_phenology"
# 选项2: 使用T物候（对比分析）
# PHENO_DIR = ROOT / "Phenology_Output_1" / "T_phenology"
# 选项3: 使用SIF物候（需单独提取）
# PHENO_DIR = ROOT / "Phenology_Output_1" / "SIF_phenology"

# 土地覆盖数据
LANDCOVER_FILE = ROOT / "Landcover" / "MCD12Q1" / "MCD12Q1_IGBP_2018.tif"

# 输出目录
OUTPUT_ROOT = ROOT / "Wang2025_Analysis"

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
SOS_MAX_DOY = 150        # SOS最大DOY（约束条件）

# 统计分析参数
TREND_PVALUE_THRESHOLD = 0.05      # 趋势显著性阈值
MOVING_WINDOW_SIZE = 15            # 滑动窗口大小（年）

# ==================== 文件命名格式配置 ====================

# TR (蒸腾) 文件命名格式 - ERA5-Land格式
# 格式: ERA5L_ET_transp_Daily_mm_YYYYMMDD.tif
# 示例: ERA5L_ET_transp_Daily_mm_19820101.tif
TR_FILE_PREFIX = "ERA5L_ET_transp_Daily_mm_"
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
    'SOS': [60, 150],             # DOY
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

    # ERA5-Land格式: ERA5L_ET_transp_Daily_mm_YYYYMMDD.tif
    file_path = TR_DAILY_DIR / f"ERA5L_ET_transp_Daily_mm_{yyyymmdd}.tif"

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

    # 日尺度格式: SMrz_YYYYMMDD.tif
    file_path = SM_DAILY_DIR / f"SMrz_{yyyymmdd}.tif"

    if file_path.exists():
        return file_path

    return None

def get_pheno_file_path(var_name, year):
    """获取物候文件路径"""
    filename = PHENO_FILE_FORMAT[var_name].format(year=year)
    return PHENO_DIR / filename

def get_annual_file_path(var_name, year):
    """获取年度数据文件路径"""
    if var_name == 'SIF':
        return SIF_ANNUAL_DIR / f"SIF_annual_{year}.tif"
    elif var_name == 'SM':
        return SM_ROOT_DIR / f"SMrz_{year}.tif"
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

    if not PHENO_DIR.exists():
        errors.append(f"物候数据目录不存在: {PHENO_DIR}")

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
