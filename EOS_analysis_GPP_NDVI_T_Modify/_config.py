#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
配置文件: Wang (2025) 分析流程
NDVI物候 + GLEAM蒸腾版本
"""

from pathlib import Path
import os

# ==================== 数据路径配置 ====================

# 根目录（请修改为您的数据根目录）
ROOT = Path(r"I:\F\Data4")

# GLEAM数据路径
GLEAM_ROOT = ROOT / "Meteorological Data" / "GLEAM"

# TR (蒸腾) - 来自GLEAM Et（替代ERA5-Land）
TR_DAILY_DIR = GLEAM_ROOT / "Et" / "Et_Daily_1"

# 土壤水分路径（GLEAM） - 推荐使用深层SMrz
SM_DAILY_DIR = GLEAM_ROOT / "SMrz" / "SMrz_Daily_1"    # 深层土壤水分（日尺度）- 推荐
# SM_DAILY_DIR = GLEAM_ROOT / "SMs" / "SMs_Daily"      # 表层土壤水分（日尺度）- 备选

# NDVI数据路径（替代GPP）
NDVI_DAILY_DIR = ROOT / "GIMMS_NDVI" / "GIMMS_NDVI_daily_interpolated"  # NDVI日数据（已插值）

# GPP数据路径（重定向到NDVI，保持向后兼容）
GPP_DAILY_DIR = NDVI_DAILY_DIR  # 重定向到NDVI日数据（02-06模块使用）
GPP_8DAY_DIR = ROOT / "GLASS_GPP" / "GLASS_GPP_8days_1"  # 8天数据保留原路径（未使用）

# SIF数据路径（如需要可配置）
SIF_DAILY_DIR = ROOT / "SIF_Data" / "CSIF_daily"       # 日尺度SIF（暂未使用）
SIF_ANNUAL_DIR = ROOT / "SIF_Data" / "CSIF_annual"     # 年度SIF总量（暂未使用）

# 物候数据路径 - 使用NDVI物候
PHENO_DIR = ROOT / "Phenology_Output_1" / "NDVI_phenology"

# 土地覆盖数据
LANDCOVER_FILE = ROOT / "Landcover" / "MCD12Q1" / "MCD12Q1_IGBP_2018.tif"

# 输出目录（NDVI+GLEAM版本）
OUTPUT_ROOT = ROOT / "Wang2025_Analysis_EOS_GPP_NDVI_GLEAM"
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

# TR (蒸腾) 文件命名格式 - GLEAM格式
# 格式: Et_YYYYMMDD.tif
# 示例: Et_19820101.tif
TR_FILE_FORMAT = "Et_{date}.tif"  # {date} = YYYYMMDD

# 土壤水分文件命名格式 - GLEAM格式
# 格式: SMrz_YYYYMMDD.tif (日尺度)
SM_FILE_FORMAT = "SMrz_{date}.tif"  # 日尺度，{date} = YYYYMMDD

# NDVI文件命名格式
# 日尺度（插值后）: NDVI_YYYYMMDD.tif
NDVI_DAILY_FORMAT = "NDVI_{date}.tif"  # {date} = YYYYMMDD, 如: NDVI_19820101.tif

# GPP文件命名格式（重定向到NDVI格式，保持向后兼容）
GPP_DAILY_FORMAT = NDVI_DAILY_FORMAT  # 重定向到NDVI格式（02-06模块使用）
GPP_8DAY_FORMAT = "GLASS_GPP_{year}{doy:03d}.tif"  # 8天格式保留（未使用）

# 物候文件命名格式 - NDVI物候
PHENO_FILE_FORMAT = {
    'SOS': 'sos_ndvi_{year}.tif',      # NDVI物候
    'POS': 'pos_doy_ndvi_{year}.tif',  # POS使用pos_doy
    'EOS': 'eos_ndvi_{year}.tif',
    'POS_VALUE': 'pos_value_ndvi_{year}.tif',  # 峰值数值
    'QUALITY': 'quality_flags_ndvi_{year}.tif'  # 质量标记
}

# 年度数据命名格式（暂未使用）
ANNUAL_FILE_FORMAT = {
    'SIF': 'SIF_annual_{year}.tif',
    'SM': 'SMrz_{year}.tif',
    'TRc': 'TRc_{year}.tif'
}

# ==================== 输出文件命名格式（集中配置） ====================
# 修改此处可全局更改输出文件名中的"GPP"为"NDVI"

# 中间变量名称（用于输出文件名）
MIDDLE_VAR_NAME = "NDVI"  # 可选: "GPP" 或 "NDVI"

# ==================== 通用数据路径（根据MIDDLE_VAR_NAME自动选择） ====================
# 07模块等需要灵活切换数据源的模块使用这些通用配置
# 注意：这里不使用 GPP_DAILY_DIR（它被重定向到NDVI供其他模块使用）
if MIDDLE_VAR_NAME == "NDVI":
    DAILY_DATA_DIR = NDVI_DAILY_DIR
    DAILY_DATA_FORMAT = NDVI_DAILY_FORMAT
    PHENO_SOURCE_DIR = ROOT / "Phenology_Output_1" / "NDVI_phenology"
    PHENO_SOURCE_FORMAT = {
        'SOS': 'sos_ndvi_{year}.tif',
        'POS': 'pos_doy_ndvi_{year}.tif',
        'EOS': 'eos_ndvi_{year}.tif',
        'POS_VALUE': 'pos_value_ndvi_{year}.tif',
        'QUALITY': 'quality_flags_ndvi_{year}.tif'
    }
elif MIDDLE_VAR_NAME == "GPP":
    # 直接指定真实的GPP路径（不使用重定向的GPP_DAILY_DIR）
    DAILY_DATA_DIR = ROOT / "GLASS_GPP" / "GLASS_GPP_daily_interpolated"
    DAILY_DATA_FORMAT = "GPP_{date}.tif"
    PHENO_SOURCE_DIR = ROOT / "Phenology_Output_1" / "GPP_phenology"
    PHENO_SOURCE_FORMAT = {
        'SOS': 'sos_gpp_{year}.tif',
        'POS': 'pos_doy_gpp_{year}.tif',
        'EOS': 'eos_gpp_{year}.tif',
        'POS_VALUE': 'pos_value_gpp_{year}.tif',
        'QUALITY': 'quality_flags_gpp_{year}.tif'
    }
else:
    # 默认使用NDVI
    DAILY_DATA_DIR = NDVI_DAILY_DIR
    DAILY_DATA_FORMAT = NDVI_DAILY_FORMAT
    PHENO_SOURCE_DIR = PHENO_DIR
    PHENO_SOURCE_FORMAT = PHENO_FILE_FORMAT

# 02模块输出（累积值）
OUTPUT_CUMULATIVE_FORMAT = f"{MIDDLE_VAR_NAME}c_{{year}}.tif"  # NDVIc_2000.tif

# 02模块气候态
OUTPUT_CLIMATOLOGY_FORMAT = {
    'daily': f"{MIDDLE_VAR_NAME}_daily_climatology.tif",  # NDVI_daily_climatology.tif
    'cumulative_av': f"{MIDDLE_VAR_NAME}c_av.tif",        # NDVIc_av.tif
}

# 03c模块分解输出
OUTPUT_DECOMP_FORMAT = {
    'window_change': f"{MIDDLE_VAR_NAME}_window_change_{{year}}.tif",
    'fixed_window': f"{MIDDLE_VAR_NAME}_fixed_window_{{year}}.tif",
    'eos_change': f"{MIDDLE_VAR_NAME}_eos_change_{{year}}.tif",
    'pos_change': f"{MIDDLE_VAR_NAME}_pos_change_{{year}}.tif",
    'fixed_rate': f"Fixed_{MIDDLE_VAR_NAME}rate_{{year}}.tif",
    'fixed_window_length': f"Fixed_{MIDDLE_VAR_NAME}Window_Length.tif",
}

# 04c模块缓存文件
OUTPUT_CACHE_FORMAT = {
    'season': f"{MIDDLE_VAR_NAME}_{{season}}_{{year}}.tif",
    'lsp': f"{MIDDLE_VAR_NAME}_LSP_{{tag}}_{{year}}.tif",
}

# ==================== 绘图参数配置 ====================

# 图形尺寸
FIG_DPI = 300
FIG_FORMAT = 'png'  # 或 'pdf', 'svg'

# 绘图输出目录（06/07模块共用）
FIGURES_DIR = OUTPUT_ROOT / "Figures_All"

# 07模块输出文件名格式（物候异常分组生长曲线）
# 使用MIDDLE_VAR_NAME实现变量替换
OUTPUT_COMPOSITE_CURVES_FORMAT = {
    'english': f"phenology_composite_{MIDDLE_VAR_NAME}_curves_EN.png",
    'chinese': f"phenology_composite_{MIDDLE_VAR_NAME}_curves_CN.png",
    'decomposition': f"phenology_composite_{MIDDLE_VAR_NAME}_curves_with_decomp.png",
}

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
    根据日期获取TR文件路径（GLEAM格式）

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


def get_NDVI_file_path(date_obj):
    """
    获取NDVI文件路径（替代GPP）

    Parameters:
    -----------
    date_obj : datetime
        日期对象
    """
    yyyymmdd = date_obj.strftime("%Y%m%d")
    file_path = NDVI_DAILY_DIR / NDVI_DAILY_FORMAT.format(date=yyyymmdd)
    if file_path.exists():
        return file_path
    return None


# 保留GPP函数以兼容旧代码，但实际使用NDVI
def get_GPP_file_path(date_obj, daily=True):
    """
    获取GPP/NDVI文件路径（此版本返回NDVI路径）

    Parameters:
    -----------
    date_obj : datetime
        日期对象
    daily : bool
        True为日尺度，False为8天尺度
    """
    if daily:
        # 返回NDVI路径
        return get_NDVI_file_path(date_obj)
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

    if not NDVI_DAILY_DIR.exists():
        errors.append(f"NDVI数据目录不存在: {NDVI_DAILY_DIR}")

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
    print("当前配置信息 (NDVI物候 + GLEAM蒸腾版本)")
    print("="*70)
    print(f"\n数据路径:")
    print(f"  根目录: {ROOT}")
    print(f"  TR数据(GLEAM): {TR_DAILY_DIR}")
    print(f"  NDVI数据: {NDVI_DAILY_DIR}")
    print(f"  物候数据(NDVI): {PHENO_DIR}")
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
