#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
配置文件: Wang (2025) 分析流程
NDVI物候 + GLEAM蒸腾版本 — EOS分析
"""

from pathlib import Path
import os

# ==================== 数据路径配置 ====================

# 根目录（请修改为您的数据根目录）
ROOT = Path(r"I:\F\Data4")

# GLEAM数据路径
GLEAM_ROOT = ROOT / "Meteorological Data" / "GLEAM"

# TR (蒸腾) 数据源 — 固定常量
GLEAM_TR_DAILY_DIR = GLEAM_ROOT / "Et" / "Et_Daily_1"
ERA5_TR_DAILY_DIR = ROOT / "Meteorological Data" / "ERA5_Land" / "ET_components" / "ET_transp" / "ET_transp_Daily" / "ET_transp_Daily_2"

# TR (蒸腾) - 当前使用（默认GLEAM，运行时由apply_run_config切换）
TR_DAILY_DIR = GLEAM_TR_DAILY_DIR

# 土壤水分路径（GLEAM） - 推荐使用深层SMrz
SM_DAILY_DIR = GLEAM_ROOT / "SMrz" / "SMrz_Daily_1"    # 深层土壤水分（日尺度）- 推荐
# SM_DAILY_DIR = GLEAM_ROOT / "SMs" / "SMs_Daily"      # 表层土壤水分（日尺度）- 备选

# NDVI数据路径（GIMMS NDVI3g）
NDVI_DAILY_DIR = ROOT / "GIMMS_NDVI" / "GIMMS_NDVI_daily_interpolated"

# GPP数据路径（GLASS GPP）— 固定常量，始终指向真实GPP数据
GPP_DAILY_DIR = ROOT / "GLASS_GPP" / "GLASS_GPP_daily_interpolated"
GPP_8DAY_DIR = ROOT / "GLASS_GPP" / "GLASS_GPP_8days_1"

# GPP物候数据路径
GPP_PHENO_DIR = ROOT / "Phenology_Output_1" / "GPP_phenology"

# SIF数据路径（如需要可配置）
SIF_DAILY_DIR = ROOT / "SIF_Data" / "CSIF_daily"       # 日尺度SIF（暂未使用）
SIF_ANNUAL_DIR = ROOT / "SIF_Data" / "CSIF_annual"     # 年度SIF总量（暂未使用）

# NDVI物候数据路径（固定常量）
NDVI_PHENO_DIR = ROOT / "Phenology_Output_1" / "NDVI_phenology"

# 土地覆盖数据
LANDCOVER_FILE = ROOT / "Landcover" / "MCD12Q1" / "MCD12Q1_IGBP_2018.tif"

# 输出目录（EOS分析版本）
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

# TR (蒸腾) 文件命名格式 — 固定常量
GLEAM_TR_FILE_FORMAT = "Et_{date}.tif"  # GLEAM格式: Et_YYYYMMDD.tif
ERA5_TR_FILE_FORMAT = "ERA5L_ET_transp_Daily_mm_{date}.tif"  # ERA5-Land格式

# TR文件命名格式 - 当前使用（默认GLEAM，运行时由apply_run_config切换）
TR_FILE_FORMAT = GLEAM_TR_FILE_FORMAT

# 土壤水分文件命名格式 - GLEAM格式
# 格式: SMrz_YYYYMMDD.tif (日尺度)
SM_FILE_FORMAT = "SMrz_{date}.tif"  # 日尺度，{date} = YYYYMMDD

# NDVI文件命名格式
# 日尺度（插值后）: NDVI_YYYYMMDD.tif
NDVI_DAILY_FORMAT = "NDVI_{date}.tif"  # {date} = YYYYMMDD, 如: NDVI_19820101.tif

# GPP文件命名格式（GLASS GPP）— 固定常量
GPP_DAILY_FORMAT = "GPP_{date}.tif"  # {date} = YYYYMMDD, 如: GPP_19820101.tif
GPP_8DAY_FORMAT = "GLASS_GPP_{year}{doy:03d}.tif"

# GPP物候文件命名格式
GPP_PHENO_FORMAT = {
    'SOS': 'sos_gpp_{year}.tif',
    'POS': 'pos_doy_gpp_{year}.tif',
    'EOS': 'eos_gpp_{year}.tif',
    'POS_VALUE': 'pos_value_gpp_{year}.tif',
    'QUALITY': 'quality_flags_gpp_{year}.tif'
}

# NDVI物候文件命名格式（固定常量）
NDVI_PHENO_FORMAT = {
    'SOS': 'sos_ndvi_{year}.tif',
    'POS': 'pos_doy_ndvi_{year}.tif',
    'EOS': 'eos_ndvi_{year}.tif',
    'POS_VALUE': 'pos_value_ndvi_{year}.tif',
    'QUALITY': 'quality_flags_ndvi_{year}.tif'
}

# 年度数据命名格式（暂未使用）
ANNUAL_FILE_FORMAT = {
    'SIF': 'SIF_annual_{year}.tif',
    'SM': 'SMrz_{year}.tif',
    'TRc': 'TRc_{year}.tif'
}

# ==================== 数据模式配置 ====================
# 控制02-06模块处理的数据类型:
#   "NDVI" - 仅处理NDVI数据
#   "GPP"  - 仅处理GPP数据
#   "BOTH" - 依次处理NDVI和GPP（自动遍历）
DATA_MODE = "BOTH"

# 控制TR数据源:
#   "GLEAM" - 仅使用GLEAM蒸腾数据
#   "ERA5"  - 仅使用ERA5-Land蒸腾数据
#   "BOTH"  - 依次使用两种TR数据源（自动遍历）
TR_MODE = "BOTH"

# ==================== 可切换别名（02-06模块导入使用） ====================
# 默认值为NDVI模式（向后兼容）
# 运行时通过 apply_run_config() 切换到 GPP 或回到 NDVI
MIDDLE_VAR_NAME = "NDVI"

# 当前模式的日尺度数据路径（VAR = NDVI 或 GPP，随模式切换）
VAR_DAILY_DIR = NDVI_DAILY_DIR
VAR_DAILY_FORMAT = NDVI_DAILY_FORMAT

# 当前模式的物候路径
PHENO_DIR = NDVI_PHENO_DIR
PHENO_FILE_FORMAT = NDVI_PHENO_FORMAT.copy()

# ==================== 通用数据路径（07模块兼容） ====================
DAILY_DATA_DIR = NDVI_DAILY_DIR
DAILY_DATA_FORMAT = NDVI_DAILY_FORMAT
PHENO_SOURCE_DIR = NDVI_PHENO_DIR
PHENO_SOURCE_FORMAT = NDVI_PHENO_FORMAT

# ==================== 07模块多组合运行配置 ====================
# 五种数据源组合：物候数据来源 × 日尺度数据来源
RUN_CONFIGS_07 = [
    {
        'run_name': 'NDVI_pheno_NDVI_daily',       # 输出子目录名
        'var_label': 'NDVI',                        # 图表标签
        'pheno_dir': NDVI_PHENO_DIR,
        'pheno_format': NDVI_PHENO_FORMAT,
        'daily_dir': NDVI_DAILY_DIR,
        'daily_format': NDVI_DAILY_FORMAT,
        'data_type': 'ndvi',
    },
    {
        'run_name': 'GPP_pheno_GPP_daily',
        'var_label': 'GPP',
        'pheno_dir': GPP_PHENO_DIR,
        'pheno_format': GPP_PHENO_FORMAT,
        'daily_dir': GPP_DAILY_DIR,
        'daily_format': GPP_DAILY_FORMAT,
        'data_type': 'gpp',
    },
    {
        'run_name': 'NDVI_pheno_GPP_daily',
        'var_label': 'GPP(NDVI_pheno)',
        'pheno_dir': NDVI_PHENO_DIR,
        'pheno_format': NDVI_PHENO_FORMAT,
        'daily_dir': GPP_DAILY_DIR,
        'daily_format': GPP_DAILY_FORMAT,
        'data_type': 'gpp',
    },
    {
        'run_name': 'NDVI_pheno_GLEAM_TR_daily',
        'var_label': 'GLEAM_TR(NDVI_pheno)',
        'pheno_dir': NDVI_PHENO_DIR,
        'pheno_format': NDVI_PHENO_FORMAT,
        'daily_dir': GLEAM_TR_DAILY_DIR,
        'daily_format': GLEAM_TR_FILE_FORMAT,
        'data_type': 'tr',
    },
    {
        'run_name': 'NDVI_pheno_ERA5_TR_daily',
        'var_label': 'ERA5_TR(NDVI_pheno)',
        'pheno_dir': NDVI_PHENO_DIR,
        'pheno_format': NDVI_PHENO_FORMAT,
        'daily_dir': ERA5_TR_DAILY_DIR,
        'daily_format': ERA5_TR_FILE_FORMAT,
        'data_type': 'tr',
    },
    {
        'run_name': 'GPP_pheno_GLEAM_TR_daily',
        'var_label': 'GLEAM_TR(GPP_pheno)',
        'pheno_dir': GPP_PHENO_DIR,
        'pheno_format': GPP_PHENO_FORMAT,
        'daily_dir': GLEAM_TR_DAILY_DIR,
        'daily_format': GLEAM_TR_FILE_FORMAT,
        'data_type': 'tr',
    },
    {
        'run_name': 'GPP_pheno_ERA5_TR_daily',
        'var_label': 'ERA5_TR(GPP_pheno)',
        'pheno_dir': GPP_PHENO_DIR,
        'pheno_format': GPP_PHENO_FORMAT,
        'daily_dir': ERA5_TR_DAILY_DIR,
        'daily_format': ERA5_TR_FILE_FORMAT,
        'data_type': 'tr',
    },
]

# ==================== 02-06模块运行配置 ====================
# 根据 DATA_MODE × TR_MODE 构建运行列表（最多2×2=4种组合）

# 基础配置模板（植被数据 × TR数据源）
_NDVI_GLEAM_CONFIG = {
    'MIDDLE_VAR_NAME': 'NDVI',
    'VAR_DAILY_DIR': NDVI_DAILY_DIR,
    'VAR_DAILY_FORMAT': NDVI_DAILY_FORMAT,
    'PHENO_DIR': NDVI_PHENO_DIR,
    'PHENO_FILE_FORMAT': NDVI_PHENO_FORMAT,
    'TR_DAILY_DIR': GLEAM_TR_DAILY_DIR,
    'TR_FILE_FORMAT': GLEAM_TR_FILE_FORMAT,
    'TR_SOURCE_NAME': 'GLEAM',
}
_NDVI_ERA5_CONFIG = {
    'MIDDLE_VAR_NAME': 'NDVI',
    'VAR_DAILY_DIR': NDVI_DAILY_DIR,
    'VAR_DAILY_FORMAT': NDVI_DAILY_FORMAT,
    'PHENO_DIR': NDVI_PHENO_DIR,
    'PHENO_FILE_FORMAT': NDVI_PHENO_FORMAT,
    'TR_DAILY_DIR': ERA5_TR_DAILY_DIR,
    'TR_FILE_FORMAT': ERA5_TR_FILE_FORMAT,
    'TR_SOURCE_NAME': 'ERA5',
}
_GPP_GLEAM_CONFIG = {
    'MIDDLE_VAR_NAME': 'GPP',
    'VAR_DAILY_DIR': GPP_DAILY_DIR,
    'VAR_DAILY_FORMAT': GPP_DAILY_FORMAT,
    'PHENO_DIR': GPP_PHENO_DIR,
    'PHENO_FILE_FORMAT': GPP_PHENO_FORMAT,
    'TR_DAILY_DIR': GLEAM_TR_DAILY_DIR,
    'TR_FILE_FORMAT': GLEAM_TR_FILE_FORMAT,
    'TR_SOURCE_NAME': 'GLEAM',
}
_GPP_ERA5_CONFIG = {
    'MIDDLE_VAR_NAME': 'GPP',
    'VAR_DAILY_DIR': GPP_DAILY_DIR,
    'VAR_DAILY_FORMAT': GPP_DAILY_FORMAT,
    'PHENO_DIR': GPP_PHENO_DIR,
    'PHENO_FILE_FORMAT': GPP_PHENO_FORMAT,
    'TR_DAILY_DIR': ERA5_TR_DAILY_DIR,
    'TR_FILE_FORMAT': ERA5_TR_FILE_FORMAT,
    'TR_SOURCE_NAME': 'ERA5',
}

# 按 DATA_MODE × TR_MODE 构建运行列表
_all_configs = {
    ('NDVI', 'GLEAM'): _NDVI_GLEAM_CONFIG,
    ('NDVI', 'ERA5'):  _NDVI_ERA5_CONFIG,
    ('GPP',  'GLEAM'): _GPP_GLEAM_CONFIG,
    ('GPP',  'ERA5'):  _GPP_ERA5_CONFIG,
}

_data_modes = ['NDVI', 'GPP'] if DATA_MODE == "BOTH" else [DATA_MODE]
_tr_modes = ['GLEAM', 'ERA5'] if TR_MODE == "BOTH" else [TR_MODE]

RUN_CONFIGS_02_06 = [_all_configs[(d, t)] for d in _data_modes for t in _tr_modes]



def _build_output_formats(var_name):
    """根据变量名构建输出文件格式字典"""
    return {
        'OUTPUT_CUMULATIVE_FORMAT': f"{var_name}c_{{year}}.tif",
        'OUTPUT_CLIMATOLOGY_FORMAT': {
            'daily': f"{var_name}_daily_climatology.tif",
            'cumulative_av': f"{var_name}c_av.tif",
        },
        'OUTPUT_DECOMP_FORMAT': {
            'window_change': f"{var_name}_window_change_{{year}}.tif",
            'fixed_window': f"{var_name}_fixed_window_{{year}}.tif",
            'eos_change': f"{var_name}_eos_change_{{year}}.tif",
            'pos_change': f"{var_name}_pos_change_{{year}}.tif",
            'fixed_rate': f"Fixed_{var_name}rate_{{year}}.tif",
            'fixed_window_length': f"Fixed_{var_name}Window_Length.tif",
        },
        'OUTPUT_CACHE_FORMAT': {
            'season': f"{var_name}_{{season}}_{{year}}.tif",
            'lsp': f"{var_name}_LSP_{{tag}}_{{year}}.tif",
        },
    }


def apply_run_config(cfg, caller_globals=None):
    """
    将运行配置应用到调用模块的全局命名空间。

    用法（在02-06模块的 __main__ 中）:
        from _config import RUN_CONFIGS_02_06, apply_run_config
        for cfg in RUN_CONFIGS_02_06:
            apply_run_config(cfg, globals())
            # 重置局部别名（如 OUTPUT_DIR = TRC_ANNUAL_DIR）
            main()

    原理:
        from _config import VAR_DAILY_DIR 会在导入时创建局部绑定，
        后续修改 _config.VAR_DAILY_DIR 不会影响已导入的局部变量。
        此函数通过直接修改调用者的 globals() 字典来更新这些绑定。
    """
    import _config

    var_name = cfg['MIDDLE_VAR_NAME']

    tr_source = cfg.get('TR_SOURCE_NAME', 'GLEAM')

    # 模式相关输出根目录（包含TR源标识）
    mode_root = OUTPUT_ROOT / f"{var_name}_{tr_source}"

    # 数据路径 + 输出目录更新
    updates = {
        'MIDDLE_VAR_NAME': var_name,
        'VAR_DAILY_DIR': cfg['VAR_DAILY_DIR'],
        'VAR_DAILY_FORMAT': cfg['VAR_DAILY_FORMAT'],
        'PHENO_DIR': cfg['PHENO_DIR'],
        'PHENO_FILE_FORMAT': cfg['PHENO_FILE_FORMAT'],
        # TR数据源切换
        'TR_DAILY_DIR': cfg.get('TR_DAILY_DIR', GLEAM_TR_DAILY_DIR),
        'TR_FILE_FORMAT': cfg.get('TR_FILE_FORMAT', GLEAM_TR_FILE_FORMAT),
        # 模式相关输出目录（NDVI/ 或 GPP/ 子目录）
        'TRC_ANNUAL_DIR': mode_root / "TRc_annual",
        'CLIMATOLOGY_DIR': mode_root / "Climatology",
        'DECOMPOSITION_FIXED_DIR': mode_root / "Decomposition_FixedWindow",
        'STATISTICAL_FIXED_DIR': mode_root / "Statistical_Analysis_FixedWindow",
        'FIGURES_DIR': mode_root / "Figures_All",
    }
    # 输出文件格式更新
    updates.update(_build_output_formats(var_name))

    # 1) 更新 _config 模块自身（影响 get_var_daily_file_path 等函数）
    for key, value in updates.items():
        setattr(_config, key, value)

    # 2) 更新调用者模块的全局变量
    if caller_globals is not None:
        for key, value in updates.items():
            if key in caller_globals:
                caller_globals[key] = value

    # 3) 确保模式相关输出目录存在
    for dir_key in ('TRC_ANNUAL_DIR', 'CLIMATOLOGY_DIR', 'DECOMPOSITION_FIXED_DIR',
                    'STATISTICAL_FIXED_DIR', 'FIGURES_DIR'):
        updates[dir_key].mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"  数据模式切换: {var_name} + {tr_source}_TR")
    print(f"  日尺度数据: {cfg['VAR_DAILY_DIR']}")
    print(f"  TR数据:   {cfg.get('TR_DAILY_DIR', 'default')}")
    print(f"  物候数据:   {cfg['PHENO_DIR']}")
    print(f"  输出根目录: {mode_root}")
    print(f"{'='*60}")


# ==================== 输出文件命名格式（默认NDVI，运行时由apply_run_config更新） ====================
_default_formats = _build_output_formats(MIDDLE_VAR_NAME)
OUTPUT_CUMULATIVE_FORMAT = _default_formats['OUTPUT_CUMULATIVE_FORMAT']
OUTPUT_CLIMATOLOGY_FORMAT = _default_formats['OUTPUT_CLIMATOLOGY_FORMAT']
OUTPUT_DECOMP_FORMAT = _default_formats['OUTPUT_DECOMP_FORMAT']
OUTPUT_CACHE_FORMAT = _default_formats['OUTPUT_CACHE_FORMAT']

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
    根据日期获取TR文件路径（当前TR数据源）

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
    获取NDVI文件路径

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


def get_var_daily_file_path(date_obj, daily=True):
    """
    获取当前模式的日尺度数据文件路径（可切换的 VAR_DAILY_DIR）。
    NDVI模式下返回NDVI文件，GPP模式下返回GPP文件。

    Parameters:
    -----------
    date_obj : datetime
        日期对象
    """
    yyyymmdd = date_obj.strftime("%Y%m%d")
    file_path = VAR_DAILY_DIR / VAR_DAILY_FORMAT.format(date=yyyymmdd)
    if file_path.exists():
        return file_path
    return None


def get_GPP_file_path(date_obj, daily=True):
    """
    获取真实GPP数据文件路径（固定常量，始终指向GLASS GPP数据）。
    不随 DATA_MODE 切换。

    Parameters:
    -----------
    date_obj : datetime
        日期对象
    daily : bool
        True为日尺度GPP，False为8天GPP
    """
    if daily:
        yyyymmdd = date_obj.strftime("%Y%m%d")
        file_path = GPP_DAILY_DIR / GPP_DAILY_FORMAT.format(date=yyyymmdd)
        if file_path.exists():
            return file_path
        return None
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

    if not ERA5_TR_DAILY_DIR.exists():
        errors.append(f"ERA5 TR数据目录不存在: {ERA5_TR_DAILY_DIR}")

    if not NDVI_DAILY_DIR.exists():
        errors.append(f"NDVI数据目录不存在: {NDVI_DAILY_DIR}")

    if not GPP_DAILY_DIR.exists():
        errors.append(f"GPP数据目录不存在: {GPP_DAILY_DIR}")

    if not SM_DAILY_DIR.exists():
        errors.append(f"土壤水分目录不存在: {SM_DAILY_DIR}")

    if not NDVI_PHENO_DIR.exists():
        errors.append(f"NDVI物候目录不存在: {NDVI_PHENO_DIR}")

    if not GPP_PHENO_DIR.exists():
        errors.append(f"GPP物候目录不存在: {GPP_PHENO_DIR}")

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
    print("当前配置信息 (EOS分析 — NDVI物候 + 双TR数据源版本)")
    print("="*70)
    print(f"\n数据模式: DATA_MODE = {DATA_MODE}, TR_MODE = {TR_MODE}")
    run_labels = [f"{c['MIDDLE_VAR_NAME']}+{c.get('TR_SOURCE_NAME', 'GLEAM')}" for c in RUN_CONFIGS_02_06]
    print(f"  02-06将处理: {' → '.join(run_labels)}")
    print(f"\n数据路径:")
    print(f"  根目录: {ROOT}")
    print(f"  TR数据(GLEAM): {GLEAM_TR_DAILY_DIR}")
    print(f"  TR数据(ERA5):  {ERA5_TR_DAILY_DIR}")
    print(f"  NDVI数据(GIMMS): {NDVI_DAILY_DIR}")
    print(f"  GPP数据(GLASS): {GPP_DAILY_DIR}")
    print(f"  NDVI物候: {NDVI_PHENO_DIR}")
    print(f"  GPP物候: {GPP_PHENO_DIR}")
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
