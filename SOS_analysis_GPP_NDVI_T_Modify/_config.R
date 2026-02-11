#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ============================================================
# _config.R: 05系列R脚本统一配置
# 与Python端 _config.py 保持一致
#
# 通过环境变量控制运行模式:
#   WANG_DATA_MODE: "NDVI"(默认) / "GPP" — 植被数据类型
#   WANG_TR_MODE:   "GLEAM"(默认) / "ERA5" — 蒸腾数据源
#
# 使用方法: 在R脚本开头 source("_config.R")
# ============================================================

# ==================== 数据模式 ====================
# 从环境变量读取，默认为 "NDVI"
DATA_MODE <- Sys.getenv("WANG_DATA_MODE", unset = "NDVI")
if (!DATA_MODE %in% c("NDVI", "GPP")) {
  stop(sprintf("[_config.R] 不支持的数据模式: '%s'（仅支持 'NDVI' 或 'GPP'）", DATA_MODE))
}

# TR数据源模式（从环境变量读取，默认 "GLEAM"）
TR_MODE <- Sys.getenv("WANG_TR_MODE", unset = "GLEAM")
if (!TR_MODE %in% c("GLEAM", "ERA5")) {
  stop(sprintf("[_config.R] 不支持的TR模式: '%s'（仅支持 'GLEAM' 或 'ERA5'）", TR_MODE))
}

# ==================== 根目录（自动检测） ====================
if (.Platform$OS.type == "windows") {
  ROOT <- "I:/F/Data4"
} else {
  if (dir.exists("/mnt/i/F/Data4")) {
    ROOT <- "/mnt/i/F/Data4"
  } else {
    ROOT <- "I:/F/Data4"
  }
}

# ==================== 输出根目录 ====================
OUTPUT_ROOT <- file.path(ROOT, "Wang2025_Analysis_SOS_GPP_NDVI_GLEAM")

# 模式子目录（与Python端 apply_run_config 输出结构一致）
# 格式: OUTPUT_ROOT/<DATA_MODE>_<TR_MODE>，如 NDVI_GLEAM, GPP_ERA5
MODE_ROOT <- file.path(OUTPUT_ROOT, paste0(DATA_MODE, "_", TR_MODE))

# ==================== 核心变量名 ====================
MIDDLE_VAR_NAME <- DATA_MODE  # "NDVI" 或 "GPP"

# ==================== 模式相关路径 ====================
if (DATA_MODE == "NDVI") {
  # 物候数据 — GIMMS NDVI
  PHENO_DIR <- file.path(ROOT, "Phenology_Output_1", "NDVI_phenology")
  # 日尺度变量数据 — GIMMS NDVI
  VAR_DAILY_DIR <- file.path(ROOT, "GIMMS_NDVI", "GIMMS_NDVI_daily_interpolated")
  VAR_DAILY_PATTERN <- "NDVI_{date}.tif"
} else {
  # 物候数据 — GLASS GPP
  PHENO_DIR <- file.path(ROOT, "Phenology_Output_1", "GPP_phenology")
  # 日尺度变量数据 — GLASS GPP
  VAR_DAILY_DIR <- file.path(ROOT, "GLASS_GPP", "GLASS_GPP_daily_interpolated")
  VAR_DAILY_PATTERN <- "GPP_{date}.tif"
}

# 分解结果目录（03c模块输出，位于模式子目录下）
DECOMP_DIR <- file.path(MODE_ROOT, "Decomposition_FixedWindow")

# 文件名模式（模式相关）
FIXED_RATE_PATTERN <- sprintf("Fixed_%srate_%%d.tif", MIDDLE_VAR_NAME)
SOS_PATTERN <- sprintf("sos_%s_%%d.tif", tolower(MIDDLE_VAR_NAME))
POS_PATTERN <- sprintf("pos_doy_%s_%%d.tif", tolower(MIDDLE_VAR_NAME))

# ==================== 共享路径（不随模式变化） ====================
MASK_FILE <- file.path(OUTPUT_ROOT, "masks", "combined_mask.tif")
TEMPLATE_FILE <- file.path(OUTPUT_ROOT, "masks", "template_grid.tif")

# 气象数据路径（ERA5-Land）
PRECIP_DAILY_DIR <- file.path(ROOT, "Meteorological Data", "ERA5_Land", "Pre", "Pre_Daily", "Pre_Daily_2")
TA_DAILY_DIR <- file.path(ROOT, "Meteorological Data", "ERA5_Land", "Tem", "Tem_Daily", "Tem_Daily_2")
SW_DAILY_DIR <- file.path(ROOT, "Meteorological Data", "ERA5_Land", "DSW", "DSW_Daily", "DSW_Daily_2")

# 气象文件命名模式
PRECIP_DAILY_PATTERN <- "ERA5L_PrecipDaily_mm_{date}.tif"
TA_DAILY_PATTERN <- "ERA5L_T2mDaily_C_{date}.tif"
SW_DAILY_PATTERN <- "ERA5L_SWDaily_MJ_{date}.tif"

# ==================== 时间范围 ====================
YEAR_START <- 1982
YEAR_END <- 2018

# ==================== 全局常量 ====================
NODATA_OUT <- -9999
NODATA_ABS_MAX <- 1e20

# ==================== 05b共享数据目录 ====================
# 05b输出的Derived数据被05c/05d使用
SEM_DATA_DUAL_FIXED_DIR <- file.path(MODE_ROOT, "SEM_Data_Dual_Fixed")
DERIVED_DIR <- file.path(SEM_DATA_DUAL_FIXED_DIR, "Derived")

# 确保共享目录存在
dir.create(DERIVED_DIR, showWarnings = FALSE, recursive = TRUE)

# ==================== 配置信息输出 ====================
cat(sprintf("\n[_config.R] 数据模式: %s, TR数据源: %s\n", DATA_MODE, TR_MODE))
cat(sprintf("  ROOT:             %s\n", ROOT))
cat(sprintf("  MODE_ROOT:        %s\n", MODE_ROOT))
cat(sprintf("  TR_MODE:          %s\n", TR_MODE))
cat(sprintf("  MIDDLE_VAR_NAME:  %s\n", MIDDLE_VAR_NAME))
cat(sprintf("  PHENO_DIR:        %s\n", PHENO_DIR))
cat(sprintf("  VAR_DAILY_DIR:    %s\n", VAR_DAILY_DIR))
cat(sprintf("  DECOMP_DIR:       %s\n", DECOMP_DIR))
cat(sprintf("  DERIVED_DIR:      %s\n", DERIVED_DIR))
