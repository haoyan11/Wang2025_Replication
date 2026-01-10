#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ===============================================================
# Module 05c: 稳健SEM分析（Pooled Pixel-Year）
# ===============================================================
#
# 核心逻辑：
# 1) combined_mask 作为基础掩膜；
# 2) 全年份数据交集 + 年份比例阈值筛选有效像元；
# 3) 像元-年份 pooled 数据拟合单一 SEM；
# 4) Cluster bootstrap（按像元）给出置信区间；
# 5) 可选去趋势并输出原始/去趋势两套结果。
#
# Version: 1.2.2
# Author: Wang2025 Replication Project
# Date: 2025-01-10
#
# Changelog:
#   v1.2.2 - 修复Bug Fix 4 (2025-01-10)
#            clean_outliers函数错误过滤Fixed_Trate负值
#            影响: 从421543行恢复到~800000行（保留干旱年份数据）
#   v1.2.1 - 修复Bug Fix 3 (2025-01-05)
#            TR_fixed_window允许负值（allow_negative=TRUE）
#            影响: 从8574像元恢复到~26000像元
#   v1.2.0 - 初始版本，双时间尺度pooled SEM
# ===============================================================

suppressPackageStartupMessages({
  library(raster)
  library(lavaan)
  library(data.table)
  library(future.apply)
})

# ==================== 配置 ====================

# 根目录（自动检测）
if (.Platform$OS.type == "windows") {
  ROOT <- "I:/F/Data4"
} else {
  if (dir.exists("/mnt/i/F/Data4")) {
    ROOT <- "/mnt/i/F/Data4"
  } else {
    ROOT <- "I:/F/Data4"
  }
}

# 路径配置
OUTPUT_ROOT <- file.path(ROOT, "Wang2025_Analysis")
PHENO_DIR <- file.path(ROOT, "Phenology_Output_1", "GPP_phenology_EPSG4326")
DECOMP_DIR <- file.path(OUTPUT_ROOT, "Decomposition_FixedWindow")
MASK_FILE <- file.path(OUTPUT_ROOT, "masks", "combined_mask.tif")
TEMPLATE_FILE <- file.path(OUTPUT_ROOT, "masks", "template_grid.tif")

# 缓存数据路径（来自05b）
DERIVED_DIR <- file.path(OUTPUT_ROOT, "SEM_Data_Dual_Fixed", "Derived")

# 输出目录（与05b命名风格保持一致：SEM_Data_*/SEM_Results_*）
DATA_DIR <- file.path(OUTPUT_ROOT, "SEM_Data_Pooled_SOS")
OUTPUT_DIR <- file.path(OUTPUT_ROOT, "SEM_Results_Pooled_SOS")

# 创建输出目录
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

DATA_DIR_BASE <- DATA_DIR
OUTPUT_DIR_BASE <- OUTPUT_DIR

set_output_dirs <- function(suffix = "") {
  if (!is.null(suffix) && nzchar(suffix)) {
    DATA_DIR <<- paste0(DATA_DIR_BASE, suffix)
    OUTPUT_DIR <<- paste0(OUTPUT_DIR_BASE, suffix)
  } else {
    DATA_DIR <<- DATA_DIR_BASE
    OUTPUT_DIR <<- OUTPUT_DIR_BASE
  }
  dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
  dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
}

# 年份范围
YEAR_START <- 1982
YEAR_END <- 2018

# 常量
NODATA_OUT <- -9999
NODATA_ABS_MAX <- 1e20
FILTER_SEM_OUTLIERS <- TRUE  # 输出结果异常值过滤（类似04c）
SEM_COEF_ABS_MAX <- 5        # 系数绝对值阈值（标准化系数）
SEM_P_MIN <- 0               # p值下限
SEM_P_MAX <- 1               # p值上限
SEM_R2_MIN <- 0              # R²下限
SEM_R2_MAX <- 1              # R²上限
DETREND_ENABLE <- FALSE      # 是否启用去趋势
RUN_BOTH_DETREND <- TRUE     # 同时输出不去趋势+去趋势
DETREND_BY_PIXEL <- TRUE     # 去趋势按像元（推荐）

# 有效性阈值（与05b一致）
# 注：pooled 模式不使用窗口内日尺度有效比例过滤（仅按年份比例筛选）
MIN_VALID_YEAR_FRAC <- 0.60

# Bootstrap配置
N_BOOTSTRAP <- 800  # Bootstrap重采样次数
N_CORES <- 10  # 并行核心数（与05b一致）
MEDIATION_DENOM_EPS <- 1e-6  # 中介比例分母稳定性阈值

cat("\n======================================================================\n")
cat("稳健SEM分析 - 混合池方法（Pooled Pixel-Year Analysis）\n")
cat("======================================================================\n")
cat(sprintf("年份范围: %d-%d\n", YEAR_START, YEAR_END))
cat(sprintf("Bootstrap次数: %d\n", N_BOOTSTRAP))
cat(sprintf("并行核心数: %d\n", N_CORES))
cat("----------------------------------------------------------------------\n")

# ==================== 辅助函数 ====================

check_raster_alignment <- function(template_r, year) {
  sample_files <- list(
    list(MASK_FILE, "mask", TRUE),
    list(file.path(DECOMP_DIR, sprintf("TR_fixed_window_%d.tif", year)), "TR_fixed_window", FALSE),
    list(file.path(DECOMP_DIR, "Fixed_Window_Length.tif"), "Fixed_Window_Length", FALSE),
    list(file.path(PHENO_DIR, sprintf("sos_gpp_%d.tif", year)), "SOS", FALSE),
    list(file.path(PHENO_DIR, sprintf("pos_doy_gpp_%d.tif", year)), "POS", FALSE),
    list(file.path(DERIVED_DIR, sprintf("GPP_season_%d.tif", year)), "GPP_season", TRUE),
    list(file.path(DERIVED_DIR, sprintf("P_pre_%d.tif", year)), "P_pre", TRUE)
  )

  for (item in sample_files) {
    f <- item[[1]]
    label <- item[[2]]
    optional <- item[[3]]
    if (!file.exists(f)) {
      if (optional) {
        cat(sprintf("  ⚠ 缺少样本文件: %s (%s)\n", label, basename(f)))
        next
      }
      stop(sprintf("缺少样本文件: %s (%s)", label, basename(f)))
    }
    ok <- compareRaster(template_r, raster(f), stopiffalse = FALSE)
    if (!ok) {
      stop(sprintf("网格不一致: %s (%s)", label, basename(f)))
    } else {
      cat(sprintf("  ✓ 网格一致: %s (%s)\n", label, basename(f)))
    }
  }
}

# 检查文件是否存在
check_file <- function(file_path, description) {
  if (!file.exists(file_path)) {
    stop(sprintf("✗ 文件不存在: %s\n  路径: %s", description, file_path))
  }
  cat(sprintf("  ✓ %s\n", description))
}

# Bug Fix 3: NODATA过滤函数（类似05b的sanitize_values）
get_nodata <- function(r, fallback = NODATA_OUT) {
  nd <- NAvalue(r)
  if (is.na(nd) || !is.finite(nd)) fallback else nd
}

sanitize_values <- function(vals, na_value = -9999, allow_negative = TRUE, max_abs = NODATA_ABS_MAX) {
  # 显式NODATA值转换
  if (!is.null(na_value) && !is.na(na_value)) {
    vals[vals == na_value] <- NA
  }
  # 非有限值
  vals[!is.finite(vals)] <- NA
  # 极端异常值（绝对值过大）
  vals[abs(vals) > max_abs] <- NA
  # 可选：负值处理
  if (!allow_negative) {
    vals[vals < 0] <- NA
  }
  return(vals)
}

# 按像元去趋势（线性，逐变量）
detrend_by_pixel <- function(dt, vars, year_col = "year") {
  dt[, (vars) := {
    y <- get(year_col)
    y_mean <- mean(y)
    y_var <- sum((y - y_mean) ^ 2)
    lapply(.SD, function(x) {
      x_mean <- mean(x)
      if (!is.finite(y_var) || y_var <= 0) {
        return(x - x_mean)
      }
      slope <- sum((y - y_mean) * (x - x_mean)) / y_var
      intercept <- x_mean - slope * y_mean
      x - (intercept + slope * y)
    })
  }, by = pixel, .SDcols = vars]
  dt
}

filter_sem_param_table <- function(params) {
  info <- c(coef_extreme = 0, p_invalid = 0)
  if (!FILTER_SEM_OUTLIERS) {
    return(list(params = params, info = info))
  }

  if ("std.all" %in% names(params)) {
    bad_coef <- is.finite(params$std.all) & (abs(params$std.all) > SEM_COEF_ABS_MAX)
    info["coef_extreme"] <- sum(bad_coef)
    params$std.all[bad_coef] <- NA_real_
  }

  if ("pvalue" %in% names(params)) {
    bad_p <- is.finite(params$pvalue) & (params$pvalue < SEM_P_MIN | params$pvalue > SEM_P_MAX)
    info["p_invalid"] <- sum(bad_p)
    params$pvalue[bad_p] <- NA_real_
  }

  list(params = params, info = info)
}

filter_sem_r2 <- function(r2_vals) {
  info <- c(r2_invalid = 0)
  if (!FILTER_SEM_OUTLIERS) {
    return(list(r2 = r2_vals, info = info))
  }

  bad_r2 <- is.finite(r2_vals) & (r2_vals < SEM_R2_MIN | r2_vals > SEM_R2_MAX)
  info["r2_invalid"] <- sum(bad_r2)
  r2_vals[bad_r2] <- NA_real_
  list(r2 = r2_vals, info = info)
}

filter_boot_matrix <- function(boot_mat) {
  info <- c(coef_extreme = 0)
  if (!FILTER_SEM_OUTLIERS) {
    return(list(mat = boot_mat, info = info))
  }

  bad_coef <- is.finite(boot_mat) & (abs(boot_mat) > SEM_COEF_ABS_MAX)
  info["coef_extreme"] <- sum(bad_coef)
  boot_mat[bad_coef] <- NA_real_
  list(mat = boot_mat, info = info)
}

# 中介比例安全计算（防止分母接近0导致爆炸）
safe_mediation_ratio <- function(c_val, d_val, e_val, eps = MEDIATION_DENOM_EPS) {
  denom <- e_val + c_val * d_val
  if (!is.finite(denom) || abs(denom) < eps) {
    return(NA_real_)
  }
  (c_val * d_val) / denom
}

# ==================== 步骤1：验证缓存文件 ====================

cat("\n=== 步骤1：验证缓存文件 ===\n")
cat("检查05b生成的缓存文件...\n")

years <- YEAR_START:YEAR_END

# 检查必需的文件（注意：mask文件作为基础掩膜，最终仍需数据交集筛选）
cat("  注意：使用combined_mask作为基础掩膜，有效像元由数据交集+年份比例阈值确定\n")
check_file(file.path(DECOMP_DIR, "Fixed_Window_Length.tif"), "固定窗口长度")

# 检查所有年份的缓存文件
missing_files <- 0
for (year in years) {
  files_to_check <- c(
    sprintf("P_pre_%d.tif", year),
    sprintf("T_pre_%d.tif", year),
    sprintf("SW_pre_%d.tif", year),
    sprintf("P_season_%d.tif", year),
    sprintf("T_season_%d.tif", year),
    sprintf("SW_season_%d.tif", year),
    sprintf("GPP_season_%d.tif", year)
  )

  for (f in files_to_check) {
    file_path <- file.path(DERIVED_DIR, f)
    if (!file.exists(file_path)) {
      cat(sprintf("  ✗ 缺失: %s\n", f))
      missing_files <- missing_files + 1
    }
  }

  # 检查TR_fixed_window
  tr_file <- file.path(DECOMP_DIR, sprintf("TR_fixed_window_%d.tif", year))
  if (!file.exists(tr_file)) {
    cat(sprintf("  ✗ 缺失: TR_fixed_window_%d.tif\n", year))
    missing_files <- missing_files + 1
  }

  # 检查SOS
  sos_file <- file.path(PHENO_DIR, sprintf("sos_gpp_%d.tif", year))
  if (!file.exists(sos_file)) {
    cat(sprintf("  ✗ 缺失: sos_gpp_%d.tif\n", year))
    missing_files <- missing_files + 1
  }
}

if (missing_files > 0) {
  stop(sprintf("\n✗ 错误：缺失 %d 个文件！\n请先运行05b代码生成所有缓存文件。", missing_files))
}

cat(sprintf("✓ 所有 %d 年的缓存文件验证通过\n", length(years)))

# 网格一致性硬检查（模板/样本栅格）
cat("\n[网格一致性检查]\n")
template_path <- if (file.exists(TEMPLATE_FILE)) {
  TEMPLATE_FILE
} else if (file.exists(file.path(DECOMP_DIR, sprintf("TR_fixed_window_%d.tif", YEAR_START)))) {
  file.path(DECOMP_DIR, sprintf("TR_fixed_window_%d.tif", YEAR_START))
} else {
  MASK_FILE
}
template_r <- raster(template_path)
cat(sprintf("  模板文件: %s\n", template_path))
check_raster_alignment(template_r, YEAR_START)

# ==================== 步骤2：提取像元数据（Bug Fix 1: 使用全部37年数据）====================

cat("\n=== 步骤2：动态确定有效像元（掩膜 + 年份比例阈值）===\n")
cat("说明：使用combined_mask作为基础掩膜，并在全时段交集中筛选有效像元\n")
cat("      这确保像元在所有变量和大部分年份都有有效值\n")
cat("      （与05b阈值逻辑保持一致）\n\n")

# Bug Fix 1优化: 两阶段筛选，避免读取所有像元×所有年份
# 阶段1: 用1982年快速筛选候选像元（排除99%无效像元）
# 阶段2: 只对候选像元检查全部37年是否完整

cat("  【阶段1】用1982年快速筛选候选像元...\n")

# 读取fixed_len（所有年份共用）
fixed_len_r <- raster(file.path(DECOMP_DIR, "Fixed_Window_Length.tif"))
fixed_len_vals_all <- getValues(fixed_len_r)
fixed_len_nodata <- get_nodata(fixed_len_r)
n_cells_total <- length(fixed_len_vals_all)

cat(sprintf("    栅格总像元数: %d\n", n_cells_total))

# 读取掩膜并限制有效像元
mask_r <- raster(MASK_FILE)
mask_vals_all <- getValues(mask_r)
if (length(mask_vals_all) != n_cells_total) {
  stop("combined_mask.tif 的像元数与模板不一致，无法继续。")
}
mask_valid <- is.finite(mask_vals_all) & (mask_vals_all > 0)
cat(sprintf("    掩膜有效像元数: %d (%.2f%%)\n",
            sum(mask_valid), 100 * sum(mask_valid) / n_cells_total))

# 读取1982年数据快速筛选
year_test <- years[1]  # 1982
tr_test_r <- raster(file.path(DECOMP_DIR, sprintf("TR_fixed_window_%d.tif", year_test)))
sos_test_r <- raster(file.path(PHENO_DIR, sprintf("sos_gpp_%d.tif", year_test)))
gpp_test_r <- raster(file.path(DERIVED_DIR, sprintf("GPP_season_%d.tif", year_test)))
p_pre_test_r <- raster(file.path(DERIVED_DIR, sprintf("P_pre_%d.tif", year_test)))
t_pre_test_r <- raster(file.path(DERIVED_DIR, sprintf("T_pre_%d.tif", year_test)))
sw_pre_test_r <- raster(file.path(DERIVED_DIR, sprintf("SW_pre_%d.tif", year_test)))
p_season_test_r <- raster(file.path(DERIVED_DIR, sprintf("P_season_%d.tif", year_test)))
t_season_test_r <- raster(file.path(DERIVED_DIR, sprintf("T_season_%d.tif", year_test)))
sw_season_test_r <- raster(file.path(DERIVED_DIR, sprintf("SW_season_%d.tif", year_test)))

# 栅格对齐一致性检查（若不一致则直接停止）
cat("    对齐检查: 固定窗口/物候/缓存栅格...\n")
compareRaster(
  fixed_len_r, tr_test_r, sos_test_r, gpp_test_r,
  p_pre_test_r, t_pre_test_r, sw_pre_test_r,
  p_season_test_r, t_season_test_r, sw_season_test_r,
  stopiffalse = TRUE
)
cat("    ✓ 栅格对齐检查通过\n")

# Bug Fix 3: 应用NODATA过滤
tr_test <- sanitize_values(getValues(tr_test_r), get_nodata(tr_test_r), allow_negative = TRUE)  # TR_fixed_window 可以是负值！
sos_test <- sanitize_values(getValues(sos_test_r), get_nodata(sos_test_r), allow_negative = FALSE)
gpp_test <- sanitize_values(getValues(gpp_test_r), get_nodata(gpp_test_r), allow_negative = FALSE)
p_pre_test <- sanitize_values(getValues(p_pre_test_r), get_nodata(p_pre_test_r), allow_negative = FALSE)
t_pre_test <- sanitize_values(getValues(t_pre_test_r), get_nodata(t_pre_test_r))
sw_pre_test <- sanitize_values(getValues(sw_pre_test_r), get_nodata(sw_pre_test_r), allow_negative = FALSE)
p_season_test <- sanitize_values(getValues(p_season_test_r), get_nodata(p_season_test_r), allow_negative = FALSE)
t_season_test <- sanitize_values(getValues(t_season_test_r), get_nodata(t_season_test_r))
sw_season_test <- sanitize_values(getValues(sw_season_test_r), get_nodata(sw_season_test_r), allow_negative = FALSE)
fixed_len_test <- sanitize_values(fixed_len_vals_all, fixed_len_nodata, allow_negative = FALSE)

# 变量特定约束
sos_test[sos_test < 1 | sos_test > 365] <- NA
# tr_test 可以是负值，不过滤！
gpp_test[gpp_test <= 0] <- NA
fixed_len_test[fixed_len_test <= 0] <- NA

# 【诊断】逐个变量统计有效像元数
cat("\n    【诊断】1982年各变量有效像元统计：\n")
cat(sprintf("      TR:         %6d (%.2f%%)\n", sum(!is.na(tr_test)), 100*sum(!is.na(tr_test))/n_cells_total))
cat(sprintf("      SOS:        %6d (%.2f%%)\n", sum(!is.na(sos_test)), 100*sum(!is.na(sos_test))/n_cells_total))
cat(sprintf("      GPP:        %6d (%.2f%%)\n", sum(!is.na(gpp_test)), 100*sum(!is.na(gpp_test))/n_cells_total))
cat(sprintf("      P_pre:      %6d (%.2f%%)\n", sum(!is.na(p_pre_test)), 100*sum(!is.na(p_pre_test))/n_cells_total))
cat(sprintf("      T_pre:      %6d (%.2f%%)\n", sum(!is.na(t_pre_test)), 100*sum(!is.na(t_pre_test))/n_cells_total))
cat(sprintf("      SW_pre:     %6d (%.2f%%)\n", sum(!is.na(sw_pre_test)), 100*sum(!is.na(sw_pre_test))/n_cells_total))
cat(sprintf("      P_season:   %6d (%.2f%%)\n", sum(!is.na(p_season_test)), 100*sum(!is.na(p_season_test))/n_cells_total))
cat(sprintf("      T_season:   %6d (%.2f%%)\n", sum(!is.na(t_season_test)), 100*sum(!is.na(t_season_test))/n_cells_total))
cat(sprintf("      SW_season:  %6d (%.2f%%)\n", sum(!is.na(sw_season_test)), 100*sum(!is.na(sw_season_test))/n_cells_total))
cat(sprintf("      Fixed_len:  %6d (%.2f%%)\n", sum(!is.na(fixed_len_test)), 100*sum(!is.na(fixed_len_test))/n_cells_total))

# 核心3变量交集（类似05b）
valid_core <- !is.na(tr_test) & !is.na(sos_test) & !is.na(gpp_test) &
              !is.na(fixed_len_test) & is.finite(tr_test / fixed_len_test) &
              mask_valid
cat(sprintf("\n      核心3变量交集: %6d (%.2f%%) ← 这应该接近05b的28000\n",
            sum(valid_core), 100*sum(valid_core)/n_cells_total))

# 1982年的候选像元（在该年所有变量都有效）
valid_1982 <- !is.na(tr_test) & !is.na(sos_test) & !is.na(gpp_test) &
              !is.na(p_pre_test) & !is.na(t_pre_test) & !is.na(sw_pre_test) &
              !is.na(p_season_test) & !is.na(t_season_test) & !is.na(sw_season_test) &
              !is.na(fixed_len_test) & is.finite(tr_test / fixed_len_test) &
              mask_valid

candidate_cells <- which(valid_1982)
n_candidates <- length(candidate_cells)

cat(sprintf("    1982年有效候选像元: %d (占比 %.2f%%)\n",
            n_candidates, 100 * n_candidates / n_cells_total))

# 清理测试数据
rm(tr_test_r, sos_test_r, gpp_test_r, p_pre_test_r, t_pre_test_r,
   sw_pre_test_r, p_season_test_r, t_season_test_r, sw_season_test_r,
   tr_test, sos_test, gpp_test, p_pre_test, t_pre_test, sw_pre_test,
   p_season_test, t_season_test, sw_season_test, fixed_len_test, valid_1982)
gc()

# 【阶段2】只对候选像元检查全部37年
cat(sprintf("\n  【阶段2】对候选像元检查全部37年数据完整性...\n"))
cat(sprintf("    需要检查: %d 像元 × %d 年 = %d 数据点\n",
            n_candidates, length(years), n_candidates * length(years)))
cat(sprintf("    （相比全量检查减少了 %.1f%% 的计算量）\n\n",
            100 * (1 - n_candidates / n_cells_total)))

# 只读取候选像元的数据
year_list <- vector("list", length(years))

for (i in seq_along(years)) {
  year <- years[i]

  if (i %% 5 == 0 || i == 1 || i == length(years)) {
    cat(sprintf("    进度: %d/%d 年份 (%.1f%%) - %d\n",
                i, length(years), 100 * i / length(years), year))
  }

  # 读取该年的所有变量
  tr_fixed_r <- raster(file.path(DECOMP_DIR, sprintf("TR_fixed_window_%d.tif", year)))
  sos_r <- raster(file.path(PHENO_DIR, sprintf("sos_gpp_%d.tif", year)))
  gpp_season_r <- raster(file.path(DERIVED_DIR, sprintf("GPP_season_%d.tif", year)))
  p_pre_r <- raster(file.path(DERIVED_DIR, sprintf("P_pre_%d.tif", year)))
  t_pre_r <- raster(file.path(DERIVED_DIR, sprintf("T_pre_%d.tif", year)))
  sw_pre_r <- raster(file.path(DERIVED_DIR, sprintf("SW_pre_%d.tif", year)))
  p_season_r <- raster(file.path(DERIVED_DIR, sprintf("P_season_%d.tif", year)))
  t_season_r <- raster(file.path(DERIVED_DIR, sprintf("T_season_%d.tif", year)))
  sw_season_r <- raster(file.path(DERIVED_DIR, sprintf("SW_season_%d.tif", year)))

  # Bug Fix 3: 只提取候选像元的值并过滤NODATA
  tr_vals <- sanitize_values(getValues(tr_fixed_r)[candidate_cells], get_nodata(tr_fixed_r), allow_negative = TRUE)  # TR_fixed_window 可以是负值！
  sos_vals <- sanitize_values(getValues(sos_r)[candidate_cells], get_nodata(sos_r), allow_negative = FALSE)
  gpp_vals <- sanitize_values(getValues(gpp_season_r)[candidate_cells], get_nodata(gpp_season_r), allow_negative = FALSE)
  p_pre_vals <- sanitize_values(getValues(p_pre_r)[candidate_cells], get_nodata(p_pre_r), allow_negative = FALSE)
  t_pre_vals <- sanitize_values(getValues(t_pre_r)[candidate_cells], get_nodata(t_pre_r))
  sw_pre_vals <- sanitize_values(getValues(sw_pre_r)[candidate_cells], get_nodata(sw_pre_r), allow_negative = FALSE)
  p_season_vals <- sanitize_values(getValues(p_season_r)[candidate_cells], get_nodata(p_season_r), allow_negative = FALSE)
  t_season_vals <- sanitize_values(getValues(t_season_r)[candidate_cells], get_nodata(t_season_r))
  sw_season_vals <- sanitize_values(getValues(sw_season_r)[candidate_cells], get_nodata(sw_season_r), allow_negative = FALSE)

  # 额外的变量特定约束
  sos_vals[sos_vals < 1 | sos_vals > 365] <- NA
  # tr_vals 可以是负值，不过滤！
  gpp_vals[gpp_vals <= 0] <- NA

  # 计算Fixed_Trate（只对候选像元）
  fixed_len_vals <- sanitize_values(fixed_len_vals_all[candidate_cells], fixed_len_nodata, allow_negative = FALSE)
  fixed_len_vals[fixed_len_vals <= 0] <- NA
  fixed_trate_vals <- tr_vals / fixed_len_vals

  # 构建该年的数据（只包含候选像元）
  year_data <- data.table(
    pixel = candidate_cells,
    year = year,
    Fixed_Trate = fixed_trate_vals,
    SOS = sos_vals,
    GPP_season = gpp_vals,
    P_pre = p_pre_vals,
    T_pre = t_pre_vals,
    SW_pre = sw_pre_vals,
    P_season = p_season_vals,
    T_season = t_season_vals,
    SW_season = sw_season_vals
  )

  # 缓存到列表，最后一次性合并
  year_list[[i]] <- year_data
}

long_df_temp <- rbindlist(year_list, use.names = TRUE)
rm(year_list)

cat("\n  ✓ 候选像元数据读取完成\n")
cat(sprintf("    总行数: %d\n", nrow(long_df_temp)))

# 使用像元分组计数确定在所有37年都有效的像元（避免dcast超宽表）
cat("\n  正在筛选满足有效年份比例的像元...\n")
required_years <- ceiling(MIN_VALID_YEAR_FRAC * length(years))
cat(sprintf("  要求有效年份数: %d (%.0f%%)\n",
            required_years, 100 * MIN_VALID_YEAR_FRAC))

row_complete <- complete.cases(long_df_temp[, .(
  Fixed_Trate, SOS, GPP_season,
  P_pre, T_pre, SW_pre,
  P_season, T_season, SW_season
)])
long_df_temp[, row_complete := row_complete]

pixel_counts <- long_df_temp[, .(n_valid = sum(row_complete)), by = pixel]
complete_pixels <- pixel_counts[n_valid >= required_years, pixel]
n_pixels <- length(complete_pixels)

cat(sprintf("  ✓ 有效像元数: %d\n", n_pixels))
cat(sprintf("    占候选像元比例: %.2f%% (%d/%d)\n",
            100 * n_pixels / n_candidates, n_pixels, n_candidates))
cat(sprintf("    占总像元比例: %.2f%%\n", 100 * n_pixels / n_cells_total))
cat(sprintf("    年份数: %d\n", length(years)))
cat(sprintf("    最终数据行数(含NA行前): %d\n", n_pixels * length(years)))

# 提取完整像元的数据
setkey(long_df_temp, pixel)
complete_pixels_dt <- data.table(pixel = complete_pixels)
long_df <- long_df_temp[complete_pixels_dt, on = "pixel", nomatch = 0L]
long_df <- long_df[row_complete == TRUE]
long_df[, row_complete := NULL]

# 清理临时数据
rm(long_df_temp, pixel_counts, complete_pixels_dt)
gc()

cat(sprintf("\n  ✓ 数据提取完成\n"))
cat(sprintf("    总行数: %d\n", nrow(long_df)))
cat(sprintf("    总列数: %d\n", ncol(long_df)))

# ==================== 步骤3：数据清洗 ====================

cat("\n=== 步骤3：数据清洗 ===\n")

# 统计原始数据
n_before <- nrow(long_df)
cat(sprintf("  原始数据行数: %d\n", n_before))

# 移除包含NA或无效值的行
# 清洗NA值
long_df_clean <- long_df[complete.cases(long_df)]
n_after_na <- nrow(long_df_clean)
cat(sprintf("  移除NA后: %d 行 (移除 %d 行, %.2f%%)\n",
            n_after_na, n_before - n_after_na,
            100 * (n_before - n_after_na) / n_before))

# 清洗极端值（可选，根据实际情况调整）
# 这里采用保守策略，只移除明显异常的值
clean_outliers <- function(dt) {
  # Fixed_Trate可以是负值（干旱年份），只过滤Inf/NaN
  dt <- dt[is.finite(Fixed_Trate)]  # Bug Fix 4: 允许负值！
  # SOS应该在1-365之间
  dt <- dt[SOS >= 1 & SOS <= 365]
  # GPP应该为正
  dt <- dt[GPP_season > 0 & is.finite(GPP_season)]
  # 降水应该非负
  dt <- dt[P_pre >= 0 & P_season >= 0]

  return(dt)
}

long_df_clean <- clean_outliers(long_df_clean)
n_final <- nrow(long_df_clean)

cat(sprintf("  清洗极端值后: %d 行 (移除 %d 行, %.2f%%)\n",
            n_final, n_after_na - n_final,
            100 * (n_after_na - n_final) / n_after_na))

cat(sprintf("\n  ✓ 最终有效数据: %d 行\n", n_final))
cat(sprintf("    数据完整率: %.2f%%\n", 100 * n_final / n_before))

# 统计每个变量的基本信息
cat("\n  变量基本统计:\n")
vars_to_summarize <- c("Fixed_Trate", "SOS", "GPP_season",
                       "P_pre", "T_pre", "SW_pre",
                       "P_season", "T_season", "SW_season")

for (var in vars_to_summarize) {
  vals <- long_df_clean[[var]]
  cat(sprintf("    %s: 均值=%.3f, 标准差=%.3f, 范围=[%.3f, %.3f]\n",
              var, mean(vals, na.rm = TRUE), sd(vals, na.rm = TRUE),
              min(vals, na.rm = TRUE), max(vals, na.rm = TRUE)))
}

run_sem_pipeline <- function(run_label, detrend_enable) {
  cat(sprintf("\n======================================================================\n"))
  cat(sprintf("运行: %s\n", run_label))
  cat("======================================================================\n")

  long_df_use <- copy(long_df_clean)
  if (detrend_enable) {
    cat("去趋势: 开启\n")
    if (DETREND_BY_PIXEL) {
      long_df_use <- detrend_by_pixel(long_df_use, vars_to_summarize)
    } else {
      y <- long_df_use$year
      y_mean <- mean(y)
      y_var <- sum((y - y_mean) ^ 2)
      for (var in vars_to_summarize) {
        x <- long_df_use[[var]]
        x_mean <- mean(x)
        if (!is.finite(y_var) || y_var <= 0) {
          long_df_use[[var]] <- x - x_mean
          next
        }
        slope <- sum((y - y_mean) * (x - x_mean)) / y_var
        intercept <- x_mean - slope * y_mean
        long_df_use[[var]] <- x - (intercept + slope * y)
      }
    }
  } else {
    cat("去趋势: 关闭\n")
  }

  # 保存清洗后的数据
  write.csv(long_df_use, file.path(DATA_DIR, "sem_pooled_raw.csv"), row.names = FALSE)
  cat(sprintf("\n  ✓ 原始数据已保存: %s\n", file.path(DATA_DIR, "sem_pooled_raw.csv")))

  # ==================== 步骤4：标准化 ====================
  cat("\n=== 步骤4：Z-score标准化 ===\n")
  cat("对所有变量进行Z-score标准化（均值=0，标准差=1）...\n")

  long_df_std <- copy(long_df_use)

  for (var in vars_to_summarize) {
    vals <- long_df_std[[var]]
    vals_std <- (vals - mean(vals, na.rm = TRUE)) / sd(vals, na.rm = TRUE)
    long_df_std[[var]] <- vals_std
  }

  cat("  ✓ 标准化完成\n")

  # 验证标准化
  cat("\n  标准化后验证（应接近0和1）:\n")
  for (var in vars_to_summarize) {
    vals <- long_df_std[[var]]
    cat(sprintf("    %s: 均值=%.6f, 标准差=%.6f\n",
                var, mean(vals, na.rm = TRUE), sd(vals, na.rm = TRUE)))
  }

  # 保存标准化数据
  write.csv(long_df_std, file.path(DATA_DIR, "sem_pooled_standardized.csv"), row.names = FALSE)
  cat(sprintf("\n  ✓ 标准化数据已保存: %s\n", file.path(DATA_DIR, "sem_pooled_standardized.csv")))

  # ==================== 步骤5：构建SEM模型 ====================
  cat("\n=== 步骤5：构建SEM模型 ===\n")
  cat("双时间尺度完整路径SEM模型...\n")

  sem_model <- '
    # ===== 第一层：季前气候 → SOS（物候响应）=====
    SOS ~ a1*P_pre + a2*T_pre + a3*SW_pre

    # ===== 第二层：SOS + 季前气候 + 生长季气候 → GPP（碳固定）=====
    GPP_season ~ b*SOS + f1*P_pre + f2*T_pre + f3*SW_pre +
                 c1*P_season + c2*T_season + c3*SW_season

    # ===== 第三层：SOS + GPP + 季前气候 + 生长季气候 → Fixed_Trate（蒸腾）=====
    Fixed_Trate ~ g*SOS + d*GPP_season +
                  h1*P_pre + h2*T_pre + h3*SW_pre +
                  e1*P_season + e2*T_season + e3*SW_season

    # ===== 间接效应分解 =====

    # 季前气候的间接效应路径
    P_pre_via_SOS_GPP  := a1 * b * d
    T_pre_via_SOS_GPP  := a2 * b * d
    SW_pre_via_SOS_GPP := a3 * b * d

    P_pre_via_SOS  := a1 * g
    T_pre_via_SOS  := a2 * g
    SW_pre_via_SOS := a3 * g

    P_pre_via_GPP  := f1 * d
    T_pre_via_GPP  := f2 * d
    SW_pre_via_GPP := f3 * d

    # 季前气候的总间接效应
    P_pre_indirect  := a1*b*d + a1*g + f1*d
    T_pre_indirect  := a2*b*d + a2*g + f2*d
    SW_pre_indirect := a3*b*d + a3*g + f3*d

    # 生长季气候通过GPP的间接效应
    P_season_via_GPP  := c1 * d
    T_season_via_GPP  := c2 * d
    SW_season_via_GPP := c3 * d

    # GPP的中介比例
    P_GPP_mediation := (c1*d) / (e1 + c1*d)
    T_GPP_mediation := (c2*d) / (e2 + c2*d)
    SW_GPP_mediation := (c3*d) / (e3 + c3*d)
  '

  cat("\n模型结构:\n")
  cat("  第一层: 季前气候 → SOS\n")
  cat("  第二层: SOS + 季前气候 + 生长季气候 → GPP\n")
  cat("  第三层: SOS + GPP + 季前气候 + 生长季气候 → Fixed_Trate\n")
  cat("\n  包含完整的直接效应和间接效应路径\n")

  # ==================== 步骤6：拟合SEM模型 ====================
  cat("\n=== 步骤6：拟合SEM模型 ===\n")
  cat(sprintf("使用 %d 个观测点拟合SEM...\n", nrow(long_df_std)))

  fit <- sem(sem_model, data = as.data.frame(long_df_std), estimator = "MLR")

  cat("\n✓ SEM模型拟合完成\n")

  # 提取结果
  cat("\n模型拟合摘要:\n")
  summary(fit, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

  # 保存拟合结果
  fit_summary <- summary(fit, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
  sink(file.path(OUTPUT_DIR, "SEM_pooled_summary.txt"))
  print(fit_summary)
  sink()

  # 提取参数估计
  params <- parameterEstimates(fit, standardized = TRUE)
  param_filter <- filter_sem_param_table(params)
  params <- param_filter$params
  if (FILTER_SEM_OUTLIERS) {
    cat(sprintf("  输出异常值过滤（全域SEM）：系数极端=%d, p值异常=%d\n",
                param_filter$info["coef_extreme"], param_filter$info["p_invalid"]))
  }

  # 修正中介比例（分母过小置为NA）
  get_label_est <- function(params_df, label) {
    v <- params_df$est[params_df$label == label]
    if (length(v) == 0) return(NA_real_)
    v[1]
  }

  d_val <- get_label_est(params, "d")
  p_ratio <- safe_mediation_ratio(get_label_est(params, "c1"), d_val, get_label_est(params, "e1"))
  t_ratio <- safe_mediation_ratio(get_label_est(params, "c2"), d_val, get_label_est(params, "e2"))
  sw_ratio <- safe_mediation_ratio(get_label_est(params, "c3"), d_val, get_label_est(params, "e3"))

  update_ratio <- function(params_df, label, value) {
    idx <- which(params_df$label == label)
    if (length(idx) > 0) {
      params_df$est[idx] <- value
      params_df$std.all[idx] <- NA_real_
    }
    params_df
  }

  params <- update_ratio(params, "P_GPP_mediation", p_ratio)
  params <- update_ratio(params, "T_GPP_mediation", t_ratio)
  params <- update_ratio(params, "SW_GPP_mediation", sw_ratio)

  if (any(is.na(c(p_ratio, t_ratio, sw_ratio)))) {
    cat(sprintf("  ⚠️ 中介比例分母过小，已置NA: P=%s, T=%s, SW=%s\n",
                ifelse(is.na(p_ratio), "NA", "OK"),
                ifelse(is.na(t_ratio), "NA", "OK"),
                ifelse(is.na(sw_ratio), "NA", "OK")))
  }
  write.csv(params, file.path(OUTPUT_DIR, "SEM_pooled_parameters.csv"), row.names = FALSE)

  # 提取拟合指标
  fit_measures <- fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "tli",
                                     "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr"))
  write.csv(t(as.data.frame(fit_measures)), file.path(OUTPUT_DIR, "SEM_pooled_fitindices.csv"))

  # 提取R²
  r2_vals <- inspect(fit, "r2")
  r2_filter <- filter_sem_r2(r2_vals)
  r2_vals <- r2_filter$r2
  if (FILTER_SEM_OUTLIERS) {
    cat(sprintf("  输出异常值过滤（全域SEM）：R²异常=%d\n",
                r2_filter$info["r2_invalid"]))
  }
  write.csv(as.data.frame(r2_vals), file.path(OUTPUT_DIR, "SEM_pooled_R2.csv"))

  cat(sprintf("\n  ✓ 结果已保存到: %s\n", OUTPUT_DIR))

  # ==================== 步骤7：Cluster Bootstrap ====================
  cat("\n=== 步骤7：Cluster Bootstrap 置信区间 ===\n")
  cat(sprintf("Bootstrap重采样: %d 次\n", N_BOOTSTRAP))
  cat(sprintf("并行核心数: %d\n", N_CORES))
  cat("重采样单元: 像元（保留时间结构）\n\n")

  # 进度条支持（可选）
  use_progressr <- requireNamespace("progressr", quietly = TRUE)
  if (use_progressr) {
    progressr::handlers("txtprogressbar")
    progressr::handlers(global = TRUE)
    options(progressr.enable = TRUE)
    options(progressr.intrusive = TRUE)
    Sys.setenv(R_PROGRESSR_ENABLE = "TRUE")
    cat("  进度显示: progressr (txtprogressbar)\n")
  } else {
    cat("  进度显示: 每50次打印一次\n")
  }

  # 设置并行
  plan(multisession, workers = N_CORES)

  unique_pixels <- unique(long_df_std$pixel)
  n_pix_total <- length(unique_pixels)

  cat(sprintf("  总像元数: %d\n", n_pix_total))
  cat(sprintf("  每个像元年份数: %d\n", length(unique(long_df_std$year))))
  cat("\n  开始Bootstrap重采样...\n")

  # 为join加速设置键
  setkey(long_df_std, pixel)

  if (use_progressr) {
    boot_results <- progressr::with_progress({
      p <- progressr::progressor(along = seq_len(N_BOOTSTRAP))
      future_lapply(seq_len(N_BOOTSTRAP), function(i) {
        p()

        # 有放回抽样像元（保留每个像元的时间结构）
        boot_pixels <- sample(unique_pixels, n_pix_total, replace = TRUE)

        # Bug Fix 2: 提取这些像元的所有年份数据（保留重复采样）
        # 错误方法: boot_df <- long_df_std[pixel %in% boot_pixels]  # %in% 会去重！
        # 正确方法: 使用lapply逐个提取，保留重复
        boot_i <- data.table(pixel = boot_pixels)
        boot_df <- long_df_std[boot_i, on = "pixel", allow.cartesian = TRUE, nomatch = 0L]
        boot_df <- as.data.frame(boot_df[, .(
          Fixed_Trate, SOS, GPP_season,
          P_pre, T_pre, SW_pre,
          P_season, T_season, SW_season
        )])

        # 拟合SEM
        fit_boot <- try(sem(sem_model, data = boot_df, estimator = "ML", se = "none"), silent = TRUE)

        if (inherits(fit_boot, "try-error")) {
          return(list(ok = FALSE, reason = "try-error", msg = as.character(fit_boot)))
        }
        if (!lavInspect(fit_boot, "converged")) {
          return(list(ok = FALSE, reason = "not-converged", msg = "not converged"))
        }

        # 提取标准化系数
        pe <- try(parameterEstimates(fit_boot, standardized = TRUE), silent = TRUE)
        if (inherits(pe, "try-error")) {
          return(list(ok = FALSE, reason = "param-error", msg = as.character(pe)))
        }

        # 提取有标签的参数（路径系数和间接效应）
        labeled_params <- pe[pe$label != "" | pe$op == ":=", ]
        coefs <- ifelse(labeled_params$op == ":=", labeled_params$est, labeled_params$std.all)
        names(coefs) <- ifelse(labeled_params$label != "",
                               labeled_params$label,
                               labeled_params$lhs)

        # 修正中介比例（分母过小置为NA）
        boot_get_label_est <- function(pe_df, label) {
          v <- pe_df$est[pe_df$label == label]
          if (length(v) == 0) return(NA_real_)
          v[1]
        }
        d_boot <- boot_get_label_est(pe, "d")
        coefs["P_GPP_mediation"] <- safe_mediation_ratio(boot_get_label_est(pe, "c1"), d_boot, boot_get_label_est(pe, "e1"))
        coefs["T_GPP_mediation"] <- safe_mediation_ratio(boot_get_label_est(pe, "c2"), d_boot, boot_get_label_est(pe, "e2"))
        coefs["SW_GPP_mediation"] <- safe_mediation_ratio(boot_get_label_est(pe, "c3"), d_boot, boot_get_label_est(pe, "e3"))

        # 提取R²
        r2_vals <- lavInspect(fit_boot, "r2")

        list(ok = TRUE, coefs = coefs, r2 = r2_vals)
      }, future.seed = TRUE)
    })
  } else {
    boot_results <- future_lapply(seq_len(N_BOOTSTRAP), function(i) {
      if (i %% 50 == 0 || i == 1) {
        cat(sprintf("    Bootstrap迭代: %d/%d (%.1f%%)\n", i, N_BOOTSTRAP, 100 * i / N_BOOTSTRAP))
      }

      # 有放回抽样像元（保留每个像元的时间结构）
      boot_pixels <- sample(unique_pixels, n_pix_total, replace = TRUE)

      # Bug Fix 2: 提取这些像元的所有年份数据（保留重复采样）
      # 错误方法: boot_df <- long_df_std[pixel %in% boot_pixels]  # %in% 会去重！
      # 正确方法: 使用lapply逐个提取，保留重复
      boot_i <- data.table(pixel = boot_pixels)
      boot_df <- long_df_std[boot_i, on = "pixel", allow.cartesian = TRUE, nomatch = 0L]
      boot_df <- as.data.frame(boot_df[, .(
        Fixed_Trate, SOS, GPP_season,
        P_pre, T_pre, SW_pre,
        P_season, T_season, SW_season
      )])

      # 拟合SEM
      fit_boot <- try(sem(sem_model, data = boot_df, estimator = "ML", se = "none"), silent = TRUE)

      if (inherits(fit_boot, "try-error")) {
        return(list(ok = FALSE, reason = "try-error", msg = as.character(fit_boot)))
      }
      if (!lavInspect(fit_boot, "converged")) {
        return(list(ok = FALSE, reason = "not-converged", msg = "not converged"))
      }

      # 提取标准化系数
      pe <- try(parameterEstimates(fit_boot, standardized = TRUE), silent = TRUE)
      if (inherits(pe, "try-error")) {
        return(list(ok = FALSE, reason = "param-error", msg = as.character(pe)))
      }

      # 提取有标签的参数（路径系数和间接效应）
      labeled_params <- pe[pe$label != "" | pe$op == ":=", ]
      coefs <- ifelse(labeled_params$op == ":=", labeled_params$est, labeled_params$std.all)
      names(coefs) <- ifelse(labeled_params$label != "",
                             labeled_params$label,
                             labeled_params$lhs)

      # 修正中介比例（分母过小置为NA）
      boot_get_label_est <- function(pe_df, label) {
        v <- pe_df$est[pe_df$label == label]
        if (length(v) == 0) return(NA_real_)
        v[1]
      }
      d_boot <- boot_get_label_est(pe, "d")
      coefs["P_GPP_mediation"] <- safe_mediation_ratio(boot_get_label_est(pe, "c1"), d_boot, boot_get_label_est(pe, "e1"))
      coefs["T_GPP_mediation"] <- safe_mediation_ratio(boot_get_label_est(pe, "c2"), d_boot, boot_get_label_est(pe, "e2"))
      coefs["SW_GPP_mediation"] <- safe_mediation_ratio(boot_get_label_est(pe, "c3"), d_boot, boot_get_label_est(pe, "e3"))

      # 提取R²
      r2_vals <- lavInspect(fit_boot, "r2")

      list(ok = TRUE, coefs = coefs, r2 = r2_vals)
    }, future.seed = TRUE)
  }

  # 关闭并行
  plan(sequential)

  cat("\n  ✓ Bootstrap重采样完成\n")

  # 移除失败的迭代并统计原因
  valid_boots <- Filter(function(x) isTRUE(x$ok), boot_results)
  failed_boots <- Filter(function(x) !isTRUE(x$ok), boot_results)
  n_valid <- length(valid_boots)
  n_failed <- length(failed_boots)

  cat(sprintf("\n  有效Bootstrap样本: %d/%d (%.1f%%)\n",
              n_valid, N_BOOTSTRAP, 100 * n_valid / N_BOOTSTRAP))

  if (n_failed > 0) {
    reason_counts <- sort(table(vapply(failed_boots, function(x) x$reason, character(1))), decreasing = TRUE)
    cat("  失败原因统计:\n")
    for (reason in names(reason_counts)) {
      cat(sprintf("    - %s: %d\n", reason, reason_counts[[reason]]))
    }

    boot_log_path <- file.path(OUTPUT_DIR, "SEM_pooled_bootstrap_failures.log")
    fail_msgs <- vapply(failed_boots, function(x) x$msg, character(1))
    fail_msgs <- fail_msgs[!is.na(fail_msgs) & nzchar(fail_msgs)]

    log_lines <- c(
      sprintf("Bootstrap失败样本: %d/%d (%.1f%%)", n_failed, N_BOOTSTRAP, 100 * n_failed / N_BOOTSTRAP),
      "原因统计:",
      sprintf("  %s", paste0(names(reason_counts), ": ", as.integer(reason_counts)))
    )
    if (length(fail_msgs) > 0) {
      log_lines <- c(log_lines, "", "前20条错误信息:", sprintf("  %s", head(fail_msgs, 20)))
    }
    writeLines(log_lines, boot_log_path)
    cat(sprintf("  ✓ 失败日志已保存: %s\n", boot_log_path))
  }

  if (n_valid < N_BOOTSTRAP * 0.5) {
    stop("超过50%的Bootstrap样本拟合失败！模型可能存在问题。")
  }

  # ==================== 步骤8：汇总Bootstrap结果 ====================
  cat("\n=== 步骤8：汇总Bootstrap结果 ===\n")
  cat("计算95%置信区间...\n")

  # 提取所有系数
  all_coef_names <- unique(unlist(lapply(valid_boots, function(x) names(x$coefs))))
  boot_mat <- matrix(NA, nrow = n_valid, ncol = length(all_coef_names))
  colnames(boot_mat) <- all_coef_names

  for (i in seq_along(valid_boots)) {
    coefs <- valid_boots[[i]]$coefs
    for (j in seq_along(coefs)) {
      col_idx <- which(colnames(boot_mat) == names(coefs)[j])
      if (length(col_idx) > 0) {
        boot_mat[i, col_idx] <- coefs[j]
      }
    }
  }

  boot_filter <- filter_boot_matrix(boot_mat)
  boot_mat <- boot_filter$mat
  if (FILTER_SEM_OUTLIERS) {
    cat(sprintf("  输出异常值过滤（Bootstrap）：系数极端=%d\n",
                boot_filter$info["coef_extreme"]))
  }

  outlier_filter_df <- data.frame(
    metric = c("coef_extreme", "p_invalid", "r2_invalid", "boot_coef_extreme"),
    count = as.integer(c(param_filter$info["coef_extreme"],
                         param_filter$info["p_invalid"],
                         r2_filter$info["r2_invalid"],
                         boot_filter$info["coef_extreme"]))
  )
  write.csv(outlier_filter_df,
            file.path(OUTPUT_DIR, "SEM_pooled_outlier_filtering.csv"),
            row.names = FALSE)

  # 计算统计量
  boot_summary <- data.frame(
    Parameter = colnames(boot_mat),
    Mean = apply(boot_mat, 2, mean, na.rm = TRUE),
    SD = apply(boot_mat, 2, sd, na.rm = TRUE),
    CI_lower = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE),
    CI_upper = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE),
    Median = apply(boot_mat, 2, median, na.rm = TRUE)
  )

  # 保存Bootstrap结果
  write.csv(boot_summary, file.path(OUTPUT_DIR, "SEM_pooled_bootstrap_summary.csv"), row.names = FALSE)
  write.csv(boot_mat, file.path(OUTPUT_DIR, "SEM_pooled_bootstrap_coefficients.csv"), row.names = FALSE)

  cat("  ✓ Bootstrap结果已保存\n")

  # 打印关键结果
  cat("\n关键路径系数（Bootstrap均值 ± 95% CI）:\n")

  key_paths <- c("a1", "a2", "a3",  # 季前 → SOS
                 "b", "f1", "f2", "f3", "c1", "c2", "c3",  # → GPP
                 "g", "d", "h1", "h2", "h3", "e1", "e2", "e3")  # → Fixed_Trate

  for (path in key_paths) {
    if (path %in% boot_summary$Parameter) {
      row <- boot_summary[boot_summary$Parameter == path, ]
      cat(sprintf("  %s: %.3f [%.3f, %.3f]\n",
                  path, row$Mean, row$CI_lower, row$CI_upper))
    }
  }

  # ==================== 步骤9：生成报告 ====================
  cat("\n=== 步骤9：生成分析报告 ===\n")

  # 拟合优度判定（参考阈值）
  fit_ok <- is.finite(fit_measures["cfi"]) && fit_measures["cfi"] >= 0.95 &&
            is.finite(fit_measures["rmsea"]) && fit_measures["rmsea"] <= 0.06 &&
            is.finite(fit_measures["srmr"]) && fit_measures["srmr"] <= 0.08
  fit_quality_label <- if (fit_ok) {
    "模型拟合良好"
  } else {
    "模型拟合不佳（未达常用阈值，需谨慎解释）"
  }

  # 创建综合报告
  report <- c(
    "======================================================================",
    "稳健SEM分析报告 - 混合池方法",
    "======================================================================",
    "",
    sprintf("分析日期: %s", Sys.time()),
    sprintf("数据范围: %d-%d", YEAR_START, YEAR_END),
    sprintf("像元数量: %d", n_pix_total),
    sprintf("总观测数: %d", nrow(long_df_std)),
    sprintf("Bootstrap次数: %d (有效: %d)", N_BOOTSTRAP, n_valid),
    if (FILTER_SEM_OUTLIERS) {
      sprintf("异常值过滤统计: coef_extreme=%d, p_invalid=%d, r2_invalid=%d, boot_coef_extreme=%d",
              outlier_filter_df$count[outlier_filter_df$metric == "coef_extreme"],
              outlier_filter_df$count[outlier_filter_df$metric == "p_invalid"],
              outlier_filter_df$count[outlier_filter_df$metric == "r2_invalid"],
              outlier_filter_df$count[outlier_filter_df$metric == "boot_coef_extreme"])
    } else {
      "异常值过滤统计: 关闭"
    },
    "",
    "----------------------------------------------------------------------",
    "模型结构:",
    "----------------------------------------------------------------------",
    "第一层: 季前气候(P_pre, T_pre, SW_pre) → SOS",
    "第二层: SOS + 季前气候 + 生长季气候(P_season, T_season, SW_season) → GPP",
    "第三层: SOS + GPP + 季前气候 + 生长季气候 → Fixed_Trate",
    "",
    "----------------------------------------------------------------------",
    "模型拟合指标:",
    "----------------------------------------------------------------------",
    sprintf("χ² = %.2f, df = %d, p = %.4f",
            fit_measures["chisq"], fit_measures["df"], fit_measures["pvalue"]),
    sprintf("CFI = %.3f", fit_measures["cfi"]),
    sprintf("TLI = %.3f", fit_measures["tli"]),
    sprintf("RMSEA = %.3f [%.3f, %.3f]",
            fit_measures["rmsea"], fit_measures["rmsea.ci.lower"], fit_measures["rmsea.ci.upper"]),
    sprintf("SRMR = %.3f", fit_measures["srmr"]),
    "",
    "----------------------------------------------------------------------",
    "R²值:",
    "----------------------------------------------------------------------",
    sprintf("SOS: %.3f", r2_vals["SOS"]),
    sprintf("GPP_season: %.3f", r2_vals["GPP_season"]),
    sprintf("Fixed_Trate: %.3f", r2_vals["Fixed_Trate"]),
    "",
    "----------------------------------------------------------------------",
    "关键发现:",
    "----------------------------------------------------------------------",
    "1. 样本量极大，统计效力强",
    sprintf("   - 观测点数: %d", nrow(long_df_std)),
    sprintf("   - 相比Annual Mean样本量提升: %.0f倍", nrow(long_df_std) / 37),
    "",
    "2. Bootstrap提供稳健的置信区间",
    sprintf("   - 重采样次数: %d", n_valid),
    "   - 详见 SEM_pooled_bootstrap_summary.csv",
    "",
    sprintf("3. %s", fit_quality_label),
    sprintf("   - CFI = %.3f (参考阈值 >0.95)", fit_measures["cfi"]),
    sprintf("   - RMSEA = %.3f (参考阈值 <0.06)", fit_measures["rmsea"]),
    sprintf("   - SRMR = %.3f (参考阈值 <0.08)", fit_measures["srmr"]),
    "",
    "----------------------------------------------------------------------",
    "输出文件:",
    "----------------------------------------------------------------------",
    sprintf("数据目录: %s", DATA_DIR),
    "  1. sem_pooled_raw.csv - 原始混合数据",
    "  2. sem_pooled_standardized.csv - 标准化混合数据",
    "",
    sprintf("结果目录: %s", OUTPUT_DIR),
    "  3. SEM_pooled_fitsummary.txt - 完整拟合结果",
    "  4. SEM_pooled_parameters.csv - 所有参数估计",
    "  5. SEM_pooled_fitindices.csv - 拟合指标",
    "  6. SEM_pooled_R2.csv - R²值",
    "  7. SEM_pooled_bootstrap_summary.csv - Bootstrap结果摘要",
    "  8. SEM_pooled_bootstrap_coefficients.csv - 原始Bootstrap系数",
    "  9. SEM_pooled_analysis_report.txt - 本报告",
    "",
    "======================================================================",
    "分析完成！",
    "======================================================================"
  )

  writeLines(report, file.path(OUTPUT_DIR, "SEM_pooled_analysis_report.txt"))
  cat(report, sep = "\n")

  cat(sprintf("\n✓ 完整分析报告已保存: %s\n", file.path(OUTPUT_DIR, "SEM_pooled_analysis_report.txt")))

  cat("\n======================================================================\n")
  cat("✓✓✓ 所有分析完成！✓✓✓\n")
  cat(sprintf("输出目录: %s\n", OUTPUT_DIR))
  cat("======================================================================\n\n")
}

if (RUN_BOTH_DETREND) {
  run_list <- list(
    list(label = "原始结果", detrend = FALSE, suffix = ""),
    list(label = "去趋势结果", detrend = TRUE, suffix = "_detrended")
  )
} else {
  run_suffix <- if (DETREND_ENABLE) "_detrended" else ""
  run_label <- if (DETREND_ENABLE) "去趋势结果" else "原始结果"
  run_list <- list(list(label = run_label, detrend = DETREND_ENABLE, suffix = run_suffix))
}

for (run in run_list) {
  set_output_dirs(run$suffix)
  run_sem_pipeline(run$label, run$detrend)
}

set_output_dirs("")
