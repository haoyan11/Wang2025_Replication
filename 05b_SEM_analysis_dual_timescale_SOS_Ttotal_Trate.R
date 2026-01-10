#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Module 05b: 双时间尺度SEM分析 - 基于N04路径结构 + 04c固定窗口方法
#
# 设计思路：
# 1. 借鉴N04_Cal_SEM的双时间尺度路径结构
# 2. 使用04c代码的固定窗口方法（Fixed_Trate + 固定窗口[SOSav, POSav]）
# 3. 计算季前气候因子（SOSav前3个月）
# 4. 计算固定窗口生长季（SOSav-POSav）的GPP和气象因子
# 5. 构建完整路径：季前气候 → SOS → Fixed_Trate → Ttotal
# 6. 气候因子：降水(P)、气温(T)、短波辐射(SW)
#
# 路径结构（参考N04）：
#   季前气候(P,T,SW) → SOS → Fixed_Trate → Ttotal
#      ↓                ↓         ↑          ↑
#      └────────────────┴─────────┴──────────┘ (直接路径)
#   生长季气候(P,T,SW) ──────────┴──────────┘ (同期路径)
#
# ⚠️ 关键修改（对比原05b版本）：
# 1. Ttotal → TRc_{year}.tif（当年SOS-POS累计蒸腾）
# 2. Fixed_Trate = TR_fixed_window / Fixed_Window_Length（与04c一致）
# 3. 生长季气候：当年SOS-POS窗口 → 多年平均固定窗口[SOSav, POSav]
# 4. 生长季GPP：当年SOS-POS窗口 → 多年平均固定窗口[SOSav, POSav]
# 5. 季前气候：保持不变（基于SOSav前3个月）
#
# Version: 2.0.0 (04c Fixed Window Method)
# Author: Wang2025 Replication Project

suppressPackageStartupMessages({
  library(raster)
  library(lavaan)
  library(semPlot)
  library(parallel)  # 添加并行计算支持
})

# ===【并行化配置】===
# 并行工作进程数（固定10核）
N_CORES <- 10
cat(sprintf("\n[并行化配置] 使用 %d 个CPU核心进行并行计算\n", N_CORES))

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

# 自动生成的路径
OUTPUT_ROOT <- file.path(ROOT, "Wang2025_Analysis")
PHENO_DIR <- file.path(ROOT, "Phenology_Output_1", "GPP_phenology_EPSG4326")
# ===【修改1】Fixed_Trate数据源（04c固定窗口方法）===
# 原路径: Decomposition/TRproduct_{year}.tif（03a模块输出）
# 新路径: Decomposition_FixedWindow/TR_fixed_window_{year}.tif（03c模块输出）
# 说明：Fixed_Trate需临时计算 = TR_fixed_window / Fixed_Window_Length
DECOMP_DIR <- file.path(OUTPUT_ROOT, "Decomposition_FixedWindow")
TRC_DIR <- file.path(OUTPUT_ROOT, "TRc_annual")
MASK_FILE <- file.path(OUTPUT_ROOT, "masks", "combined_mask.tif")
TEMPLATE_FILE <- file.path(OUTPUT_ROOT, "masks", "template_grid.tif")

# 日尺度数据路径
GPP_DAILY_DIR <- file.path(ROOT, "GLASS_GPP", "GLASS_GPP_daily_interpolated")
PRECIP_DAILY_DIR <- file.path(ROOT, "Meteorological Data", "ERA5_Land", "Pre", "Pre_Daily", "Pre_Daily_2")
TA_DAILY_DIR <- file.path(ROOT, "Meteorological Data", "ERA5_Land", "Tem", "Tem_Daily", "Tem_Daily_2")
SW_DAILY_DIR <- file.path(ROOT, "Meteorological Data", "ERA5_Land", "DSW", "DSW_Daily", "DSW_Daily_2")

# 输出目录（与05a命名风格保持一致：SEM_Data_*/SEM_Results_*）
DATA_DIR <- file.path(OUTPUT_ROOT, "SEM_Data_Dual_Fixed_Ttotal_Trate")
DERIVED_DIR <- file.path(DATA_DIR, "Derived")
OUTPUT_DIR <- file.path(OUTPUT_ROOT, "SEM_Results_Dual_Fixed_Ttotal_Trate")
PIXELWISE_DIR <- file.path(OUTPUT_DIR, "Pixelwise")

# 像元模式复用已有物候气候态缓存
USE_EXISTING_CLIMATOLOGY <- TRUE

# 创建输出目录
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DERIVED_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PIXELWISE_DIR, showWarnings = FALSE, recursive = TRUE)

DATA_DIR_BASE <- DATA_DIR
DERIVED_DIR_BASE <- DERIVED_DIR
OUTPUT_DIR_BASE <- OUTPUT_DIR
PIXELWISE_DIR_BASE <- PIXELWISE_DIR

set_output_dirs <- function(suffix = "") {
  if (!is.null(suffix) && nzchar(suffix)) {
    DATA_DIR <<- paste0(DATA_DIR_BASE, suffix)
    OUTPUT_DIR <<- paste0(OUTPUT_DIR_BASE, suffix)
    PIXELWISE_DIR <<- file.path(OUTPUT_DIR, "Pixelwise")
  } else {
    DATA_DIR <<- DATA_DIR_BASE
    OUTPUT_DIR <<- OUTPUT_DIR_BASE
    PIXELWISE_DIR <<- PIXELWISE_DIR_BASE
  }
  DERIVED_DIR <<- DERIVED_DIR_BASE
  dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
  dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
  dir.create(PIXELWISE_DIR, showWarnings = FALSE, recursive = TRUE)
}

# 年份范围
YEAR_START <- 1982
YEAR_END <- 2018

# 全局常量
NODATA_OUT <- -9999
NODATA_ABS_MAX <- 1e20

# ==================== 控制选项 ====================
# SEM分析模式：
#   "pixel_time_series" - 每个像元的时间序列单独分析（计算量大，推荐用于最终结果）
#   "annual_mean"       - 全域年度均值合并分析（快速，适合测试）
#   "pixel_year"        - 像元-年份对合并分析（中等，适合大规模）
SEM_SAMPLE_MODE <- "pixel_time_series"  # 像元尺度先计算，再全域汇总

USE_CACHE <- TRUE                  # 启用缓存以加快处理速度
MIN_VALID_FRAC <- 0.60            # SOS-POS窗口内有效数据最低比例
MIN_VALID_YEAR_FRAC <- 0.60       # 像元时间序列最少有效年份比例
WRITE_PIXELWISE_RASTERS <- FALSE  # 是否输出像元级栅格（会生成大量文件）
DETREND_ENABLE <- FALSE           # 是否启用去趋势（仅对annual_mean生效）
DETREND_METHOD <- "linear"        # 目前仅支持linear
DETREND_PIXEL_ENABLE <- TRUE      # 是否启用像元级去趋势（仅pixel_time_series生效）
RUN_BOTH_DETREND <- TRUE          # annual_mean 同时输出不去趋势+去趋势
PIXEL_BLOCK_PARALLEL <- TRUE      # 是否启用像元块并行
PIXEL_BLOCK_CORES <- 10  # 像元块并行核心数
PIXEL_BLOCK_CHUNK <- 5            # 并行分块大小（块/批）
FILTER_SEM_OUTLIERS <- TRUE       # 输出结果异常值过滤（类似04c）
SEM_COEF_ABS_MAX <- 5             # 系数绝对值阈值（标准化系数）
SEM_P_MIN <- 0                    # p值下限
SEM_P_MAX <- 1                    # p值上限
SEM_R2_MIN <- 0                   # R²下限
SEM_R2_MAX <- 1                   # R²上限

# 日尺度文件命名模式
GPP_DAILY_PATTERN <- "GPP_{date}.tif"
PRECIP_DAILY_PATTERN <- "ERA5L_PrecipDaily_mm_{date}.tif"
TA_DAILY_PATTERN <- "ERA5L_T2mDaily_C_{date}.tif"
SW_DAILY_PATTERN <- "ERA5L_SWDaily_MJ_{date}.tif"

cat("\n======================================================================\n")
cat("双时间尺度SEM分析 - 固定窗口方法（04c Fixed Window）\n")
cat("======================================================================\n")

cat("\n[环境检测]\n")
cat(sprintf("  操作系统类型: %s\n", .Platform$OS.type))
cat(sprintf("  R版本: %s\n", R.version.string))
cat(sprintf("  根目录: %s\n", ROOT))

# ==================== 路径验证 ====================
cat("\n[路径检查]\n")
critical_paths <- list(
  "Root directory" = ROOT,
  "Phenology directory" = PHENO_DIR,
  "Decomposition FixedWindow directory" = DECOMP_DIR,  # 修改：使用03c输出目录
  "Mask file" = MASK_FILE,
  "GPP daily directory" = GPP_DAILY_DIR,
  "Precip daily directory" = PRECIP_DAILY_DIR,
  "Ta daily directory" = TA_DAILY_DIR,
  "SW daily directory" = SW_DAILY_DIR
)

all_exist <- TRUE
for (name in names(critical_paths)) {
  path <- critical_paths[[name]]
  exists <- file.exists(path) || dir.exists(path)
  status <- ifelse(exists, "[OK]", "[MISSING]")
  cat(sprintf("  %s %s\n    %s\n", status, name, path))
  if (!exists) all_exist <- FALSE
}

if (!all_exist) {
  cat("\n⚠️ 警告：部分路径不存在！\n")
  cat("  请检查：\n")
  cat("  1. ROOT路径是否正确\n")
  cat("  2. 是否已运行01-03模块生成必需数据\n")
  cat("  3. 路径分隔符是否正确（Windows用反斜杠\\或正斜杠/）\n\n")
  stop("路径验证失败，请修正后重试")
}
cat("✓ 所有关键路径存在\n")

cat(sprintf("\n年份范围: %d-%d\n", YEAR_START, YEAR_END))
cat(sprintf("分析模式: %s\n", SEM_SAMPLE_MODE))
cat("----------------------------------------------------------------------\n")

# ==================== 辅助函数 ====================

# 生成年份内的日期序列
date_seq_year <- function(year) {
  seq.Date(as.Date(sprintf("%d-01-01", year)), as.Date(sprintf("%d-12-31", year)), by = "day")
}

# 构建日尺度文件路径
build_daily_path <- function(dir_path, pattern, date_str) {
  file.path(dir_path, gsub("\\{date\\}", date_str, pattern))
}

# 诊断函数：检查日尺度文件
diagnose_daily_files <- function(year, dir_path, pattern, sample_dates = 5) {
  cat(sprintf("\n  诊断: %s\n", dir_path))
  cat(sprintf("    模式: %s\n", pattern))

  dates <- date_seq_year(year)
  sample_idx <- seq(1, length(dates), length.out = min(sample_dates, length(dates)))

  for (i in sample_idx) {
    d <- dates[i]
    date_str <- strftime(d, "%Y%m%d")
    f <- build_daily_path(dir_path, pattern, date_str)
    exists <- file.exists(f)
    status <- ifelse(exists, "✓", "✗")
    cat(sprintf("    %s %s (DOY %03d)\n", status, basename(f), as.integer(strftime(d, "%j"))))
    if (!exists && i == sample_idx[1]) {
      # 显示完整路径帮助调试
      cat(sprintf("      完整路径: %s\n", f))
    }
  }
}

# 数据清洗：移除无效值
sanitize_values <- function(vals, nodata, allow_negative = TRUE) {
  # 首先处理非有限值（NA, Inf, -Inf, NaN）
  valid <- is.finite(vals) & (abs(vals) < NODATA_ABS_MAX)

  # 处理nodata值
  if (!is.null(nodata) && !is.na(nodata) && is.finite(nodata)) {
    valid <- valid & (vals != nodata)
  } else {
    valid <- valid & (vals != NODATA_OUT)
  }

  # 处理负值
  if (!allow_negative) {
    valid <- valid & (vals >= 0)
  }

  # 清理无效值
  vals[!valid] <- NA_real_

  # 额外检查：确保没有Inf/-Inf残留
  vals[is.infinite(vals)] <- NA_real_

  vals
}

set_nodata_if_missing <- function(r, fallback = NODATA_OUT) {
  nodata <- NAvalue(r)
  if (is.na(nodata) || !is.finite(nodata)) {
    NAvalue(r) <- fallback
  }
  r
}

filter_sem_outputs <- function(vals, p_vals = NULL, r2_vals = NULL) {
  info <- c(coef_extreme = 0, p_invalid = 0, r2_invalid = 0)
  if (!FILTER_SEM_OUTLIERS) {
    return(list(vals = vals, p_vals = p_vals, r2_vals = r2_vals, info = info))
  }

  if (!is.null(vals)) {
    bad_coef <- is.finite(vals) & (abs(vals) > SEM_COEF_ABS_MAX)
    info["coef_extreme"] <- sum(bad_coef)
    vals[bad_coef] <- NA_real_
  }

  if (!is.null(p_vals)) {
    bad_p <- is.finite(p_vals) & (p_vals < SEM_P_MIN | p_vals > SEM_P_MAX)
    info["p_invalid"] <- sum(bad_p)
    p_vals[bad_p] <- NA_real_
  }

  if (!is.null(r2_vals)) {
    bad_r2 <- is.finite(r2_vals) & (r2_vals < SEM_R2_MIN | r2_vals > SEM_R2_MAX)
    info["r2_invalid"] <- sum(bad_r2)
    r2_vals[bad_r2] <- NA_real_
  }

  list(vals = vals, p_vals = p_vals, r2_vals = r2_vals, info = info)
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

# ==================== 去趋势工具 ====================
detrend_series <- function(x, t) {
  ok <- is.finite(x) & is.finite(t)
  if (sum(ok) < 3) {
    return(x)
  }
  fit <- lm(x ~ t, subset = ok)
  res <- rep(NA_real_, length(x))
  res[ok] <- residuals(fit)
  res
}

apply_detrend_sem_data <- function(df, year_col = "year") {
  out <- df
  t <- out[[year_col]]
  for (col in names(out)) {
    if (col == year_col) next
    out[[col]] <- detrend_series(out[[col]], t)
  }
  out
}

# Z-score标准化
scale_vec <- function(x) {
  # 移除NA值计算均值和标准差
  x_clean <- x[is.finite(x)]
  if (length(x_clean) < 2) {
    return(rep(NA_real_, length(x)))
  }
  m <- mean(x_clean, na.rm = TRUE)
  s <- sd(x_clean, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(rep(NA_real_, length(x)))
  }
  (x - m) / s
}

# 计算VIF（方差膨胀因子）
calc_vif_vector <- function(X) {
  p <- ncol(X)
  vif_vals <- rep(NA_real_, p)
  for (i in seq_len(p)) {
    y <- X[, i]
    X_others <- X[, -i, drop = FALSE]
    if (ncol(X_others) == 0) {
      vif_vals[i] <- 1
      next
    }
    beta <- tryCatch(qr.solve(X_others, y), error = function(e) rep(NA_real_, ncol(X_others)))
    if (any(is.na(beta))) {
      vif_vals[i] <- NA_real_
      next
    }
    y_pred <- X_others %*% beta
    rss <- sum((y - y_pred)^2)
    tss <- sum((y - mean(y))^2)
    r2 <- ifelse(tss > 0, 1 - rss / tss, 0)
    vif_vals[i] <- ifelse(r2 < 0.9999, 1 / (1 - r2), Inf)
  }
  vif_vals
}

# 线性回归：计算系数和p值
regress_beta_p <- function(X, y) {
  XtX <- crossprod(X)
  if (qr(XtX)$rank < ncol(XtX)) {
    return(list(beta = rep(NA_real_, ncol(X)), p = rep(NA_real_, ncol(X))))
  }
  beta <- solve(XtX, crossprod(X, y))
  resid <- y - X %*% beta
  df <- nrow(X) - ncol(X)
  if (df <= 0) {
    return(list(beta = rep(NA_real_, ncol(X)), p = rep(NA_real_, ncol(X))))
  }
  sigma2 <- sum(resid^2) / df
  se <- sqrt(diag(sigma2 * solve(XtX)))
  tval <- beta / se
  pval <- 2 * (1 - pt(abs(tval), df))
  list(beta = as.vector(beta), p = as.vector(pval))
}

# 计算R2（不含截距）
calc_r2 <- function(X, y) {
  if (is.null(dim(X)) || ncol(X) == 0) {
    return(NA_real_)
  }
  beta <- tryCatch(qr.solve(X, y), error = function(e) rep(NA_real_, ncol(X)))
  if (any(is.na(beta))) {
    return(NA_real_)
  }
  y_pred <- X %*% beta
  rss <- sum((y - y_pred)^2)
  tss <- sum((y - mean(y))^2)
  if (!is.finite(tss) || tss <= 0) {
    return(NA_real_)
  }
  1 - rss / tss
}

# 相关系数和p值
corr_with_p <- function(x, y) {
  r <- cor(x, y)
  n <- length(x)
  if (!is.finite(r) || n < 3) {
    return(list(r = NA_real_, p = NA_real_))
  }
  r <- max(min(r, 0.999999), -0.999999)
  tval <- r * sqrt((n - 2) / (1 - r^2))
  pval <- 2 * (1 - pt(abs(tval), n - 2))
  list(r = r, p = pval)
}

# 栅格掩膜
mask_raster <- function(r, mask_r) {
  mask(r, mask_r)
}

check_raster_alignment <- function(template_r, mask_r, year) {
  if (!compareRaster(template_r, mask_r, stopiffalse = FALSE)) {
    stop("模板与掩膜网格不一致，请检查 template_grid.tif / combined_mask.tif")
  }

  sample_files <- list(
    list(file.path(DECOMP_DIR, sprintf("TR_fixed_window_%d.tif", year)), "TR_fixed_window"),
    list(file.path(DECOMP_DIR, "Fixed_Window_Length.tif"), "Fixed_Window_Length"),
    list(file.path(PHENO_DIR, sprintf("sos_gpp_%d.tif", year)), "SOS"),
    list(file.path(PHENO_DIR, sprintf("pos_doy_gpp_%d.tif", year)), "POS"),
    list(file.path(TRC_DIR, sprintf("TRc_%d.tif", year)), "TRc")
  )

  for (item in sample_files) {
    f <- item[[1]]
    label <- item[[2]]
    if (!file.exists(f)) {
      cat(sprintf("  ⚠ 缺少样本文件: %s (%s)\n", label, basename(f)))
      next
    }
    ok <- compareRaster(template_r, raster(f), stopiffalse = FALSE)
    if (!ok) {
      stop(sprintf("网格不一致: %s (%s)", label, basename(f)))
    } else {
      cat(sprintf("  ✓ 网格一致: %s (%s)\n", label, basename(f)))
    }
  }
}

# 读取栅格块
read_block_values <- function(file_path, row, nrows, allow_negative = TRUE) {
  r <- raster(file_path)
  vals <- getValues(r, row = row, nrows = nrows)
  nodata <- NAvalue(r)
  sanitize_values(vals, nodata, allow_negative)
}

# 闰年判断 + 无闰日DOY映射到真实日期
is_leap_year <- function(year) {
  (year %% 4 == 0 & year %% 100 != 0) | (year %% 400 == 0)
}

doy_to_date_noleap <- function(year, doy) {
  if (is.na(doy) || doy < 1 || doy > 365) {
    return(NA)
  }
  offset <- if (is_leap_year(year) && doy >= 60) 1 else 0
  as.Date(paste0(year, "-01-01")) + (doy - 1 + offset)
}

# ==================== 核心计算函数 ====================

# 计算季前气候因子（基于像元多年平均SOS）
calc_preseason_climate <- function(year, sos_climatology_r, daily_dir, pattern, mask_r,
                                    cache_file, var_name, allow_negative = TRUE) {
  if (file.exists(cache_file)) {
    cat(sprintf("    [缓存] %s\n", var_name))
    return(raster(cache_file))
  }

  cat(sprintf("    计算 %s (基于像元多年平均SOS前3个月)...\n", var_name))

  template <- mask_r
  sos_clim_vals <- getValues(sos_climatology_r)
  sos_clim_vals <- sanitize_values(sos_clim_vals, NAvalue(sos_climatology_r), allow_negative = FALSE)

  # 初始化输出数组
  mean_vals <- rep(NA_real_, ncell(template))

  # 确定有效像元
  valid_cells <- which(!is.na(sos_clim_vals))
  n_valid <- length(valid_cells)

  cat(sprintf("      有效像元数: %d\n", n_valid))

  # 统计季前窗口分布
  sos_clim_valid <- sos_clim_vals[valid_cells]
  cat(sprintf("      SOS气候态范围: %.0f - %.0f DOY\n",
              min(sos_clim_valid, na.rm = TRUE),
              max(sos_clim_valid, na.rm = TRUE)))

  # ===【优化：使用缓存避免重复读取文件】===
  # 逐像元计算（使用日期缓存大幅提升效率）
  sos_int <- round(sos_clim_valid)
  sos_int <- pmax(1, pmin(365, sos_int))  # 修复BUG：SOS下限应为1而非90（季前窗口可跨年到前一年）

  # 确定需要读取的日期范围
  min_sos <- min(sos_int, na.rm = TRUE)
  max_sos <- max(sos_int, na.rm = TRUE)
  earliest_start <- min_sos - 90

  # 构建需要读取的DOY列表
  if (earliest_start < 1) {
    doys_prev_year <- (365 + earliest_start):365
    doys_this_year <- 1:(max_sos - 1)
  } else {
    doys_prev_year <- integer(0)
    doys_this_year <- earliest_start:(max_sos - 1)
  }

  total_days <- length(doys_prev_year) + length(doys_this_year)
  cat(sprintf("      需读取日期范围: %d 天\n", total_days))

  # 创建日期数据缓存（环境变量）
  daily_cache <- new.env(hash = TRUE, size = total_days)

  # 预读取所有需要的日期数据（每个文件只读一次）
  cat("      [1/2] 正在预读取日期数据...\n")
  files_read <- 0

  # 读取上年数据
  if (length(doys_prev_year) > 0) {
    for (doy in doys_prev_year) {
      date_obj <- doy_to_date_noleap(year - 1, doy)
      if (is.na(date_obj)) {
        next
      }
      date_str <- format(date_obj, "%Y%m%d")
      file_path <- build_daily_path(daily_dir, pattern, date_str)

      if (file.exists(file_path)) {
        tryCatch({
          daily_r <- raster(file_path)
          daily_vals <- getValues(daily_r)
          daily_vals <- sanitize_values(daily_vals, NAvalue(daily_r), allow_negative)
          cache_key <- paste0("doy_", doy, "_year_", year - 1)
          assign(cache_key, daily_vals, envir = daily_cache)
          files_read <- files_read + 1
        }, error = function(e) {})
      }
    }
  }

  # 读取本年数据
  for (doy in doys_this_year) {
    date_obj <- doy_to_date_noleap(year, doy)
    if (is.na(date_obj)) {
      next
    }
    date_str <- format(date_obj, "%Y%m%d")
    file_path <- build_daily_path(daily_dir, pattern, date_str)

    if (file.exists(file_path)) {
      tryCatch({
        daily_r <- raster(file_path)
        daily_vals <- getValues(daily_r)
        daily_vals <- sanitize_values(daily_vals, NAvalue(daily_r), allow_negative)
        cache_key <- paste0("doy_", doy, "_year_", year)
        assign(cache_key, daily_vals, envir = daily_cache)
        files_read <- files_read + 1
      }, error = function(e) {})
    }
  }

  cat(sprintf("      ✓ 已缓存 %d 个日期文件\n", files_read))

  # 逐像元计算季前均值（从缓存读取，极快）
  cat("      [2/2] 正在计算像元季前均值...\n")

  for (i in seq_along(valid_cells)) {
    cell_idx <- valid_cells[i]
    sos_doy <- sos_int[i]

    # 确定季前窗口
    preseason_start_doy <- sos_doy - 90
    preseason_end_doy <- sos_doy - 1

    # 构建窗口DOY列表
    if (preseason_start_doy < 1) {
      window_doys <- c((365 + preseason_start_doy):365, 1:preseason_end_doy)
      window_years <- c(rep(year - 1, 365 - (365 + preseason_start_doy) + 1),
                       rep(year, preseason_end_doy))
    } else {
      window_doys <- preseason_start_doy:preseason_end_doy
      window_years <- rep(year, length(window_doys))
    }

    # 从缓存提取该像元在窗口内的所有值
    cell_vals <- numeric(0)
    for (j in seq_along(window_doys)) {
      cache_key <- paste0("doy_", window_doys[j], "_year_", window_years[j])
      if (exists(cache_key, envir = daily_cache)) {
        daily_vals <- get(cache_key, envir = daily_cache)
        if (cell_idx <= length(daily_vals)) {
          val <- daily_vals[cell_idx]
          if (!is.na(val)) {
            cell_vals <- c(cell_vals, val)
          }
        }
      }
    }

    # 计算均值
    if (length(cell_vals) > 0) {
      mean_vals[cell_idx] <- mean(cell_vals)
    }

    # 进度显示
    if (i %% max(1, floor(n_valid / 5)) == 0 || i == n_valid) {
      cat(sprintf("      进度: %d/%d 像元 (%.1f%%)\n",
                  i, n_valid, 100 * i / n_valid))
    }
  }

  cat(sprintf("      ✓ 完成: %d 像元计算完毕\n", sum(!is.na(mean_vals))))

  # 构建输出栅格
  out_r <- setValues(template, mean_vals)
  out_r <- mask_raster(out_r, mask_r)

  # 保存缓存
  writeRaster(out_r, cache_file, overwrite = TRUE, datatype = "FLT4S")

  out_r
}

# ===【修改2】计算生长季气候因子（固定窗口[SOSav, POSav]均值）===
# 原版本：使用当年SOS-POS窗口
# 新版本：使用多年平均固定窗口[SOSav, POSav]
calc_season_climate_fixed <- function(year, sos_climatology_r, pos_climatology_r, daily_dir, pattern, mask_r,
                                      cache_file, var_name, allow_negative = TRUE) {
  if (file.exists(cache_file)) {
    cat(sprintf("    [缓存] %s\n", var_name))
    return(raster(cache_file))
  }

  cat(sprintf("    计算 %s (固定窗口[SOSav, POSav]均值)...\n", var_name))

  template <- mask_r
  # 使用多年平均物候（气候态）
  sos_vals <- getValues(sos_climatology_r)
  pos_vals <- getValues(pos_climatology_r)

  sos_vals <- sanitize_values(sos_vals, NAvalue(sos_climatology_r), allow_negative = FALSE)
  pos_vals <- sanitize_values(pos_vals, NAvalue(pos_climatology_r), allow_negative = FALSE)

  sos_int <- as.integer(round(sos_vals))
  pos_int <- as.integer(round(pos_vals))

  # 限制到1-365
  sos_int <- pmax(1, pmin(365, sos_int))
  pos_int <- pmax(1, pmin(365, pos_int))

  valid_pheno <- !is.na(sos_int) & !is.na(pos_int) & (sos_int < pos_int)  # 修复：改为严格小于，与POSav>SOSav验证保持一致

  # 初始化
  sum_vals <- rep(0, ncell(template))
  cnt_vals <- rep(0L, ncell(template))

  if (!any(valid_pheno)) {
    out_r <- setValues(template, rep(NA_real_, ncell(template)))
    out_r <- mask_raster(out_r, mask_r)
    writeRaster(out_r, cache_file, overwrite = TRUE)
    return(out_r)
  }

  # 确定DOY范围
  min_doy <- min(sos_int[valid_pheno], na.rm = TRUE)
  max_doy <- max(pos_int[valid_pheno], na.rm = TRUE)

  cat(sprintf("      物候DOY范围: %d-%d\n", min_doy, max_doy))

  # 逐日累加（使用无闰日DOY映射到真实日期）
  for (doy in seq(min_doy, max_doy)) {
    date_obj <- doy_to_date_noleap(year, doy)
    if (is.na(date_obj)) {
      next
    }
    date_str <- format(date_obj, "%Y%m%d")

    file_path <- build_daily_path(daily_dir, pattern, date_str)

    if (!file.exists(file_path)) next

    daily_result <- tryCatch({
      daily_r <- raster(file_path)
      daily_vals <- getValues(daily_r)
      daily_vals <- sanitize_values(daily_vals, NAvalue(daily_r), allow_negative)

      valid_daily <- !is.na(daily_vals)
      in_window <- valid_pheno & (sos_int <= doy) & (pos_int >= doy)
      use_mask <- in_window & valid_daily

      list(daily_vals = daily_vals, use_mask = use_mask)
    }, error = function(e) {
      NULL
    })

    if (!is.null(daily_result) && any(daily_result$use_mask)) {
      sum_vals[daily_result$use_mask] <- sum_vals[daily_result$use_mask] +
        daily_result$daily_vals[daily_result$use_mask]
      cnt_vals[daily_result$use_mask] <- cnt_vals[daily_result$use_mask] + 1L
    }
  }

  window_len <- pos_int - sos_int + 1
  min_required <- ceiling(MIN_VALID_FRAC * window_len)
  cat(sprintf("      有效像元: %d (cnt >= %.0f%%窗口长度)\n",
              sum(cnt_vals >= min_required, na.rm = TRUE),
              MIN_VALID_FRAC * 100))

  # 计算均值
  mean_vals <- rep(NA_real_, ncell(template))
  good <- valid_pheno & (cnt_vals >= min_required)
  mean_vals[good] <- sum_vals[good] / cnt_vals[good]

  out_r <- setValues(template, mean_vals)
  out_r <- mask_raster(out_r, mask_r)

  writeRaster(out_r, cache_file, overwrite = TRUE, datatype = "FLT4S")

  out_r
}

# ==================== 数据准备 ====================
# ===【修改3】数据准备函数：Fixed_Trate + 固定窗口===
prepare_dual_timescale_data <- function(year, sos_climatology_r, pos_climatology_r,
                                       fixed_window_length_r, mask_r) {
  cat(sprintf("\n年份: %d\n", year))

  # ===【修改3a】读取TR_fixed_window并计算Fixed_Trate（替代TRc）===
  tr_fixed_file <- file.path(DECOMP_DIR, sprintf("TR_fixed_window_%d.tif", year))
  if (!file.exists(tr_fixed_file)) {
    cat(sprintf("  跳过: TR_fixed_window文件不存在 (%s)\n", basename(tr_fixed_file)))
    return(NULL)
  }

  # 读取当年物候数据（用于SEM中的SOS变量）
  sos_file <- file.path(PHENO_DIR, sprintf("sos_gpp_%d.tif", year))

  if (!file.exists(sos_file)) {
    cat(sprintf("  跳过: SOS文件不存在\n"))
    return(NULL)
  }

  tr_fixed_r <- raster(tr_fixed_file)
  ttotal_file <- file.path(TRC_DIR, sprintf("TRc_%d.tif", year))
  if (!file.exists(ttotal_file)) {
    cat(sprintf("  ✗ Ttotal文件不存在: %s\n", basename(ttotal_file)))
    return(NULL)
  }
  ttotal_r <- raster(ttotal_file)
  sos_r <- raster(sos_file)

  tr_fixed_r <- set_nodata_if_missing(tr_fixed_r)
  ttotal_r <- set_nodata_if_missing(ttotal_r)
  sos_r <- set_nodata_if_missing(sos_r)

  # 应用掩膜
  sos_r <- mask_raster(sos_r, mask_r)
  tr_fixed_r <- mask_raster(tr_fixed_r, mask_r)
  ttotal_r <- mask_raster(ttotal_r, mask_r)

  # 计算Fixed_Trate = TR_fixed_window / Fixed_Window_Length
  cat("  [计算Fixed_Trate] TR_fixed_window / Fixed_Window_Length...\n")
  fixed_trate_r <- tr_fixed_r / fixed_window_length_r
  fixed_trate_r <- mask_raster(fixed_trate_r, mask_r)

  # 缓存文件路径
  cache_prefix <- file.path(DERIVED_DIR, sprintf("%s_%d.tif", c(
    "P_pre", "T_pre", "SW_pre",
    "P_season", "T_season", "SW_season",
    "GPP_season"
  ), year))

  names(cache_prefix) <- c("P_pre", "T_pre", "SW_pre",
                           "P_season", "T_season", "SW_season",
                           "GPP_season")

  # ===== 计算季前气候因子（使用SOS气候态） =====
  cat("  [1/7] 季前气候因子（基于像元多年平均SOS）:\n")

  # ===【并行化优化】季前3个变量并行计算===
  if (N_CORES > 1) {
    cat(sprintf("    使用 %d 个核心并行计算...\n", N_CORES))

    # 创建并行集群（Windows兼容）
    cl <- makeCluster(N_CORES)
    on.exit(stopCluster(cl), add = TRUE)

    # 导出必要的变量和函数到集群
    clusterExport(cl, c("year", "sos_climatology_r", "mask_r", "cache_prefix",
                        "calc_preseason_climate", "build_daily_path", "sanitize_values",
                        "mask_raster", "NODATA_OUT", "NODATA_ABS_MAX",
                        "MIN_VALID_FRAC", "is_leap_year", "doy_to_date_noleap",
                        "PRECIP_DAILY_DIR", "PRECIP_DAILY_PATTERN",
                        "TA_DAILY_DIR", "TA_DAILY_PATTERN",
                        "SW_DAILY_DIR", "SW_DAILY_PATTERN"),
                  envir = environment())

    # 加载raster包到每个工作进程
    clusterEvalQ(cl, library(raster))

    # 并行计算3个季前变量
    preseason_list <- parLapply(cl, list(
      list(dir = PRECIP_DAILY_DIR, pattern = PRECIP_DAILY_PATTERN, var = "P_pre", neg = FALSE),
      list(dir = TA_DAILY_DIR, pattern = TA_DAILY_PATTERN, var = "T_pre", neg = TRUE),
      list(dir = SW_DAILY_DIR, pattern = SW_DAILY_PATTERN, var = "SW_pre", neg = FALSE)
    ), function(params) {
      calc_preseason_climate(year, sos_climatology_r, params$dir, params$pattern,
                            mask_r, cache_prefix[params$var], params$var, params$neg)
    })

    stopCluster(cl)
    on.exit()  # 清除on.exit

    p_pre_r <- preseason_list[[1]]
    t_pre_r <- preseason_list[[2]]
    sw_pre_r <- preseason_list[[3]]

  } else {
    # 单核模式（原始串行计算）
    p_pre_r <- calc_preseason_climate(year, sos_climatology_r, PRECIP_DAILY_DIR, PRECIP_DAILY_PATTERN,
                                       mask_r, cache_prefix["P_pre"], "P_pre", allow_negative = FALSE)

    t_pre_r <- calc_preseason_climate(year, sos_climatology_r, TA_DAILY_DIR, TA_DAILY_PATTERN,
                                       mask_r, cache_prefix["T_pre"], "T_pre", allow_negative = TRUE)

    sw_pre_r <- calc_preseason_climate(year, sos_climatology_r, SW_DAILY_DIR, SW_DAILY_PATTERN,
                                        mask_r, cache_prefix["SW_pre"], "SW_pre", allow_negative = FALSE)
  }

  # ===【修改3b】计算固定窗口生长季气候因子 =====
  cat("  [2/7] 生长季气候因子（固定窗口[SOSav, POSav]）:\n")
  cat("  [3/7] 生长季GPP（固定窗口[SOSav, POSav]）:\n")

  # ===【并行化优化】生长季4个变量并行计算===
  if (N_CORES > 1) {
    cat(sprintf("    使用 %d 个核心并行计算...\n", N_CORES))

    # 创建并行集群（Windows兼容）
    cl <- makeCluster(min(N_CORES, 4))  # 最多4个变量
    on.exit(stopCluster(cl), add = TRUE)

    # 导出必要的变量和函数到集群
    clusterExport(cl, c("year", "sos_climatology_r", "pos_climatology_r", "mask_r", "cache_prefix",
                        "calc_season_climate_fixed", "build_daily_path", "sanitize_values",
                        "mask_raster", "NODATA_OUT", "NODATA_ABS_MAX",
                        "MIN_VALID_FRAC", "is_leap_year", "doy_to_date_noleap",
                        "PRECIP_DAILY_DIR", "PRECIP_DAILY_PATTERN",
                        "TA_DAILY_DIR", "TA_DAILY_PATTERN",
                        "SW_DAILY_DIR", "SW_DAILY_PATTERN",
                        "GPP_DAILY_DIR", "GPP_DAILY_PATTERN"),
                  envir = environment())

    # 加载raster包到每个工作进程
    clusterEvalQ(cl, library(raster))

    # 并行计算4个生长季变量（3个气候+1个GPP）
    season_list <- parLapply(cl, list(
      list(dir = PRECIP_DAILY_DIR, pattern = PRECIP_DAILY_PATTERN, var = "P_season", neg = FALSE),
      list(dir = TA_DAILY_DIR, pattern = TA_DAILY_PATTERN, var = "T_season", neg = TRUE),
      list(dir = SW_DAILY_DIR, pattern = SW_DAILY_PATTERN, var = "SW_season", neg = FALSE),
      list(dir = GPP_DAILY_DIR, pattern = GPP_DAILY_PATTERN, var = "GPP_season", neg = FALSE)
    ), function(params) {
      calc_season_climate_fixed(year, sos_climatology_r, pos_climatology_r, params$dir, params$pattern,
                                mask_r, cache_prefix[params$var], params$var, params$neg)
    })

    stopCluster(cl)
    on.exit()  # 清除on.exit

    p_season_r <- season_list[[1]]
    t_season_r <- season_list[[2]]
    sw_season_r <- season_list[[3]]
    gpp_season_r <- season_list[[4]]

  } else {
    # 单核模式（原始串行计算）
    p_season_r <- calc_season_climate_fixed(year, sos_climatology_r, pos_climatology_r, PRECIP_DAILY_DIR, PRECIP_DAILY_PATTERN,
                                             mask_r, cache_prefix["P_season"], "P_season", allow_negative = FALSE)

    t_season_r <- calc_season_climate_fixed(year, sos_climatology_r, pos_climatology_r, TA_DAILY_DIR, TA_DAILY_PATTERN,
                                             mask_r, cache_prefix["T_season"], "T_season", allow_negative = TRUE)

    sw_season_r <- calc_season_climate_fixed(year, sos_climatology_r, pos_climatology_r, SW_DAILY_DIR, SW_DAILY_PATTERN,
                                              mask_r, cache_prefix["SW_season"], "SW_season", allow_negative = FALSE)

    gpp_season_r <- calc_season_climate_fixed(year, sos_climatology_r, pos_climatology_r, GPP_DAILY_DIR, GPP_DAILY_PATTERN,
                                               mask_r, cache_prefix["GPP_season"], "GPP_season", allow_negative = FALSE)
  }

  cat("  ✓ 数据准备完成\n")

  # ===【修改3d】返回Fixed_Trate（替代TRc）===
  list(
    ttotal = ttotal_r,
    fixed_trate = fixed_trate_r,
    sos = sos_r,
    p_pre = p_pre_r,
    t_pre = t_pre_r,
    sw_pre = sw_pre_r,
    p_season = p_season_r,
    t_season = t_season_r,
    sw_season = sw_season_r,
    gpp_season = gpp_season_r
  )
}

# ==================== 缓存准备 ====================
# ===【修改4】添加POSav和fixed_window_length_r参数===
prepare_sem_caches <- function(years, sos_climatology_r, pos_climatology_r,
                               fixed_window_length_r, mask_r) {
  cat("\n=== 准备SEM缓存数据 ===\n")
  for (year in years) {
    cat(sprintf("  年份: %d\n", year))
    tryCatch({
      prepare_dual_timescale_data(year, sos_climatology_r, pos_climatology_r,
                                 fixed_window_length_r, mask_r)
      cat(sprintf("    ✓ 完成\n"))
    }, error = function(e) {
      cat(sprintf("    ✗ 错误: %s\n", conditionMessage(e)))
      cat(sprintf("    提示: 检查年份 %d 的输入数据是否完整\n", year))
    })
  }
}

# ==================== VIF诊断 ====================
# ===【修改5】Fixed_Trate替代TRc===
calculate_vif <- function(data) {
  # 双时间尺度模型的VIF诊断
  lm_model <- lm(Ttotal ~ SOS + Fixed_Trate + P_pre + T_pre + SW_pre +
                   P_season + T_season + SW_season, data = data)
  x <- model.matrix(lm_model)[, -1, drop = FALSE]

  vif_vals <- rep(NA_real_, ncol(x))
  names(vif_vals) <- colnames(x)

  for (i in seq_len(ncol(x))) {
    y_i <- x[, i]
    x_others <- x[, -i, drop = FALSE]
    fit <- lm.fit(cbind(1, x_others), y_i)
    rss <- sum(fit$residuals^2)
    tss <- sum((y_i - mean(y_i))^2)
    r2 <- ifelse(tss > 0, 1 - rss / tss, 0)
    vif_vals[i] <- ifelse(r2 < 0.9999, 1 / (1 - r2), Inf)
  }

  vif_df <- data.frame(Variable = names(vif_vals), VIF = vif_vals)
  write.csv(vif_df, file.path(OUTPUT_DIR, "VIF_diagnostics.csv"), row.names = FALSE)

  cat("\nVIF诊断结果:\n")
  print(vif_df)
  cat("\n说明：VIF > 10 表示存在严重多重共线性\n")

  vif_vals
}

# ==================== SEM模型（像元时间序列）====================
# ===【修改6】像元时间序列SEM：Fixed_Trate替代TRc===
run_sem_pixel_time_series <- function(years, sos_climatology_r, fixed_window_length_r, mask_r) {
  cat("\n=== SEM分析（像元时间序列）===\n")
  cat("⚠️ 注意：此模式计算量极大，可能需要数小时！\n")
  cat("  建议：先用annual_mean模式测试\n\n")

  # 构建文件列表
  files <- list(
    TR_fixed_window = file.path(DECOMP_DIR, sprintf("TR_fixed_window_%d.tif", years)),
    TRc_annual = file.path(TRC_DIR, sprintf("TRc_%d.tif", years)),
    SOS = file.path(PHENO_DIR, sprintf("sos_gpp_%d.tif", years)),
    GPP_season = file.path(DERIVED_DIR, sprintf("GPP_season_%d.tif", years)),
    P_pre = file.path(DERIVED_DIR, sprintf("P_pre_%d.tif", years)),
    T_pre = file.path(DERIVED_DIR, sprintf("T_pre_%d.tif", years)),
    SW_pre = file.path(DERIVED_DIR, sprintf("SW_pre_%d.tif", years)),
    P_season = file.path(DERIVED_DIR, sprintf("P_season_%d.tif", years)),
    T_season = file.path(DERIVED_DIR, sprintf("T_season_%d.tif", years)),
    SW_season = file.path(DERIVED_DIR, sprintf("SW_season_%d.tif", years))
  )

  keep <- rep(TRUE, length(years))
  for (key in names(files)) {
    keep <- keep & file.exists(files[[key]])
  }

  years <- years[keep]
  for (key in names(files)) {
    files[[key]] <- files[[key]][keep]
  }

  min_years <- max(1, ceiling(MIN_VALID_YEAR_FRAC * length(years)))
  if (length(years) < min_years) {
    stop("像元时间序列SEM所需的完整输入年份不足")
  }

  stacks <- lapply(files, stack)
  na_values <- lapply(files, function(paths) NAvalue(raster(paths[1])))

  template <- raster(files$TR_fixed_window[1])  # 修改：TRc → TR_fixed_window
  bs <- blockSize(template)
  mask_vals_all <- getValues(mask_r)
  total_valid_cells <- sum(!is.na(mask_vals_all))
  processed_cells <- 0

  # 路径名称（双时间尺度）
  coef_names <- c(
    "SOS~P_pre",
    "SOS~T_pre",
    "SOS~SW_pre",
    "Fixed_Trate~SOS",
    "Fixed_Trate~P_pre",
    "Fixed_Trate~T_pre",
    "Fixed_Trate~SW_pre",
    "Fixed_Trate~P_season",
    "Fixed_Trate~T_season",
    "Fixed_Trate~SW_season",
    "Ttotal~SOS",
    "Ttotal~Fixed_Trate",
    "Ttotal~P_pre",
    "Ttotal~T_pre",
    "Ttotal~SW_pre",
    "Ttotal~P_season",
    "Ttotal~T_season",
    "Ttotal~SW_season",
    "P_pre_via_SOS_FixedTrate",
    "T_pre_via_SOS_FixedTrate",
    "SW_pre_via_SOS_FixedTrate"
  )

  coef_sum <- setNames(rep(0, length(coef_names)), coef_names)
  coef_sumsq <- setNames(rep(0, length(coef_names)), coef_names)
  coef_count <- setNames(rep(0, length(coef_names)), coef_names)
  sig_sum <- setNames(rep(0, length(coef_names)), coef_names)
  sig_sumsq <- setNames(rep(0, length(coef_names)), coef_names)
  sig_count <- setNames(rep(0, length(coef_names)), coef_names)

  r2_names <- c("SOS", "Fixed_Trate", "Ttotal")
  r2_sum <- setNames(rep(0, length(r2_names)), r2_names)
  r2_sumsq <- setNames(rep(0, length(r2_names)), r2_names)
  r2_count <- setNames(rep(0, length(r2_names)), r2_names)
  zero_coef <- setNames(rep(0, length(coef_names)), coef_names)
  zero_r2 <- setNames(rep(0, length(r2_names)), r2_names)
  zero_filter <- c(coef_extreme = 0, p_invalid = 0, r2_invalid = 0)
  filter_info <- zero_filter

  process_block <- function(i) {
    row <- bs$row[i]
    nrows <- bs$nrows[i]

    mask_vals <- getValues(mask_r, row = row, nrows = nrows)
    valid_mask_cells <- which(!is.na(mask_vals))
    if (length(valid_mask_cells) == 0) {
      return(list(
        coef_sum = zero_coef,
        coef_sumsq = zero_coef,
        coef_count = zero_coef,
        sig_sum = zero_coef,
        sig_sumsq = zero_coef,
        sig_count = zero_coef,
        r2_sum = zero_r2,
        r2_sumsq = zero_r2,
        r2_count = zero_r2,
        filter_info = zero_filter,
        processed_cells = 0
      ))
    }

    n_cells <- length(mask_vals)
    n_years <- length(years)

    fixed_len_block <- getValues(fixed_window_length_r, row = row, nrows = nrows)
    fixed_len_block <- sanitize_values(fixed_len_block, NAvalue(fixed_window_length_r), allow_negative = FALSE)
    fixed_len_block[fixed_len_block <= 0] <- NA

    tr_block <- t(getValues(stacks$TR_fixed_window, row = row, nrows = nrows))
    tr_block <- sanitize_values(tr_block, na_values$TR_fixed_window, allow_negative = TRUE)
    fixed_trate_block <- tr_block / matrix(fixed_len_block, nrow = n_years, ncol = n_cells, byrow = TRUE)

    ttotal_block <- t(getValues(stacks$TRc_annual, row = row, nrows = nrows))
    ttotal_block <- sanitize_values(ttotal_block, na_values$TRc_annual, allow_negative = FALSE)

    sos_block <- t(getValues(stacks$SOS, row = row, nrows = nrows))
    sos_block <- sanitize_values(sos_block, na_values$SOS, allow_negative = FALSE)
    gpp_block <- t(getValues(stacks$GPP_season, row = row, nrows = nrows))
    gpp_block <- sanitize_values(gpp_block, na_values$GPP_season, allow_negative = FALSE)
    p_pre_block <- t(getValues(stacks$P_pre, row = row, nrows = nrows))
    p_pre_block <- sanitize_values(p_pre_block, na_values$P_pre, allow_negative = FALSE)
    t_pre_block <- t(getValues(stacks$T_pre, row = row, nrows = nrows))
    t_pre_block <- sanitize_values(t_pre_block, na_values$T_pre, allow_negative = TRUE)
    sw_pre_block <- t(getValues(stacks$SW_pre, row = row, nrows = nrows))
    sw_pre_block <- sanitize_values(sw_pre_block, na_values$SW_pre, allow_negative = FALSE)
    p_season_block <- t(getValues(stacks$P_season, row = row, nrows = nrows))
    p_season_block <- sanitize_values(p_season_block, na_values$P_season, allow_negative = FALSE)
    t_season_block <- t(getValues(stacks$T_season, row = row, nrows = nrows))
    t_season_block <- sanitize_values(t_season_block, na_values$T_season, allow_negative = TRUE)
    sw_season_block <- t(getValues(stacks$SW_season, row = row, nrows = nrows))
    sw_season_block <- sanitize_values(sw_season_block, na_values$SW_season, allow_negative = FALSE)

    valid_mat <- is.finite(ttotal_block) & is.finite(fixed_trate_block) & is.finite(sos_block) &
                 is.finite(p_pre_block) & is.finite(t_pre_block) & is.finite(sw_pre_block) &
                 is.finite(p_season_block) & is.finite(t_season_block) & is.finite(sw_season_block)
    valid_counts <- colSums(valid_mat)
    valid_cells <- intersect(valid_mask_cells, which(valid_counts >= min_years))
    if (length(valid_cells) == 0) {
      return(list(
        coef_sum = zero_coef,
        coef_sumsq = zero_coef,
        coef_count = zero_coef,
        sig_sum = zero_coef,
        sig_sumsq = zero_coef,
        sig_count = zero_coef,
        r2_sum = zero_r2,
        r2_sumsq = zero_r2,
        r2_count = zero_r2,
        filter_info = zero_filter,
        processed_cells = 0
      ))
    }

    coef_sum_block <- zero_coef
    coef_sumsq_block <- zero_coef
    coef_count_block <- zero_coef
    sig_sum_block <- zero_coef
    sig_sumsq_block <- zero_coef
    sig_count_block <- zero_coef
    r2_sum_block <- zero_r2
    r2_sumsq_block <- zero_r2
    r2_count_block <- zero_r2
    filter_info_block <- zero_filter

    for (idx in valid_cells) {
      valid_years <- valid_mat[, idx]
      if (sum(valid_years) < min_years) {
        next
      }

      df <- data.frame(
        year = years[valid_years],
        Ttotal = ttotal_block[valid_years, idx],
        Fixed_Trate = fixed_trate_block[valid_years, idx],
        SOS = sos_block[valid_years, idx],
        P_pre = p_pre_block[valid_years, idx],
        T_pre = t_pre_block[valid_years, idx],
        SW_pre = sw_pre_block[valid_years, idx],
        P_season = p_season_block[valid_years, idx],
        T_season = t_season_block[valid_years, idx],
        SW_season = sw_season_block[valid_years, idx]
      )

      if (DETREND_PIXEL_ENABLE) {
        df$Ttotal <- detrend_series(df$Ttotal, df$year)
        df$Fixed_Trate <- detrend_series(df$Fixed_Trate, df$year)
        df$SOS <- detrend_series(df$SOS, df$year)
        df$P_pre <- detrend_series(df$P_pre, df$year)
        df$T_pre <- detrend_series(df$T_pre, df$year)
        df$SW_pre <- detrend_series(df$SW_pre, df$year)
        df$P_season <- detrend_series(df$P_season, df$year)
        df$T_season <- detrend_series(df$T_season, df$year)
        df$SW_season <- detrend_series(df$SW_season, df$year)
      }

      ttotal_z <- scale_vec(df$Ttotal)
      fixed_trate_z <- scale_vec(df$Fixed_Trate)
      sos_z <- scale_vec(df$SOS)
      p_pre_z <- scale_vec(df$P_pre)
      t_pre_z <- scale_vec(df$T_pre)
      sw_pre_z <- scale_vec(df$SW_pre)
      p_season_z <- scale_vec(df$P_season)
      t_season_z <- scale_vec(df$T_season)
      sw_season_z <- scale_vec(df$SW_season)

      if (any(is.na(c(ttotal_z, fixed_trate_z, sos_z, p_pre_z, t_pre_z, sw_pre_z,
                      p_season_z, t_season_z, sw_season_z)))) {
        next
      }

      X_sos <- cbind(p_pre_z, t_pre_z, sw_pre_z)
      a_res <- regress_beta_p(X_sos, sos_z)
      r2_sos <- calc_r2(X_sos, sos_z)
      a <- a_res$beta
      p_a <- a_res$p

      X_trate <- cbind(sos_z, p_pre_z, t_pre_z, sw_pre_z, p_season_z, t_season_z, sw_season_z)
      bc_res <- regress_beta_p(X_trate, fixed_trate_z)
      r2_trate <- calc_r2(X_trate, fixed_trate_z)
      bc <- bc_res$beta
      p_bc <- bc_res$p
      b <- bc[1]
      f <- bc[2:4]
      c <- bc[5:7]
      p_b <- p_bc[1]
      p_f <- p_bc[2:4]
      p_c <- p_bc[5:7]

      vars_tr <- c("SOS", "Fixed_Trate", "P_pre", "T_pre", "SW_pre",
                   "P_season", "T_season", "SW_season")
      X_tr_full <- cbind(sos_z, fixed_trate_z, p_pre_z, t_pre_z, sw_pre_z,
                         p_season_z, t_season_z, sw_season_z)
      X_tr <- X_tr_full
      vars_current <- vars_tr

      if (ncol(X_tr) > 1) {
        repeat {
          vifs <- calc_vif_vector(X_tr)
          if (all(is.finite(vifs)) && max(vifs) > 10 && ncol(X_tr) > 1) {
            rm_idx <- which.max(vifs)
            X_tr <- X_tr[, -rm_idx, drop = FALSE]
            vars_current <- vars_current[-rm_idx]
            if (ncol(X_tr) <= 1) {
              break
            }
          } else {
            break
          }
        }
      }

      d <- rep(NA_real_, length(vars_tr))
      p_d <- rep(NA_real_, length(vars_tr))

      if (ncol(X_tr) >= 1) {
        d_res <- regress_beta_p(X_tr, ttotal_z)
        r2_tr <- calc_r2(X_tr, ttotal_z)
        if (!any(is.na(d_res$beta)) && !any(is.na(d_res$p))) {
          for (j in seq_along(vars_current)) {
            idx_match <- match(vars_current[j], vars_tr)
            d[idx_match] <- d_res$beta[j]
            p_d[idx_match] <- d_res$p[j]
          }
        }
      } else {
        r2_tr <- NA_real_
      }

      if (any(is.na(c(a, b, c, f))) || any(is.na(c(p_a, p_b, p_c, p_f))) || all(is.na(d))) {
        next
      }

      vals <- c(
        a[1], a[2], a[3],
        b,
        f[1], f[2], f[3],
        c[1], c[2], c[3],
        d[1], d[2], d[3], d[4], d[5], d[6], d[7], d[8],
        a[1] * b * d[2],
        a[2] * b * d[2],
        a[3] * b * d[2]
      )

      p_vals <- c(
        p_a[1], p_a[2], p_a[3],
        p_b,
        p_f[1], p_f[2], p_f[3],
        p_c[1], p_c[2], p_c[3],
        p_d[1], p_d[2], p_d[3], p_d[4], p_d[5], p_d[6], p_d[7], p_d[8],
        NA_real_, NA_real_, NA_real_
      )

      r2_vals <- c(r2_sos, r2_trate, r2_tr)
      filter_res <- filter_sem_outputs(vals, p_vals, r2_vals)
      vals <- filter_res$vals
      p_vals <- filter_res$p_vals
      r2_vals <- filter_res$r2_vals
      filter_info_block <- filter_info_block + filter_res$info

      for (k in seq_along(coef_names)) {
        if (is.finite(vals[k])) {
          coef_sum_block[coef_names[k]] <- coef_sum_block[coef_names[k]] + vals[k]
          coef_sumsq_block[coef_names[k]] <- coef_sumsq_block[coef_names[k]] + vals[k]^2
          coef_count_block[coef_names[k]] <- coef_count_block[coef_names[k]] + 1
          if (is.finite(p_vals[k]) && p_vals[k] < 0.05) {
            sig_sum_block[coef_names[k]] <- sig_sum_block[coef_names[k]] + vals[k]
            sig_sumsq_block[coef_names[k]] <- sig_sumsq_block[coef_names[k]] + vals[k]^2
            sig_count_block[coef_names[k]] <- sig_count_block[coef_names[k]] + 1
          }
        }
      }

      for (k in seq_along(r2_names)) {
        if (is.finite(r2_vals[k])) {
          r2_sum_block[r2_names[k]] <- r2_sum_block[r2_names[k]] + r2_vals[k]
          r2_sumsq_block[r2_names[k]] <- r2_sumsq_block[r2_names[k]] + r2_vals[k]^2
          r2_count_block[r2_names[k]] <- r2_count_block[r2_names[k]] + 1
        }
      }
    }

    list(
      coef_sum = coef_sum_block,
      coef_sumsq = coef_sumsq_block,
      coef_count = coef_count_block,
      sig_sum = sig_sum_block,
      sig_sumsq = sig_sumsq_block,
      sig_count = sig_count_block,
      r2_sum = r2_sum_block,
      r2_sumsq = r2_sumsq_block,
      r2_count = r2_count_block,
      filter_info = filter_info_block,
      processed_cells = length(valid_cells)
    )
  }

  block_indices <- seq_len(bs$n)
  blocks_done <- 0

  combine_block <- function(res) {
    coef_sum <<- coef_sum + res$coef_sum
    coef_sumsq <<- coef_sumsq + res$coef_sumsq
    coef_count <<- coef_count + res$coef_count
    sig_sum <<- sig_sum + res$sig_sum
    sig_sumsq <<- sig_sumsq + res$sig_sumsq
    sig_count <<- sig_count + res$sig_count
    r2_sum <<- r2_sum + res$r2_sum
    r2_sumsq <<- r2_sumsq + res$r2_sumsq
    r2_count <<- r2_count + res$r2_count
    filter_info <<- filter_info + res$filter_info
    processed_cells <<- processed_cells + res$processed_cells
    blocks_done <<- blocks_done + 1
  }

  if (PIXEL_BLOCK_PARALLEL && PIXEL_BLOCK_CORES > 1) {
    cat(sprintf("  启用并行块处理: %d cores\n", PIXEL_BLOCK_CORES))
    cl <- makeCluster(PIXEL_BLOCK_CORES)
    clusterEvalQ(cl, library(raster))
    clusterExport(
      cl,
      c("bs", "mask_r", "fixed_window_length_r", "stacks", "na_values", "years",
        "min_years", "DETREND_PIXEL_ENABLE", "coef_names", "r2_names",
        "sanitize_values", "calc_r2", "regress_beta_p", "calc_vif_vector",
        "detrend_series", "scale_vec", "zero_coef", "zero_r2", "process_block",
        "filter_sem_outputs", "FILTER_SEM_OUTLIERS", "SEM_COEF_ABS_MAX",
        "SEM_P_MIN", "SEM_P_MAX", "SEM_R2_MIN", "SEM_R2_MAX",
        "zero_filter", "NODATA_OUT", "NODATA_ABS_MAX"),
      envir = environment()
    )

    chunk_size <- max(1, PIXEL_BLOCK_CHUNK)
    chunk_list <- split(block_indices, ceiling(seq_along(block_indices) / chunk_size))
    for (chunk_idx in seq_along(chunk_list)) {
      res_list <- parLapply(cl, chunk_list[[chunk_idx]], process_block)
      for (res in res_list) {
        combine_block(res)
      }

      if (blocks_done %% 10 == 0 || blocks_done == length(block_indices)) {
        cat(sprintf("  进度: 块 %d/%d (%.1f%%), 已处理像元 %d/%d (%.1f%%)\n",
                    blocks_done, length(block_indices), 100 * blocks_done / length(block_indices),
                    processed_cells, total_valid_cells,
                    100 * processed_cells / max(1, total_valid_cells)))
      }
    }
    stopCluster(cl)
  } else {
    for (i in block_indices) {
      res <- process_block(i)
      combine_block(res)
      if (blocks_done %% 10 == 0 || blocks_done == length(block_indices)) {
        cat(sprintf("  进度: 块 %d/%d (%.1f%%), 已处理像元 %d/%d (%.1f%%)\n",
                    blocks_done, length(block_indices), 100 * blocks_done / length(block_indices),
                    processed_cells, total_valid_cells,
                    100 * processed_cells / max(1, total_valid_cells)))
      }
    }
  }

  mean_vals <- coef_sum / coef_count
  sd_vals <- sqrt(pmax(0, coef_sumsq / coef_count - mean_vals^2))
  sig_frac <- sig_count / coef_count
  mean_sig <- sig_sum / sig_count
  sd_sig <- sqrt(pmax(0, sig_sumsq / sig_count - mean_sig^2))

  mean_vals[!is.finite(mean_vals)] <- NA_real_
  sd_vals[!is.finite(sd_vals)] <- NA_real_
  sig_frac[!is.finite(sig_frac)] <- NA_real_
  mean_sig[!is.finite(mean_sig)] <- NA_real_
  sd_sig[!is.finite(sd_sig)] <- NA_real_

  if (FILTER_SEM_OUTLIERS) {
    cat(sprintf("\n  输出异常值过滤（像元级）: 系数极端=%d, p值异常=%d, R²异常=%d\n",
                filter_info["coef_extreme"], filter_info["p_invalid"], filter_info["r2_invalid"]))
  }

  pixel_filter_df <- data.frame(
    metric = c("coef_extreme", "p_invalid", "r2_invalid"),
    count = as.integer(filter_info)
  )
  write.csv(pixel_filter_df,
            file.path(PIXELWISE_DIR, "SEM_dual_timescale_pixelwise_outlier_filtering.csv"),
            row.names = FALSE)

  summary_df <- data.frame(
    path = coef_names,
    mean = mean_vals,
    sd = sd_vals,
    n = coef_count,
    sig_n = sig_count,
    sig_frac = sig_frac,
    mean_sig = mean_sig,
    sd_sig = sd_sig
  )

  write.csv(summary_df, file.path(PIXELWISE_DIR, "SEM_dual_timescale_pixel_time_series_summary.csv"),
            row.names = FALSE)
  write.csv(summary_df, file.path(PIXELWISE_DIR, "SEM_dual_timescale_parameters.csv"),
            row.names = FALSE)

  r2_mean <- r2_sum / r2_count
  r2_sd <- sqrt(pmax(0, r2_sumsq / r2_count - r2_mean^2))
  r2_mean[!is.finite(r2_mean)] <- NA_real_
  r2_sd[!is.finite(r2_sd)] <- NA_real_

  r2_summary <- data.frame(
    variable = r2_names,
    R2 = r2_mean,
    R2_sd = r2_sd,
    n = r2_count
  )

  write.csv(r2_summary, file.path(PIXELWISE_DIR, "SEM_dual_timescale_R2_detail.csv"),
            row.names = FALSE)
  write.csv(r2_summary[, c("variable", "R2")],
            file.path(PIXELWISE_DIR, "SEM_dual_timescale_R2.csv"),
            row.names = FALSE)

  sink(file.path(PIXELWISE_DIR, "SEM_dual_timescale_summary.txt"))
  cat("像元时间序列SEM汇总\n\n")
  print(summary_df)
  if (FILTER_SEM_OUTLIERS) {
    cat("\n输出异常值过滤统计（像元级）:\n")
    print(pixel_filter_df)
  }
  cat("\nR2汇总:\n")
  print(r2_summary)
  sink()

  cat("\n像元时间序列SEM摘要:\n")
  print(summary_df)

  summary_df
}

# ==================== SEM分析（全局合并）====================
# ===【修改7】添加POSav气候态计算===
run_dual_timescale_sem <- function(years, mask_r) {
  cat("\n=== 计算物候气候态（像元多年平均：SOSav, POSav） ===\n")

  # 读取所有年份的SOS和POS数据并计算多年平均
  sos_stack_list <- list()
  pos_stack_list <- list()
  n_loaded_sos <- 0
  n_loaded_pos <- 0

  for (year in years) {
    sos_file <- file.path(PHENO_DIR, sprintf("sos_gpp_%d.tif", year))
    pos_file <- file.path(PHENO_DIR, sprintf("pos_doy_gpp_%d.tif", year))

    if (file.exists(sos_file)) {
      sos_r <- raster(sos_file)
      sos_r <- set_nodata_if_missing(sos_r)
      sos_r <- mask_raster(sos_r, mask_r)
      sos_stack_list[[length(sos_stack_list) + 1]] <- sos_r
      n_loaded_sos <- n_loaded_sos + 1
    }

    if (file.exists(pos_file)) {
      pos_r <- raster(pos_file)
      pos_r <- set_nodata_if_missing(pos_r)
      pos_r <- mask_raster(pos_r, mask_r)
      pos_stack_list[[length(pos_stack_list) + 1]] <- pos_r
      n_loaded_pos <- n_loaded_pos + 1
    }

    if ((n_loaded_sos %% 10 == 0) && (n_loaded_sos > 0)) {
      cat(sprintf("  已读取 %d/%d 年份的物候数据...\n", n_loaded_sos, length(years)))
    }
  }

  cat(sprintf("  ✓ 成功读取 %d/%d 年份的SOS数据\n", n_loaded_sos, length(years)))
  cat(sprintf("  ✓ 成功读取 %d/%d 年份的POS数据\n", n_loaded_pos, length(years)))

  if (n_loaded_sos == 0 || n_loaded_pos == 0) {
    stop("未找到足够的SOS或POS数据文件")
  }

  # 将列表转换为stack并计算多年平均
  cat("  正在计算SOS多年平均（SOSav）...\n")
  sos_stack <- stack(sos_stack_list)
  sos_climatology_r <- calc(sos_stack, fun = function(x) {
    mean(x, na.rm = TRUE)
  })
  sos_climatology_r <- mask_raster(sos_climatology_r, mask_r)

  cat("  正在计算POS多年平均（POSav）...\n")
  pos_stack <- stack(pos_stack_list)
  pos_climatology_r <- calc(pos_stack, fun = function(x) {
    mean(x, na.rm = TRUE)
  })
  pos_climatology_r <- mask_raster(pos_climatology_r, mask_r)

  # ===【修改：先统计掩膜内有效像元，再验证POS > SOS】===
  mask_vals <- getValues(mask_r)
  in_mask <- !is.na(mask_vals) & (mask_vals > 0)

  sos_clim_vals_pre <- getValues(sos_climatology_r)
  sos_clim_vals_pre <- sanitize_values(sos_clim_vals_pre, NAvalue(sos_climatology_r), FALSE)
  pos_clim_vals_pre <- getValues(pos_climatology_r)
  pos_clim_vals_pre <- sanitize_values(pos_clim_vals_pre, NAvalue(pos_climatology_r), FALSE)

  valid_sos_pre <- in_mask & !is.na(sos_clim_vals_pre)
  valid_pos_pre <- in_mask & !is.na(pos_clim_vals_pre)
  valid_overlap_pre <- valid_sos_pre & valid_pos_pre

  cat(sprintf("  掩膜内有效像元(预过滤): SOSav=%d, POSav=%d, 重叠=%d\n",
              sum(valid_sos_pre), sum(valid_pos_pre), sum(valid_overlap_pre)))

  cat("  验证POSav > SOSav...\n")
  valid_window <- valid_overlap_pre & (pos_clim_vals_pre > sos_clim_vals_pre)
  n_invalid <- sum(valid_overlap_pre & !valid_window)
  cat(sprintf("  POSav > SOSav 过滤掉的像元数: %d\n", n_invalid))
  if (n_invalid > 0) {
    cat(sprintf("  ⚠️ 过滤了 %d 个POSav≤SOSav的无效像元\n", n_invalid))
  }

  valid_window_mask <- setValues(raster(mask_r), as.integer(valid_window))
  sos_climatology_r <- mask(sos_climatology_r, valid_window_mask, maskvalue = 0)
  pos_climatology_r <- mask(pos_climatology_r, valid_window_mask, maskvalue = 0)

  # 统计气候态信息（过滤后）
  sos_clim_vals <- getValues(sos_climatology_r)
  sos_clim_vals <- sanitize_values(sos_clim_vals, NAvalue(sos_climatology_r), FALSE)
  pos_clim_vals <- getValues(pos_climatology_r)
  pos_clim_vals <- sanitize_values(pos_clim_vals, NAvalue(pos_climatology_r), FALSE)

  n_valid_sos <- sum(!is.na(sos_clim_vals))
  n_valid_pos <- sum(!is.na(pos_clim_vals))
  sos_min <- min(sos_clim_vals, na.rm = TRUE)
  sos_max <- max(sos_clim_vals, na.rm = TRUE)
  sos_mean <- mean(sos_clim_vals, na.rm = TRUE)
  pos_min <- min(pos_clim_vals, na.rm = TRUE)
  pos_max <- max(pos_clim_vals, na.rm = TRUE)
  pos_mean <- mean(pos_clim_vals, na.rm = TRUE)

  cat(sprintf("  ✓ 物候气候态计算完成:\n"))
  cat(sprintf("    SOSav - 有效像元: %d, 范围: %.1f - %.1f DOY, 平均: %.1f DOY\n",
              n_valid_sos, sos_min, sos_max, sos_mean))
  cat(sprintf("    POSav - 有效像元: %d, 范围: %.1f - %.1f DOY, 平均: %.1f DOY\n",
              n_valid_pos, pos_min, pos_max, pos_mean))

  # 保存物候气候态
  sos_clim_file <- file.path(DERIVED_DIR, "SOS_climatology.tif")
  pos_clim_file <- file.path(DERIVED_DIR, "POS_climatology.tif")
  writeRaster(sos_climatology_r, sos_clim_file, overwrite = TRUE, datatype = "FLT4S")
  writeRaster(pos_climatology_r, pos_clim_file, overwrite = TRUE, datatype = "FLT4S")
  cat(sprintf("  ✓ 物候气候态已保存: %s, %s\n", basename(sos_clim_file), basename(pos_clim_file)))

  # 读取Fixed_Window_Length（用于计算Fixed_Trate）
  cat("\n=== 读取Fixed_Window_Length ===\n")
  fixed_len_file <- file.path(DECOMP_DIR, "Fixed_Window_Length.tif")
  if (!file.exists(fixed_len_file)) {
    stop(sprintf("Fixed_Window_Length文件不存在: %s", fixed_len_file))
  }
  fixed_window_length_r <- raster(fixed_len_file)
  fixed_window_length_r <- set_nodata_if_missing(fixed_window_length_r)
  fixed_window_length_r <- mask_raster(fixed_window_length_r, mask_r)
  cat(sprintf("  ✓ Fixed_Window_Length已加载: %s\n", basename(fixed_len_file)))

  cat("\n=== 准备SEM数据 ===\n")

  sem_data <- data.frame()

  # ===【修改8】调用时传入POSav和fixed_window_length_r参数，数据框列名改为Fixed_Trate===
  for (year in years) {
    rasters <- prepare_dual_timescale_data(year, sos_climatology_r, pos_climatology_r,
                                          fixed_window_length_r, mask_r)
    if (is.null(rasters)) next

    if (SEM_SAMPLE_MODE == "annual_mean") {
      year_row <- data.frame(
        year = year,
        Ttotal = cellStats(rasters$ttotal, mean, na.rm = TRUE),
        Fixed_Trate = cellStats(rasters$fixed_trate, mean, na.rm = TRUE),
        SOS = cellStats(rasters$sos, mean, na.rm = TRUE),
        P_pre = cellStats(rasters$p_pre, mean, na.rm = TRUE),
        T_pre = cellStats(rasters$t_pre, mean, na.rm = TRUE),
        SW_pre = cellStats(rasters$sw_pre, mean, na.rm = TRUE),
        P_season = cellStats(rasters$p_season, mean, na.rm = TRUE),
        T_season = cellStats(rasters$t_season, mean, na.rm = TRUE),
        SW_season = cellStats(rasters$sw_season, mean, na.rm = TRUE)
      )
      sem_data <- rbind(sem_data, year_row)
      cat(sprintf("  年份 %d: 全域均值已添加\n", year))
    } else {
      # 提取像元值
      vals <- data.frame(
        year = year,
        Ttotal = sanitize_values(getValues(rasters$ttotal), NAvalue(rasters$ttotal), TRUE),
        Fixed_Trate = sanitize_values(getValues(rasters$fixed_trate), NAvalue(rasters$fixed_trate), TRUE),
        SOS = sanitize_values(getValues(rasters$sos), NAvalue(rasters$sos), FALSE),
        P_pre = sanitize_values(getValues(rasters$p_pre), NAvalue(rasters$p_pre), FALSE),
        T_pre = sanitize_values(getValues(rasters$t_pre), NAvalue(rasters$t_pre), TRUE),
        SW_pre = sanitize_values(getValues(rasters$sw_pre), NAvalue(rasters$sw_pre), FALSE),
        P_season = sanitize_values(getValues(rasters$p_season), NAvalue(rasters$p_season), FALSE),
        T_season = sanitize_values(getValues(rasters$t_season), NAvalue(rasters$t_season), TRUE),
        SW_season = sanitize_values(getValues(rasters$sw_season), NAvalue(rasters$sw_season), FALSE)
      )

      vals <- vals[complete.cases(vals), ]
      sem_data <- rbind(sem_data, vals)

      cat(sprintf("  年份 %d: %d 完整案例\n", year, nrow(vals)))
    }
  }

  # 数据质量检查
  cat("\n数据质量检查:\n")
  cat(sprintf("  总样本数: %d\n", nrow(sem_data)))
  cat(sprintf("  完整案例数: %d\n", sum(complete.cases(sem_data))))

  for (col in c("Ttotal", "Fixed_Trate", "SOS", "P_pre", "T_pre", "SW_pre",
                "P_season", "T_season", "SW_season")) {
    valid_n <- sum(is.finite(sem_data[[col]]))
    mean_val <- mean(sem_data[[col]], na.rm = TRUE)
    cat(sprintf("  %s: %d个有效值, 均值=%.3f\n", col, valid_n, mean_val))
  }

  # 检查是否有足够的有效数据
  min_samples <- ifelse(SEM_SAMPLE_MODE == "annual_mean", 3, 100)
  if (sum(complete.cases(sem_data)) < min_samples) {
    cat(sprintf("\n✗ 错误：完整案例数不足（<%d），无法进行SEM分析\n", min_samples))
    cat("可能原因:\n")
    cat("  1. 日尺度数据文件路径不正确\n")
    cat("  2. 日尺度数据缺失\n")
    cat("  3. 掩膜过于严格导致无有效像元\n")
    cat("\n请检查以下路径:\n")
    cat(sprintf("  GPP日尺度: %s\n", GPP_DAILY_DIR))
    cat(sprintf("  降水日尺度: %s\n", PRECIP_DAILY_DIR))
    cat(sprintf("  气温日尺度: %s\n", TA_DAILY_DIR))
    cat(sprintf("  短波辐射日尺度: %s\n", SW_DAILY_DIR))
    stop("数据准备失败：完整案例数不足")
  }

  # 可选：去趋势（仅对annual_mean生效）
  sem_data_for_sem <- sem_data
  if (DETREND_ENABLE) {
    if (SEM_SAMPLE_MODE != "annual_mean") {
      cat("⚠️ 去趋势仅对annual_mean模式生效，已跳过\n")
    } else {
      if (DETREND_METHOD != "linear") {
        cat(sprintf("⚠️ 未知去趋势方法: %s，已改用linear\n", DETREND_METHOD))
      }
      sem_data_for_sem <- apply_detrend_sem_data(sem_data_for_sem, "year")
      write.csv(sem_data_for_sem, file.path(DATA_DIR, "sem_dual_timescale_detrended.csv"), row.names = FALSE)
      cat("✓ 去趋势完成: sem_dual_timescale_detrended.csv\n")
    }
  }

  # Z-score标准化
  sem_data_std <- sem_data_for_sem
  sem_vars <- setdiff(names(sem_data_std), "year")
  sem_data_std[sem_vars] <- scale(sem_data_std[sem_vars])

  # 保存原始和标准化数据
  write.csv(sem_data, file.path(DATA_DIR, "sem_dual_timescale_raw.csv"), row.names = FALSE)
  write.csv(sem_data_std, file.path(DATA_DIR, "sem_dual_timescale_standardized.csv"), row.names = FALSE)

  cat("✓ SEM数据准备完成\n")

  # ===== 构建双时间尺度SEM模型 =====
  # ===【修改9】SEM模型：Fixed_Trate替代TRc===
  cat("\n=== 构建双时间尺度SEM模型（固定窗口方法） ===\n")
  cat("✓ 气候因子: 降水(P)、气温(T)、短波辐射(SW)\n")
  cat("✓ 时间窗口: 季前气候（SOSav前3个月）+ 生长季气候（固定窗口[SOSav, POSav]）\n")
  cat("✓ 因变量: Ttotal（年蒸腾总量）\n")

  # 完整路径模型（降水 + 温度 + 短波辐射）
  sem_model <- '
    # 第一层：季前气候 → SOS（物候响应）
    SOS ~ a1*P_pre + a2*T_pre + a3*SW_pre

    # 第二层：SOS + 季前气候 + 生长季气候 → Fixed_Trate（蒸腾速率）
    Fixed_Trate ~ b*SOS + f1*P_pre + f2*T_pre + f3*SW_pre +
            c1*P_season + c2*T_season + c3*SW_season

    # 第三层：SOS + Fixed_Trate + 季前气候 + 生长季气候 → Ttotal（年蒸腾总量）
    Ttotal ~ g*SOS + d*Fixed_Trate +
             h1*P_pre + h2*T_pre + h3*SW_pre +
             e1*P_season + e2*T_season + e3*SW_season

    # === 间接效应分解 ===

    # 季前气候的间接效应路径
    P_pre_via_SOS_FixedTrate  := a1 * b * d
    T_pre_via_SOS_FixedTrate  := a2 * b * d
    SW_pre_via_SOS_FixedTrate := a3 * b * d

    P_pre_via_SOS  := a1 * g
    T_pre_via_SOS  := a2 * g
    SW_pre_via_SOS := a3 * g

    P_pre_via_FixedTrate  := f1 * d
    T_pre_via_FixedTrate  := f2 * d
    SW_pre_via_FixedTrate := f3 * d

    # 季前气候的总间接效应
    P_pre_indirect  := a1*b*d + a1*g + f1*d
    T_pre_indirect  := a2*b*d + a2*g + f2*d
    SW_pre_indirect := a3*b*d + a3*g + f3*d

    # 生长季气候通过Fixed_Trate的间接效应
    P_season_via_FixedTrate  := c1 * d
    T_season_via_FixedTrate  := c2 * d
    SW_season_via_FixedTrate := c3 * d

    # Fixed_Trate的中介比例
    P_FixedTrate_mediation := (c1*d) / (e1 + c1*d)
    T_FixedTrate_mediation := (c2*d) / (e2 + c2*d)
    SW_FixedTrate_mediation := (c3*d) / (e3 + c3*d)
  '

  cat("\nSEM模型结构:\n")
  cat(sem_model)
  cat("\n")

  # ===== 拟合SEM模型 =====
  cat("\n=== 拟合SEM模型 ===\n")

  fit <- sem(sem_model, data = sem_data_std, estimator = "MLR")

  cat("\nSEM拟合摘要:\n")
  fit_summary <- summary(fit, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
  print(fit_summary)

  # 保存结果
  params <- parameterEstimates(fit, standardized = TRUE)
  param_filter <- filter_sem_param_table(params)
  params <- param_filter$params
  if (FILTER_SEM_OUTLIERS) {
    cat(sprintf("\n  输出异常值过滤（全域SEM）：系数极端=%d, p值异常=%d\n",
                param_filter$info["coef_extreme"], param_filter$info["p_invalid"]))
  }
  write.csv(params, file.path(OUTPUT_DIR, "SEM_dual_timescale_parameters.csv"), row.names = FALSE)

  fit_measures <- fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))
  write.csv(t(as.data.frame(fit_measures)), file.path(OUTPUT_DIR, "SEM_dual_timescale_fitindices.csv"))

  r2_vals <- inspect(fit, "r2")
  r2_filter <- filter_sem_r2(r2_vals)
  r2_vals <- r2_filter$r2
  if (FILTER_SEM_OUTLIERS) {
    cat(sprintf("  输出异常值过滤（全域SEM）：R²异常=%d\n",
                r2_filter$info["r2_invalid"]))
  }
  write.csv(as.data.frame(r2_vals), file.path(OUTPUT_DIR, "SEM_dual_timescale_R2.csv"))

  global_filter_df <- data.frame(
    metric = c("coef_extreme", "p_invalid", "r2_invalid"),
    count = as.integer(c(param_filter$info["coef_extreme"],
                         param_filter$info["p_invalid"],
                         r2_filter$info["r2_invalid"]))
  )
  write.csv(global_filter_df,
            file.path(OUTPUT_DIR, "SEM_dual_timescale_outlier_filtering.csv"),
            row.names = FALSE)

  sink(file.path(OUTPUT_DIR, "SEM_dual_timescale_summary.txt"))
  print(fit_summary)
  if (FILTER_SEM_OUTLIERS) {
    cat("\n输出异常值过滤统计（全域SEM）:\n")
    print(global_filter_df)
  }
  sink()

  # 路径图
  pdf(file.path(OUTPUT_DIR, "SEM_dual_timescale_pathdiagram.pdf"), width = 12, height = 10)
  semPaths(
    fit,
    what = "std",
    layout = "tree2",
    edge.label.cex = 0.8,
    curvePivot = TRUE,
    fade = FALSE,
    residuals = FALSE,
    intercepts = FALSE,
    nCharNodes = 0,
    sizeMan = 8,
    sizeLat = 10,
    title = TRUE,
    mar = c(3, 3, 3, 3)
  )
  dev.off()

  cat("✓ 双时间尺度SEM分析完成\n")
  cat(sprintf("  输出目录: %s\n", OUTPUT_DIR))

  fit
}

# ==================== 主函数 ====================
main <- function() {
  cat("\n======================================================================\n")
  cat("双时间尺度SEM分析 - 与N04路径结构对齐\n")
  cat("======================================================================\n")
  cat(sprintf("分析模式: %s\n", SEM_SAMPLE_MODE))
  cat(sprintf("年份范围: %d-%d\n", YEAR_START, YEAR_END))
  cat("----------------------------------------------------------------------\n")

  years <- YEAR_START:YEAR_END

  # 诊断日尺度数据文件（使用第一年作为测试）
  cat("\n[日尺度文件诊断 - 测试年份: 1982]\n")
  diagnose_daily_files(1982, GPP_DAILY_DIR, GPP_DAILY_PATTERN, sample_dates = 3)
  diagnose_daily_files(1982, PRECIP_DAILY_DIR, PRECIP_DAILY_PATTERN, sample_dates = 3)
  diagnose_daily_files(1982, TA_DAILY_DIR, TA_DAILY_PATTERN, sample_dates = 3)
  diagnose_daily_files(1982, SW_DAILY_DIR, SW_DAILY_PATTERN, sample_dates = 3)

  # 读取并诊断掩膜
  cat("\n[掩膜诊断]\n")
  mask_r <- raster(MASK_FILE)
  cat(sprintf("  掩膜文件: %s\n", MASK_FILE))
  cat(sprintf("  栅格尺寸: %d x %d (行 x 列)\n", nrow(mask_r), ncol(mask_r)))
  cat(sprintf("  总像元数: %d\n", ncell(mask_r)))

  mask_vals_orig <- getValues(mask_r)
  n_valid_orig <- sum(mask_vals_orig > 0, na.rm = TRUE)
  cat(sprintf("  原始有效像元数 (>0): %d (%.2f%%)\n",
              n_valid_orig, 100 * n_valid_orig / ncell(mask_r)))

  mask_r[mask_r <= 0] <- NA
  mask_vals <- getValues(mask_r)
  n_valid <- sum(!is.na(mask_vals))
  cat(sprintf("  处理后有效像元数: %d (%.2f%%)\n",
              n_valid, 100 * n_valid / ncell(mask_r)))

  if (n_valid == 0) {
    cat("\n⚠️ 警告：掩膜处理后没有任何有效像元！\n")
    cat("  这会导致所有计算结果为NA\n")
    cat("  建议检查掩膜文件是否正确\n")
  }

  # 网格一致性硬检查（模板/掩膜/样本栅格）
  cat("\n[网格一致性检查]\n")
  template_path <- if (file.exists(TEMPLATE_FILE)) TEMPLATE_FILE else MASK_FILE
  template_r <- raster(template_path)
  cat(sprintf("  模板文件: %s\n", template_path))
  check_raster_alignment(template_r, mask_r, YEAR_START)

  # 测试读取多个时段的日尺度文件
  cat("\n[测试读取日尺度文件]\n")
  test_dates <- c("19820301", "19820701", "19821001")
  for (date_str in test_dates) {
    test_gpp_file <- file.path(GPP_DAILY_DIR, sprintf("GPP_%s.tif", date_str))
    cat(sprintf("  测试文件: %s\n", basename(test_gpp_file)))
    cat(sprintf("  文件存在: %s\n", ifelse(file.exists(test_gpp_file), "是", "否")))

    if (file.exists(test_gpp_file)) {
      tryCatch({
        test_r <- raster(test_gpp_file)
        cat(sprintf("  栅格尺寸: %d x %d\n", nrow(test_r), ncol(test_r)))
        test_vals <- getValues(test_r)
        n_valid_test <- sum(is.finite(test_vals) & test_vals > 0, na.rm = TRUE)
        cat(sprintf("  有效像元数 (>0): %d (%.2f%%)\n",
                    n_valid_test, 100 * n_valid_test / length(test_vals)))

        # 应用掩膜后
        test_r_masked <- mask(test_r, mask_r)
        test_vals_masked <- getValues(test_r_masked)
        n_valid_masked <- sum(is.finite(test_vals_masked) & test_vals_masked > 0, na.rm = TRUE)
        cat(sprintf("  应用掩膜后有效像元数: %d (%.2f%%)\n",
                    n_valid_masked, 100 * n_valid_masked / length(test_vals_masked)))
        cat(sprintf("  掩膜后均值: %.4f\n", mean(test_vals_masked, na.rm = TRUE)))
      }, error = function(e) {
        cat(sprintf("  ✗ 读取失败: %s\n", conditionMessage(e)))
      })
    }
  }

  cat("----------------------------------------------------------------------\n")

  # 运行双时间尺度SEM（全域均值 + 像元时间序列）
  sem_mode_orig <- SEM_SAMPLE_MODE
  detrend_orig <- DETREND_ENABLE
  SEM_SAMPLE_MODE <- "annual_mean"
  run_annual_mean_once <- function(enable_detrend, label) {
    DETREND_ENABLE <<- enable_detrend
    suffix <- if (enable_detrend) "_detrended" else ""
    set_output_dirs(suffix)
    cat(sprintf("\n=== 全域均值SEM（annual_mean）[%s] ===\n", label))
    fit <- run_dual_timescale_sem(years, mask_r)

    # VIF诊断
    sem_data_file <- file.path(DATA_DIR, "sem_dual_timescale_standardized.csv")
    if (file.exists(sem_data_file)) {
      data <- read.csv(sem_data_file)
      vif_vals <- calculate_vif(data)
      if (any(vif_vals > 10, na.rm = TRUE)) {
        cat("\n⚠️ 警告：检测到VIF > 10，建议移除高VIF变量\n")
      }
    }
    fit
  }

  if (RUN_BOTH_DETREND) {
    run_annual_mean_once(FALSE, "raw")
    run_annual_mean_once(TRUE, "detrended")
  } else {
    label <- if (DETREND_ENABLE) "detrended" else "raw"
    run_annual_mean_once(DETREND_ENABLE, label)
  }

  SEM_SAMPLE_MODE <- sem_mode_orig
  set_output_dirs("")
  DETREND_ENABLE <- detrend_orig

  # ===【修改10】像元时间序列模式：添加POSav气候态计算===
    # 先计算物候气候态（SOSav, POSav）
    cat("\n=== 计算物候气候态（像元多年平均：SOSav, POSav） ===\n")
    sos_clim_file <- file.path(DERIVED_DIR, "SOS_climatology.tif")
    pos_clim_file <- file.path(DERIVED_DIR, "POS_climatology.tif")

    if (USE_EXISTING_CLIMATOLOGY && file.exists(sos_clim_file) && file.exists(pos_clim_file)) {
      cat("  ✓ 复用已存在的物候气候态缓存\n")
      sos_climatology_r <- raster(sos_clim_file)
      pos_climatology_r <- raster(pos_clim_file)
      sos_climatology_r <- set_nodata_if_missing(sos_climatology_r)
      pos_climatology_r <- set_nodata_if_missing(pos_climatology_r)
      sos_climatology_r <- mask_raster(sos_climatology_r, mask_r)
      pos_climatology_r <- mask_raster(pos_climatology_r, mask_r)
    } else {
      sos_stack_list <- list()
      pos_stack_list <- list()
      n_loaded_sos <- 0
      n_loaded_pos <- 0

      for (year in years) {
        sos_file <- file.path(PHENO_DIR, sprintf("sos_gpp_%d.tif", year))
        pos_file <- file.path(PHENO_DIR, sprintf("pos_doy_gpp_%d.tif", year))

        if (file.exists(sos_file)) {
          sos_r <- raster(sos_file)
          sos_r <- set_nodata_if_missing(sos_r)
          sos_r <- mask_raster(sos_r, mask_r)
          sos_stack_list[[length(sos_stack_list) + 1]] <- sos_r
          n_loaded_sos <- n_loaded_sos + 1
        }

        if (file.exists(pos_file)) {
          pos_r <- raster(pos_file)
          pos_r <- set_nodata_if_missing(pos_r)
          pos_r <- mask_raster(pos_r, mask_r)
          pos_stack_list[[length(pos_stack_list) + 1]] <- pos_r
          n_loaded_pos <- n_loaded_pos + 1
        }

        if ((n_loaded_sos %% 10 == 0) && (n_loaded_sos > 0)) {
          cat(sprintf("  已读取 %d/%d 年份的物候数据...\n", n_loaded_sos, length(years)))
        }
      }

      cat(sprintf("  ✓ 成功读取 %d/%d 年份的SOS数据\n", n_loaded_sos, length(years)))
      cat(sprintf("  ✓ 成功读取 %d/%d 年份的POS数据\n", n_loaded_pos, length(years)))

      if (n_loaded_sos == 0 || n_loaded_pos == 0) {
        stop("未找到足够的SOS或POS数据文件")
      }

      cat("  正在计算SOS多年平均（SOSav）...\n")
      sos_stack <- stack(sos_stack_list)
      sos_climatology_r <- calc(sos_stack, fun = function(x) { mean(x, na.rm = TRUE) })
      sos_climatology_r <- mask_raster(sos_climatology_r, mask_r)

      cat("  正在计算POS多年平均（POSav）...\n")
      pos_stack <- stack(pos_stack_list)
      pos_climatology_r <- calc(pos_stack, fun = function(x) { mean(x, na.rm = TRUE) })
      pos_climatology_r <- mask_raster(pos_climatology_r, mask_r)
    }

    # ===【修改：添加POS > SOS验证】===
    cat("  验证POSav > SOSav...\n")
    valid_window_mask <- pos_climatology_r > sos_climatology_r
    valid_window_mask[is.na(valid_window_mask)] <- FALSE
    sos_climatology_r <- mask(sos_climatology_r, valid_window_mask, maskvalue = 0)
    pos_climatology_r <- mask(pos_climatology_r, valid_window_mask, maskvalue = 0)

    n_invalid <- sum(getValues(valid_window_mask) == 0, na.rm = TRUE)
    cat(sprintf("  过滤掉的像元数: %d\n", n_invalid))
    if (n_invalid > 0) {
      cat(sprintf("  ⚠️ 过滤了 %d 个POSav≤SOSav的无效像元\n", n_invalid))
    }

    # 统计气候态信息
    sos_clim_vals <- getValues(sos_climatology_r)
    sos_clim_vals <- sanitize_values(sos_clim_vals, NAvalue(sos_climatology_r), FALSE)
    pos_clim_vals <- getValues(pos_climatology_r)
    pos_clim_vals <- sanitize_values(pos_clim_vals, NAvalue(pos_climatology_r), FALSE)

    n_valid_sos <- sum(!is.na(sos_clim_vals))
    n_valid_pos <- sum(!is.na(pos_clim_vals))
    sos_min <- min(sos_clim_vals, na.rm = TRUE)
    sos_max <- max(sos_clim_vals, na.rm = TRUE)
    sos_mean <- mean(sos_clim_vals, na.rm = TRUE)
    pos_min <- min(pos_clim_vals, na.rm = TRUE)
    pos_max <- max(pos_clim_vals, na.rm = TRUE)
    pos_mean <- mean(pos_clim_vals, na.rm = TRUE)

    cat(sprintf("  ✓ 物候气候态计算完成:\n"))
    cat(sprintf("    SOSav - 有效像元: %d, 范围: %.1f - %.1f DOY, 平均: %.1f DOY\n",
                n_valid_sos, sos_min, sos_max, sos_mean))
    cat(sprintf("    POSav - 有效像元: %d, 范围: %.1f - %.1f DOY, 平均: %.1f DOY\n",
                n_valid_pos, pos_min, pos_max, pos_mean))

    # 保存物候气候态（仅在未复用缓存时写入）
    if (!(USE_EXISTING_CLIMATOLOGY && file.exists(sos_clim_file) && file.exists(pos_clim_file))) {
      writeRaster(sos_climatology_r, sos_clim_file, overwrite = TRUE, datatype = "FLT4S")
      writeRaster(pos_climatology_r, pos_clim_file, overwrite = TRUE, datatype = "FLT4S")
      cat(sprintf("  ✓ 物候气候态已保存: %s, %s\n", basename(sos_clim_file), basename(pos_clim_file)))
    }

    # 读取Fixed_Window_Length（用于计算Fixed_Trate）
    cat("\n=== 读取Fixed_Window_Length ===\n")
    fixed_len_file <- file.path(DECOMP_DIR, "Fixed_Window_Length.tif")
    if (!file.exists(fixed_len_file)) {
      stop(sprintf("Fixed_Window_Length文件不存在: %s", fixed_len_file))
    }
    fixed_window_length_r <- raster(fixed_len_file)
    fixed_window_length_r <- set_nodata_if_missing(fixed_window_length_r)
    fixed_window_length_r <- mask_raster(fixed_window_length_r, mask_r)
    cat(sprintf("  ✓ Fixed_Window_Length已加载: %s\n", basename(fixed_len_file)))

    # 准备缓存
    prepare_sem_caches(years, sos_climatology_r, pos_climatology_r,
                      fixed_window_length_r, mask_r)

    # 运行像元时间序列SEM
    run_sem_pixel_time_series(years, sos_climatology_r, fixed_window_length_r, mask_r)

  cat("\n======================================================================\n")
  cat("✓ 双时间尺度SEM分析完成\n")
  cat(sprintf("输出目录: %s\n", OUTPUT_DIR))
  cat("======================================================================\n")
}

# ==================== 运行 ====================
# 自动运行逻辑
# - 非交互模式（Rscript）：自动运行
# - 交互模式（RStudio）：提示用户手动调用
if (!interactive()) {
  main()
} else {
  cat("\n======================================================================\n")
  cat("⚠️ 交互模式检测\n")
  cat("======================================================================\n")
  cat("请手动执行以下命令以运行双时间尺度SEM分析:\n\n")
  cat("  main()\n\n")
  cat("或者使用以下命令在非交互模式运行:\n\n")
  cat("  Rscript 05b_SEM_analysis_dual_timescale_SOS_Ttotal_Trate.R\n")
  cat("======================================================================\n\n")
}
