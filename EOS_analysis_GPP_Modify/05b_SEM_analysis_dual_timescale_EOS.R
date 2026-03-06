#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Module 05b: 双时间尺度SEM分析（EOS） - 基于N04路径结构 + 04c固定窗口方法
#
# 设计思路：
# 1. 借鉴N04_Cal_SEM的双时间尺度路径结构
# 2. 使用04c代码的固定窗口方法（Fixed_Trate + 固定窗口[POSav, EOSav]）
# 3. 计算早季气候因子（固定窗口[SOSav, POSav]）
# 4. 计算晚季气候因子/GPP（固定窗口[POSav, EOSav]）
# 5. 构建完整路径：早/晚季气候 → EOS → Fixed_GPPrate → Fixed_Trate
# 6. 气候因子：降水(P)、气温(T)、短波辐射(SW)
#
# 路径结构（参考N04）：
#   早季气候(P,T,SW) + 晚季气候(P,T,SW) → EOS → Fixed_GPPrate → Fixed_Trate
#      ↓                           ↓         ↑          ↑
#      └───────────────────────────┴─────────┴──────────┘ (直接路径)
#   晚季气候(P,T,SW) ────────────────────────┴──────────┘ (同期路径)
#
# ⚠️ 关键修改（对比原05b版本）：
# 1. TRproduct → Fixed_Trate（数据源：TR_fixed_window_{year}.tif，临时计算）
# 2. Fixed_Trate = TR_fixed_window / Fixed_Window_Length（与04c一致）
# 3. 晚季气候：当年POS-EOS窗口 → 多年平均固定窗口[POSav, EOSav]
# 4. 晚季GPP：当年POS-EOS窗口 → 多年平均固定窗口[POSav, EOSav]
# 5. 早季气候：固定窗口[SOSav, POSav]
#
# Version: 2.0.1 (Bug Fix 5)
# Author: Wang2025 Replication Project
#
# Changelog:
#   v2.0.1 - 修复Bug Fix 5 (2025-01-12)
#            annual_mean模式未正确设置SEM_SAMPLE_MODE变量
#            导致annual_mean使用了pixel_time_series的数据(956,556行 vs 37行)
#            影响: annual_mean结果完全错误，需重新运行

suppressPackageStartupMessages({
  library(raster)
  library(lavaan)
  library(semPlot)
  library(parallel)  # 添加并行计算支持
})

# ===【统一并行配置】===
PARALLEL_ENABLE <- TRUE       # 并行总开关
PARALLEL_CORES <- 10          # CPU核心数（统一参数名）
PARALLEL_CHUNK_SIZE <- 5      # 块大小（用于块并行）
AUTO_DETECT_CORES <- TRUE     # 自动检测可用核心数

if (AUTO_DETECT_CORES && PARALLEL_ENABLE) {
  available_cores <- parallel::detectCores()
  if (is.na(available_cores) || available_cores < 1) {
    available_cores <- 1
  }
  PARALLEL_CORES <- max(1, min(PARALLEL_CORES, available_cores - 1))
  cat(sprintf("✓ 检测到 %d 核心，使用 %d 核心进行并行计算\n", available_cores, PARALLEL_CORES))
}

N_CORES <- PARALLEL_CORES  # 兼容旧代码
cat(sprintf("\n[并行化配置] 使用 %d 个CPU核心进行并行计算\n", PARALLEL_CORES))

# ===【运行模式】===
RUN_MODE <- tolower(Sys.getenv("WANG_RUN_MODE", "skip"))
OVERWRITE <- identical(RUN_MODE, "overwrite")
should_write <- function(path) {
  OVERWRITE || !file.exists(path)
}
safe_write_csv <- function(df, path, ...) {
  if (!should_write(path)) {
    cat(sprintf("  [skip] %s\n", path))
    return(invisible(FALSE))
  }
  write.csv(df, path, ...)
  invisible(TRUE)
}
safe_write_lines <- function(lines, path) {
  if (!should_write(path)) {
    cat(sprintf("  [skip] %s\n", path))
    return(invisible(FALSE))
  }
  writeLines(lines, path)
  invisible(TRUE)
}
safe_write_raster <- function(r, path, ...) {
  if (!should_write(path)) {
    cat(sprintf("  [skip] %s\n", path))
    return(invisible(FALSE))
  }
  writeRaster(r, path, overwrite = OVERWRITE, ...)
  invisible(TRUE)
}

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
OUTPUT_ROOT <- file.path(ROOT, "Wang2025_Analysis_EOS_GPP_Modify")
PHENO_DIR <- file.path(ROOT, "Phenology_Output_1", "GPP_phenology")
# ===【修改1】Fixed_Trate数据源（04c固定窗口方法）===
# 原路径: Decomposition/TRproduct_{year}.tif（03a模块输出）
# 新路径: Decomposition_FixedWindow/TR_fixed_window_{year}.tif（03c模块输出）
# 说明：Fixed_Trate需临时计算 = TR_fixed_window / Fixed_Window_Length
DECOMP_DIR <- file.path(OUTPUT_ROOT, "Decomposition_FixedWindow")
MASK_FILE <- file.path(OUTPUT_ROOT, "masks", "combined_mask.tif")
TEMPLATE_FILE <- file.path(OUTPUT_ROOT, "masks", "template_grid.tif")

# 日尺度数据路径
GPP_DAILY_DIR <- file.path(ROOT, "GLASS_GPP", "GLASS_GPP_daily_interpolated")
PRECIP_DAILY_DIR <- file.path(ROOT, "Meteorological Data", "ERA5_Land", "Pre", "Pre_Daily", "Pre_Daily_2")
TA_DAILY_DIR <- file.path(ROOT, "Meteorological Data", "ERA5_Land", "Tem", "Tem_Daily", "Tem_Daily_2")
SW_DAILY_DIR <- file.path(ROOT, "Meteorological Data", "ERA5_Land", "DSW", "DSW_Daily", "DSW_Daily_2")

# 输出目录（与05a命名风格保持一致：SEM_Data_*/SEM_Results_*）
DATA_DIR <- file.path(OUTPUT_ROOT, "SEM_Data_Dual_Fixed")
DERIVED_DIR <- file.path(DATA_DIR, "Derived")
OUTPUT_DIR <- file.path(OUTPUT_ROOT, "SEM_Results_Dual_Fixed")
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

annual_outputs_ready <- function(suffix = "") {
  out_dir <- if (!is.null(suffix) && nzchar(suffix)) {
    paste0(OUTPUT_DIR_BASE, suffix)
  } else {
    OUTPUT_DIR_BASE
  }
  data_dir <- if (!is.null(suffix) && nzchar(suffix)) {
    paste0(DATA_DIR_BASE, suffix)
  } else {
    DATA_DIR_BASE
  }
  required <- c(
    file.path(out_dir, "SEM_dual_timescale_parameters.csv"),
    file.path(out_dir, "SEM_dual_timescale_fitindices.csv"),
    file.path(out_dir, "SEM_dual_timescale_R2.csv"),
    file.path(out_dir, "SEM_dual_timescale_summary.txt"),
    file.path(data_dir, "sem_dual_timescale_raw.csv"),
    file.path(data_dir, "sem_dual_timescale_standardized.csv")
  )
  all(file.exists(required))
}

pixel_outputs_ready <- function(suffix = "") {
  pix_dir <- if (!is.null(suffix) && nzchar(suffix)) {
    file.path(paste0(OUTPUT_DIR_BASE, suffix), "Pixelwise")
  } else {
    file.path(OUTPUT_DIR_BASE, "Pixelwise")
  }
  required <- c(
    file.path(pix_dir, "SEM_dual_timescale_parameters.csv"),
    file.path(pix_dir, "SEM_dual_timescale_pixelwise_bootstrap_ci.csv"),
    file.path(pix_dir, "SEM_dual_timescale_R2.csv"),
    file.path(pix_dir, "SEM_dual_timescale_summary.txt")
  )
  all(file.exists(required))
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
MIN_VALID_FRAC <- 0.60            # POS-EOS窗口内有效数据最低比例
MIN_VALID_YEAR_FRAC <- 0.60       # 像元时间序列最少有效年份比例
WRITE_PIXELWISE_RASTERS <- FALSE  # 是否输出像元级栅格（会生成大量文件）
DETREND_ENABLE <- FALSE           # 是否启用去趋势（仅对annual_mean生效）
DETREND_METHOD <- "linear"        # 目前仅支持linear
DETREND_PIXEL_ENABLE <- TRUE      # 是否启用像元级去趋势（仅pixel_time_series生效）
# 注意：两个去趋势开关作用于不同模式，设置不匹配会被跳过并打印提示
RUN_BOTH_DETREND <- TRUE          # annual_mean 同时输出不去趋势+去趋势
FILTER_SEM_OUTLIERS <- TRUE       # 输出结果异常值过滤（类似04c）
SEM_COEF_ABS_MAX <- 5             # 系数绝对值阈值（标准化系数）
SEM_P_MIN <- 0                    # p值下限
SEM_P_MAX <- 1                    # p值上限
SEM_R2_MIN <- 0                   # R²下限
SEM_R2_MAX <- 1                   # R²上限
PIXELWISE_BOOTSTRAP_ENABLE <- TRUE   # 像元系数均值 bootstrap CI
PIXELWISE_BOOTSTRAP_N <- 800         # bootstrap次数（像元均值）
PIXELWISE_BOOTSTRAP_SEED <- 202501   # bootstrap随机种子

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
    return(list(beta = rep(NA_real_, ncol(X)),
                p = rep(NA_real_, ncol(X)),
                se = rep(NA_real_, ncol(X)),
                df = NA_integer_))
  }
  beta <- solve(XtX, crossprod(X, y))
  resid <- y - X %*% beta
  df <- nrow(X) - ncol(X)
  if (df <= 0) {
    return(list(beta = rep(NA_real_, ncol(X)),
                p = rep(NA_real_, ncol(X)),
                se = rep(NA_real_, ncol(X)),
                df = df))
  }
  sigma2 <- sum(resid^2) / df
  se <- sqrt(diag(sigma2 * solve(XtX)))
  tval <- beta / se
  pval <- 2 * (1 - pt(abs(tval), df))
  list(beta = as.vector(beta), p = as.vector(pval), se = as.vector(se), df = df)
}

delta_p_two <- function(a, b, se_a, se_b) {
  if (!is.finite(a) || !is.finite(b) || !is.finite(se_a) || !is.finite(se_b)) {
    return(NA_real_)
  }
  var_ab <- (b^2) * (se_a^2) + (a^2) * (se_b^2)
  if (!is.finite(var_ab) || var_ab <= 0) {
    return(NA_real_)
  }
  z <- (a * b) / sqrt(var_ab)
  2 * (1 - stats::pnorm(abs(z)))
}

delta_p_three <- function(a, b, c, se_a, se_b, se_c) {
  if (!is.finite(a) || !is.finite(b) || !is.finite(c) ||
      !is.finite(se_a) || !is.finite(se_b) || !is.finite(se_c)) {
    return(NA_real_)
  }
  var_abc <- (b^2 * c^2) * (se_a^2) +
    (a^2 * c^2) * (se_b^2) +
    (a^2 * b^2) * (se_c^2)
  if (!is.finite(var_abc) || var_abc <= 0) {
    return(NA_real_)
  }
  z <- (a * b * c) / sqrt(var_abc)
  2 * (1 - stats::pnorm(abs(z)))
}

delta_p_ratio <- function(c_val, d_val, e_val, se_c, se_d, se_e, eps = MEDIATION_DENOM_EPS) {
  if (!is.finite(c_val) || !is.finite(d_val) || !is.finite(e_val) ||
      !is.finite(se_c) || !is.finite(se_d) || !is.finite(se_e)) {
    return(NA_real_)
  }
  denom <- e_val + c_val * d_val
  if (!is.finite(denom) || abs(denom) < eps) {
    return(NA_real_)
  }
  d_c <- (d_val * e_val) / (denom^2)
  d_d <- (c_val * e_val) / (denom^2)
  d_e <- -(c_val * d_val) / (denom^2)
  var_r <- (d_c^2) * (se_c^2) + (d_d^2) * (se_d^2) + (d_e^2) * (se_e^2)
  if (!is.finite(var_r) || var_r <= 0) {
    return(NA_real_)
  }
  z <- (c_val * d_val / denom) / sqrt(var_r)
  2 * (1 - stats::pnorm(abs(z)))
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

bootstrap_ci_mean <- function(x, n_boot, seed = NULL) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2 || n_boot < 1) {
    return(c(mean = mean(x, na.rm = TRUE), ci_low = NA_real_, ci_high = NA_real_))
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  boot_means <- replicate(n_boot, mean(x[sample.int(n, n, replace = TRUE)]))
  c(
    mean = mean(x),
    ci_low = as.numeric(quantile(boot_means, 0.025, na.rm = TRUE)),
    ci_high = as.numeric(quantile(boot_means, 0.975, na.rm = TRUE))
  )
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
    list(file.path(PHENO_DIR, sprintf("eos_gpp_%d.tif", year)), "EOS")
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

# ===【EOS】固定窗口气候因子（通用窗口[start, end]均值）===
calc_window_climate_fixed <- function(year, start_climatology_r, end_climatology_r, daily_dir, pattern, mask_r,
                                      cache_file, var_name, allow_negative = TRUE, window_label = "") {
  if (file.exists(cache_file)) {
    cat(sprintf("    [缓存] %s\n", var_name))
    return(raster(cache_file))
  }

  if (nzchar(window_label)) {
    cat(sprintf("    计算 %s (固定窗口[%s]均值)...\n", var_name, window_label))
  } else {
    cat(sprintf("    计算 %s (固定窗口均值)...\n", var_name))
  }

  template <- mask_r
  # 使用多年平均物候（气候态）
  start_vals <- getValues(start_climatology_r)
  end_vals <- getValues(end_climatology_r)

  start_vals <- sanitize_values(start_vals, NAvalue(start_climatology_r), allow_negative = FALSE)
  end_vals <- sanitize_values(end_vals, NAvalue(end_climatology_r), allow_negative = FALSE)

  start_int <- as.integer(round(start_vals))
  end_int <- as.integer(round(end_vals))

  # 限制到1-365
  start_int <- pmax(1, pmin(365, start_int))
  end_int <- pmax(1, pmin(365, end_int))

  valid_pheno <- !is.na(start_int) & !is.na(end_int) & (start_int < end_int)

  # 初始化
  sum_vals <- rep(0, ncell(template))
  cnt_vals <- rep(0L, ncell(template))

  if (!any(valid_pheno)) {
    out_r <- setValues(template, rep(NA_real_, ncell(template)))
    out_r <- mask_raster(out_r, mask_r)
    safe_write_raster(out_r, cache_file)
    return(out_r)
  }

  # 确定DOY范围
  min_doy <- min(start_int[valid_pheno], na.rm = TRUE)
  max_doy <- max(end_int[valid_pheno], na.rm = TRUE)

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
      in_window <- valid_pheno & (start_int <= doy) & (end_int >= doy)
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

  window_len <- end_int - start_int + 1
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

  safe_write_raster(out_r, cache_file, datatype = "FLT4S")

  out_r
}

# ==================== 数据准备 ====================
# ===【修改3】数据准备函数：Fixed_Trate + 固定窗口===
prepare_dual_timescale_data <- function(year, sos_climatology_r, pos_climatology_r,
                                        eos_climatology_r, fixed_window_length_r, mask_r,
                                        parallel_inner = PARALLEL_ENABLE) {
  cat(sprintf("\n年份: %d\n", year))

  # ===【修改3a】读取TR_fixed_window并计算Fixed_Trate（替代TRc）===
  tr_fixed_file <- file.path(DECOMP_DIR, sprintf("TR_fixed_window_%d.tif", year))
  if (!file.exists(tr_fixed_file)) {
    cat(sprintf("  跳过: TR_fixed_window文件不存在 (%s)\n", basename(tr_fixed_file)))
    return(NULL)
  }

  # 读取当年物候数据（用于SEM中的EOS变量）
  eos_file <- file.path(PHENO_DIR, sprintf("eos_gpp_%d.tif", year))

  if (!file.exists(eos_file)) {
    cat(sprintf("  跳过: EOS文件不存在\n"))
    return(NULL)
  }

  tr_fixed_r <- raster(tr_fixed_file)
  eos_r <- raster(eos_file)

  tr_fixed_r <- set_nodata_if_missing(tr_fixed_r)
  eos_r <- set_nodata_if_missing(eos_r)

  # 应用掩膜
  eos_r <- mask_raster(eos_r, mask_r)
  tr_fixed_r <- mask_raster(tr_fixed_r, mask_r)

  # 计算Fixed_Trate = TR_fixed_window / Fixed_Window_Length
  cat("  [计算Fixed_Trate] TR_fixed_window / Fixed_Window_Length...\n")
  fixed_trate_r <- tr_fixed_r / fixed_window_length_r
  fixed_trate_r <- mask_raster(fixed_trate_r, mask_r)

  # 缓存文件路径
  cache_prefix <- file.path(DERIVED_DIR, sprintf("%s_%d.tif", c(
    "P_pre", "T_pre", "SW_pre",
    "P_season", "T_season", "SW_season",
    "Fixed_GPPrate"
  ), year))

  names(cache_prefix) <- c("P_pre", "T_pre", "SW_pre",
                           "P_season", "T_season", "SW_season",
                           "Fixed_GPPrate")

  # ===== 计算早季气候因子（固定窗口[SOSav, POSav]） =====
  cat("  [1/7] 早季气候因子（固定窗口[SOSav, POSav]）:\n")

  # ===【并行化优化】早季3个变量并行计算===
  if (parallel_inner && PARALLEL_CORES > 1) {
    cat(sprintf("    使用 %d 个核心并行计算...\n", PARALLEL_CORES))

    # 创建并行集群（Windows兼容）
    cl <- makeCluster(PARALLEL_CORES)
    on.exit(stopCluster(cl), add = TRUE)

    # 导出必要的变量和函数到集群
    clusterExport(cl, c("year", "sos_climatology_r", "pos_climatology_r", "mask_r", "cache_prefix",
                        "calc_window_climate_fixed", "build_daily_path", "sanitize_values",
                        "mask_raster", "safe_write_raster", "should_write", "OVERWRITE",
                        "NODATA_OUT", "NODATA_ABS_MAX",
                        "MIN_VALID_FRAC", "is_leap_year", "doy_to_date_noleap",
                        "PRECIP_DAILY_DIR", "PRECIP_DAILY_PATTERN",
                        "TA_DAILY_DIR", "TA_DAILY_PATTERN",
                        "SW_DAILY_DIR", "SW_DAILY_PATTERN"),
                  envir = environment())

    # 加载raster包到每个工作进程
    clusterEvalQ(cl, library(raster))

    # 并行计算3个早季变量
    early_list <- parLapply(cl, list(
      list(dir = PRECIP_DAILY_DIR, pattern = PRECIP_DAILY_PATTERN, var = "P_pre", neg = FALSE),
      list(dir = TA_DAILY_DIR, pattern = TA_DAILY_PATTERN, var = "T_pre", neg = TRUE),
      list(dir = SW_DAILY_DIR, pattern = SW_DAILY_PATTERN, var = "SW_pre", neg = FALSE)
    ), function(params) {
      calc_window_climate_fixed(year, sos_climatology_r, pos_climatology_r,
                                params$dir, params$pattern, mask_r,
                                cache_prefix[params$var], params$var, params$neg,
                                window_label = "SOSav, POSav")
    })

    stopCluster(cl)
    on.exit()  # 清除on.exit

    p_pre_r <- early_list[[1]]
    t_pre_r <- early_list[[2]]
    sw_pre_r <- early_list[[3]]

  } else {
    # 单核模式（原始串行计算）
    p_pre_r <- calc_window_climate_fixed(year, sos_climatology_r, pos_climatology_r,
                                         PRECIP_DAILY_DIR, PRECIP_DAILY_PATTERN,
                                         mask_r, cache_prefix["P_pre"], "P_pre",
                                         allow_negative = FALSE, window_label = "SOSav, POSav")

    t_pre_r <- calc_window_climate_fixed(year, sos_climatology_r, pos_climatology_r,
                                         TA_DAILY_DIR, TA_DAILY_PATTERN,
                                         mask_r, cache_prefix["T_pre"], "T_pre",
                                         allow_negative = TRUE, window_label = "SOSav, POSav")

    sw_pre_r <- calc_window_climate_fixed(year, sos_climatology_r, pos_climatology_r,
                                          SW_DAILY_DIR, SW_DAILY_PATTERN,
                                          mask_r, cache_prefix["SW_pre"], "SW_pre",
                                          allow_negative = FALSE, window_label = "SOSav, POSav")
  }

  # ===【EOS】计算晚季气候因子与GPP =====
  cat("  [2/7] 晚季气候因子（固定窗口[POSav, EOSav]）:\n")
  cat("  [3/7] 晚季GPP（固定窗口[POSav, EOSav]）:\n")

  # ===【并行化优化】生长季4个变量并行计算===
  if (parallel_inner && PARALLEL_CORES > 1) {
    cat(sprintf("    使用 %d 个核心并行计算...\n", PARALLEL_CORES))

    # 创建并行集群（Windows兼容）
    cl <- makeCluster(min(PARALLEL_CORES, 4))  # 最多4个变量
    on.exit(stopCluster(cl), add = TRUE)

    # 导出必要的变量和函数到集群
    clusterExport(cl, c("year", "pos_climatology_r", "eos_climatology_r", "mask_r", "cache_prefix",
                        "calc_window_climate_fixed", "build_daily_path", "sanitize_values",
                        "mask_raster", "safe_write_raster", "should_write", "OVERWRITE",
                        "NODATA_OUT", "NODATA_ABS_MAX",
                        "MIN_VALID_FRAC", "is_leap_year", "doy_to_date_noleap",
                        "PRECIP_DAILY_DIR", "PRECIP_DAILY_PATTERN",
                        "TA_DAILY_DIR", "TA_DAILY_PATTERN",
                        "SW_DAILY_DIR", "SW_DAILY_PATTERN",
                        "GPP_DAILY_DIR", "GPP_DAILY_PATTERN"),
                  envir = environment())

    # 加载raster包到每个工作进程
    clusterEvalQ(cl, library(raster))

    # 并行计算3个生长季气候变量（不含GPP，GPP从DECOMP_DIR读取）
    season_list <- parLapply(cl, list(
      list(dir = PRECIP_DAILY_DIR, pattern = PRECIP_DAILY_PATTERN, var = "P_season", neg = FALSE),
      list(dir = TA_DAILY_DIR, pattern = TA_DAILY_PATTERN, var = "T_season", neg = TRUE),
      list(dir = SW_DAILY_DIR, pattern = SW_DAILY_PATTERN, var = "SW_season", neg = FALSE)
    ), function(params) {
      calc_window_climate_fixed(year, pos_climatology_r, eos_climatology_r,
                                params$dir, params$pattern, mask_r,
                                cache_prefix[params$var], params$var, params$neg,
                                window_label = "POSav, EOSav")
    })

    stopCluster(cl)
    on.exit()  # 清除on.exit

    p_season_r <- season_list[[1]]
    t_season_r <- season_list[[2]]
    sw_season_r <- season_list[[3]]

  } else {
    # 单核模式（原始串行计算）
    p_season_r <- calc_window_climate_fixed(year, pos_climatology_r, eos_climatology_r,
                                            PRECIP_DAILY_DIR, PRECIP_DAILY_PATTERN,
                                            mask_r, cache_prefix["P_season"], "P_season",
                                            allow_negative = FALSE, window_label = "POSav, EOSav")

    t_season_r <- calc_window_climate_fixed(year, pos_climatology_r, eos_climatology_r,
                                            TA_DAILY_DIR, TA_DAILY_PATTERN,
                                            mask_r, cache_prefix["T_season"], "T_season",
                                            allow_negative = TRUE, window_label = "POSav, EOSav")

    sw_season_r <- calc_window_climate_fixed(year, pos_climatology_r, eos_climatology_r,
                                             SW_DAILY_DIR, SW_DAILY_PATTERN,
                                             mask_r, cache_prefix["SW_season"], "SW_season",
                                             allow_negative = FALSE, window_label = "POSav, EOSav")
  }

  # ===【关键修改】从DECOMP_DIR读取03c生成的Fixed_GPPrate（而非重新计算）===
  fixed_gpprate_file <- file.path(DECOMP_DIR, sprintf("Fixed_GPPrate_%d.tif", year))
  if (!file.exists(fixed_gpprate_file)) {
    cat(sprintf("  警告: Fixed_GPPrate文件不存在 (%s)，将返回NULL\n", basename(fixed_gpprate_file)))
    return(NULL)
  }
  cat(sprintf("  [读取] Fixed_GPPrate from DECOMP_DIR: %s\n", basename(fixed_gpprate_file)))
  fixed_gpprate_r <- raster(fixed_gpprate_file)
  fixed_gpprate_r <- set_nodata_if_missing(fixed_gpprate_r)
  fixed_gpprate_r <- mask_raster(fixed_gpprate_r, mask_r)

  cat("  ✓ 数据准备完成\n")

  # ===【修改3d】返回Fixed_Trate和Fixed_GPPrate===
  list(
    fixed_trate = fixed_trate_r,  # 修改：trc → fixed_trate
    eos = eos_r,
    p_pre = p_pre_r,
    t_pre = t_pre_r,
    sw_pre = sw_pre_r,
    p_season = p_season_r,
    t_season = t_season_r,
    sw_season = sw_season_r,
    fixed_gpprate = fixed_gpprate_r  # 修改：gpp_season → fixed_gpprate（从DECOMP_DIR读取）
  )
}

# ==================== 缓存准备 ====================
# ===【修改4】添加POSav和fixed_window_length_r参数===
prepare_sem_caches <- function(years, sos_climatology_r, pos_climatology_r,
                               eos_climatology_r, fixed_window_length_r, mask_r) {
  cat("\n=== 准备SEM缓存数据 ===\n")
  if (PARALLEL_ENABLE && PARALLEL_CORES > 1 && length(years) > 1) {
    cat(sprintf("  使用 %d 个核心并行准备缓存...\n", PARALLEL_CORES))
    cl <- makeCluster(PARALLEL_CORES)
    clusterEvalQ(cl, library(raster))
    clusterExport(
      cl,
      c("years", "sos_climatology_r", "pos_climatology_r", "eos_climatology_r",
        "fixed_window_length_r", "mask_r",
        "prepare_dual_timescale_data", "calc_window_climate_fixed",
        "build_daily_path", "sanitize_values", "mask_raster", "safe_write_raster",
        "should_write", "OVERWRITE", "set_nodata_if_missing", "PARALLEL_CORES",
        "NODATA_OUT", "NODATA_ABS_MAX", "MIN_VALID_FRAC", "is_leap_year",
        "doy_to_date_noleap", "DECOMP_DIR", "PHENO_DIR", "DERIVED_DIR",
        "GPP_DAILY_DIR", "GPP_DAILY_PATTERN",
        "PRECIP_DAILY_DIR", "PRECIP_DAILY_PATTERN",
        "TA_DAILY_DIR", "TA_DAILY_PATTERN",
        "SW_DAILY_DIR", "SW_DAILY_PATTERN"),
      envir = environment()
    )
    res_list <- parLapply(cl, years, function(year) {
      tryCatch({
        prepare_dual_timescale_data(year, sos_climatology_r, pos_climatology_r,
                                   eos_climatology_r,
                                   fixed_window_length_r, mask_r,
                                   parallel_inner = FALSE)
        list(year = year, ok = TRUE, msg = NULL)
      }, error = function(e) {
        list(year = year, ok = FALSE, msg = conditionMessage(e))
      })
    })
    stopCluster(cl)

    for (res in res_list) {
      cat(sprintf("  年份: %d\n", res$year))
      if (isTRUE(res$ok)) {
        cat("    ✓ 完成\n")
      } else {
        cat(sprintf("    ✗ 错误: %s\n", res$msg))
        cat(sprintf("    提示: 检查年份 %d 的输入数据是否完整\n", res$year))
      }
    }
  } else {
    for (year in years) {
      cat(sprintf("  年份: %d\n", year))
      tryCatch({
        prepare_dual_timescale_data(year, sos_climatology_r, pos_climatology_r,
                                   eos_climatology_r,
                                   fixed_window_length_r, mask_r)
        cat(sprintf("    ✓ 完成\n"))
      }, error = function(e) {
        cat(sprintf("    ✗ 错误: %s\n", conditionMessage(e)))
        cat(sprintf("    提示: 检查年份 %d 的输入数据是否完整\n", year))
      })
    }
  }
}

# ==================== 中介比例安全计算 ====================
# 防止分母接近0导致mediation ratio爆炸为±Inf
# 参考05c实现，当 |e + c*d| < eps 时返回NA
MEDIATION_DENOM_EPS <- 1e-6

safe_mediation_ratio <- function(c_val, d_val, e_val, eps = MEDIATION_DENOM_EPS) {
  denom <- e_val + c_val * d_val
  if (!is.finite(denom) || abs(denom) < eps) {
    return(NA_real_)
  }
  (c_val * d_val) / denom
}

# ==================== VIF诊断 ====================
# ===【修改5】Fixed_Trate替代TRc===
calculate_vif <- function(data) {
  # 双时间尺度模型的VIF诊断
  lm_model <- lm(Fixed_Trate ~ EOS + Fixed_GPPrate + P_pre + T_pre + SW_pre +
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
  safe_write_csv(vif_df, file.path(OUTPUT_DIR, "VIF_diagnostics.csv"), row.names = FALSE)

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
    TR_fixed_window = file.path(DECOMP_DIR, sprintf("TR_fixed_window_%d.tif", years)),  # 修改：TRc → TR_fixed_window
    EOS = file.path(PHENO_DIR, sprintf("eos_gpp_%d.tif", years)),
    Fixed_GPPrate = file.path(DECOMP_DIR, sprintf("Fixed_GPPrate_%d.tif", years)),  # 从DECOMP_DIR读取03c生成的Fixed_GPPrate
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
    "EOS~P_pre",
    "EOS~T_pre",
    "EOS~SW_pre",
    "EOS~P_season",
    "EOS~T_season",
    "EOS~SW_season",
    "Fixed_GPPrate~EOS",
    "Fixed_GPPrate~P_pre",
    "Fixed_GPPrate~T_pre",
    "Fixed_GPPrate~SW_pre",
    "Fixed_GPPrate~P_season",
    "Fixed_GPPrate~T_season",
    "Fixed_GPPrate~SW_season",
    "Fixed_Trate~EOS",
    "Fixed_Trate~Fixed_GPPrate",
    "Fixed_Trate~P_pre",
    "Fixed_Trate~T_pre",
    "Fixed_Trate~SW_pre",
    "Fixed_Trate~P_season",
    "Fixed_Trate~T_season",
    "Fixed_Trate~SW_season",
    "P_pre_via_EOS_GPP",
    "T_pre_via_EOS_GPP",
    "SW_pre_via_EOS_GPP",
    "P_pre_via_EOS",
    "T_pre_via_EOS",
    "SW_pre_via_EOS",
    "P_pre_via_GPP",
    "T_pre_via_GPP",
    "SW_pre_via_GPP",
    "P_season_via_EOS_GPP",
    "T_season_via_EOS_GPP",
    "SW_season_via_EOS_GPP",
    "P_season_via_EOS",
    "T_season_via_EOS",
    "SW_season_via_EOS",
    "P_season_via_GPP",
    "T_season_via_GPP",
    "SW_season_via_GPP",
    "EOS_via_GPP",
    "P_GPP_mediation",
    "T_GPP_mediation",
    "SW_GPP_mediation"
  )

  coef_sum <- setNames(rep(0, length(coef_names)), coef_names)
  coef_sumsq <- setNames(rep(0, length(coef_names)), coef_names)
  coef_count <- setNames(rep(0, length(coef_names)), coef_names)
  sig_sum <- setNames(rep(0, length(coef_names)), coef_names)
  sig_sumsq <- setNames(rep(0, length(coef_names)), coef_names)
  sig_count <- setNames(rep(0, length(coef_names)), coef_names)

  r2_names <- c("EOS", "Fixed_GPPrate", "Fixed_Trate")
  r2_sum <- setNames(rep(0, length(r2_names)), r2_names)
  r2_sumsq <- setNames(rep(0, length(r2_names)), r2_names)
  r2_count <- setNames(rep(0, length(r2_names)), r2_names)
  zero_coef <- setNames(rep(0, length(coef_names)), coef_names)
  zero_r2 <- setNames(rep(0, length(r2_names)), r2_names)
  zero_filter <- c(coef_extreme = 0, p_invalid = 0, r2_invalid = 0)
  filter_info <- zero_filter
  coef_samples_list <- list()

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

    eos_block <- t(getValues(stacks$EOS, row = row, nrows = nrows))
    eos_block <- sanitize_values(eos_block, na_values$EOS, allow_negative = FALSE)
    gpp_block <- t(getValues(stacks$Fixed_GPPrate, row = row, nrows = nrows))
    gpp_block <- sanitize_values(gpp_block, na_values$Fixed_GPPrate, allow_negative = TRUE)  # Fixed_GPPrate是异常值，可以为负
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

    valid_mat <- is.finite(fixed_trate_block) & is.finite(eos_block) & is.finite(gpp_block) &
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
    coef_samples_block <- matrix(NA_real_, nrow = length(valid_cells), ncol = length(coef_names))
    coef_sample_rows <- 0

    for (idx in valid_cells) {
      valid_years <- valid_mat[, idx]
      if (sum(valid_years) < min_years) {
        next
      }

      df <- data.frame(
        year = years[valid_years],
        Fixed_Trate = fixed_trate_block[valid_years, idx],
        EOS = eos_block[valid_years, idx],
        Fixed_GPPrate = gpp_block[valid_years, idx],
        P_pre = p_pre_block[valid_years, idx],
        T_pre = t_pre_block[valid_years, idx],
        SW_pre = sw_pre_block[valid_years, idx],
        P_season = p_season_block[valid_years, idx],
        T_season = t_season_block[valid_years, idx],
        SW_season = sw_season_block[valid_years, idx]
      )

      if (DETREND_PIXEL_ENABLE) {
        df$Fixed_Trate <- detrend_series(df$Fixed_Trate, df$year)
        df$EOS <- detrend_series(df$EOS, df$year)
        df$Fixed_GPPrate <- detrend_series(df$Fixed_GPPrate, df$year)
        df$P_pre <- detrend_series(df$P_pre, df$year)
        df$T_pre <- detrend_series(df$T_pre, df$year)
        df$SW_pre <- detrend_series(df$SW_pre, df$year)
        df$P_season <- detrend_series(df$P_season, df$year)
        df$T_season <- detrend_series(df$T_season, df$year)
        df$SW_season <- detrend_series(df$SW_season, df$year)
      }

      fixed_trate_z <- scale_vec(df$Fixed_Trate)
      eos_z <- scale_vec(df$EOS)
      gpp_z <- scale_vec(df$Fixed_GPPrate)
      p_pre_z <- scale_vec(df$P_pre)
      t_pre_z <- scale_vec(df$T_pre)
      sw_pre_z <- scale_vec(df$SW_pre)
      p_season_z <- scale_vec(df$P_season)
      t_season_z <- scale_vec(df$T_season)
      sw_season_z <- scale_vec(df$SW_season)

      if (any(is.na(c(fixed_trate_z, eos_z, gpp_z, p_pre_z, t_pre_z, sw_pre_z,
                      p_season_z, t_season_z, sw_season_z)))) {
        next
      }

      X_eos <- cbind(p_pre_z, t_pre_z, sw_pre_z, p_season_z, t_season_z, sw_season_z)
      a_res <- regress_beta_p(X_eos, eos_z)
      r2_eos <- calc_r2(X_eos, eos_z)
      a <- a_res$beta
      p_a <- a_res$p
      se_a <- a_res$se

      X_gpp <- cbind(eos_z, p_pre_z, t_pre_z, sw_pre_z, p_season_z, t_season_z, sw_season_z)
      bc_res <- regress_beta_p(X_gpp, gpp_z)
      r2_gpp <- calc_r2(X_gpp, gpp_z)
      bc <- bc_res$beta
      p_bc <- bc_res$p
      se_bc <- bc_res$se
      b <- bc[1]
      f <- bc[2:4]
      c <- bc[5:7]
      p_b <- p_bc[1]
      p_f <- p_bc[2:4]
      p_c <- p_bc[5:7]
      se_b <- se_bc[1]
      se_f <- se_bc[2:4]
      se_c <- se_bc[5:7]

      vars_tr <- c("EOS", "Fixed_GPPrate", "P_pre", "T_pre", "SW_pre",
                   "P_season", "T_season", "SW_season")
      X_tr_full <- cbind(eos_z, gpp_z, p_pre_z, t_pre_z, sw_pre_z,
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
      se_d <- rep(NA_real_, length(vars_tr))

      if (ncol(X_tr) >= 1) {
        d_res <- regress_beta_p(X_tr, fixed_trate_z)
        r2_tr <- calc_r2(X_tr, fixed_trate_z)
        if (!any(is.na(d_res$beta)) && !any(is.na(d_res$p))) {
          for (j in seq_along(vars_current)) {
            idx_match <- match(vars_current[j], vars_tr)
            d[idx_match] <- d_res$beta[j]
            p_d[idx_match] <- d_res$p[j]
            se_d[idx_match] <- d_res$se[j]
          }
        }
      } else {
        r2_tr <- NA_real_
      }

      if (any(is.na(c(a, b, c, f))) || any(is.na(c(p_a, p_b, p_c, p_f))) || all(is.na(d))) {
        next
      }

      vals <- c(
        a[1], a[2], a[3], a[4], a[5], a[6],
        b,
        f[1], f[2], f[3],
        c[1], c[2], c[3],
        d[1], d[2], d[3], d[4], d[5], d[6], d[7], d[8],
        a[1] * b * d[2],
        a[2] * b * d[2],
        a[3] * b * d[2],
        a[1] * d[1],
        a[2] * d[1],
        a[3] * d[1],
        f[1] * d[2],
        f[2] * d[2],
        f[3] * d[2],
        a[4] * b * d[2],
        a[5] * b * d[2],
        a[6] * b * d[2],
        a[4] * d[1],
        a[5] * d[1],
        a[6] * d[1],
        c[1] * d[2],
        c[2] * d[2],
        c[3] * d[2],
        b * d[2],
        safe_mediation_ratio(c[1], d[2], d[6]),
        safe_mediation_ratio(c[2], d[2], d[7]),
        safe_mediation_ratio(c[3], d[2], d[8])
      )

      p_vals <- c(
        p_a[1], p_a[2], p_a[3], p_a[4], p_a[5], p_a[6],
        p_b,
        p_f[1], p_f[2], p_f[3],
        p_c[1], p_c[2], p_c[3],
        p_d[1], p_d[2], p_d[3], p_d[4], p_d[5], p_d[6], p_d[7], p_d[8],
        delta_p_three(a[1], b, d[2], se_a[1], se_b, se_d[2]),
        delta_p_three(a[2], b, d[2], se_a[2], se_b, se_d[2]),
        delta_p_three(a[3], b, d[2], se_a[3], se_b, se_d[2]),
        delta_p_two(a[1], d[1], se_a[1], se_d[1]),
        delta_p_two(a[2], d[1], se_a[2], se_d[1]),
        delta_p_two(a[3], d[1], se_a[3], se_d[1]),
        delta_p_two(f[1], d[2], se_f[1], se_d[2]),
        delta_p_two(f[2], d[2], se_f[2], se_d[2]),
        delta_p_two(f[3], d[2], se_f[3], se_d[2]),
        delta_p_three(a[4], b, d[2], se_a[4], se_b, se_d[2]),
        delta_p_three(a[5], b, d[2], se_a[5], se_b, se_d[2]),
        delta_p_three(a[6], b, d[2], se_a[6], se_b, se_d[2]),
        delta_p_two(a[4], d[1], se_a[4], se_d[1]),
        delta_p_two(a[5], d[1], se_a[5], se_d[1]),
        delta_p_two(a[6], d[1], se_a[6], se_d[1]),
        delta_p_two(c[1], d[2], se_c[1], se_d[2]),
        delta_p_two(c[2], d[2], se_c[2], se_d[2]),
        delta_p_two(c[3], d[2], se_c[3], se_d[2]),
        delta_p_two(b, d[2], se_b, se_d[2]),
        delta_p_ratio(c[1], d[2], d[6], se_c[1], se_d[2], se_d[6]),
        delta_p_ratio(c[2], d[2], d[7], se_c[2], se_d[2], se_d[7]),
        delta_p_ratio(c[3], d[2], d[8], se_c[3], se_d[2], se_d[8])
      )

      r2_vals <- c(r2_eos, r2_gpp, r2_tr)
      filter_res <- filter_sem_outputs(vals, p_vals, r2_vals)
      vals <- filter_res$vals
      p_vals <- filter_res$p_vals
      r2_vals <- filter_res$r2_vals
      filter_info_block <- filter_info_block + filter_res$info

      coef_sample_rows <- coef_sample_rows + 1
      coef_samples_block[coef_sample_rows, ] <- vals

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
      coef_samples = if (coef_sample_rows > 0) {
        coef_samples_block[seq_len(coef_sample_rows), , drop = FALSE]
      } else {
        NULL
      },
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
    if (!is.null(res$coef_samples)) {
      coef_samples_list[[length(coef_samples_list) + 1]] <<- res$coef_samples
    }
    processed_cells <<- processed_cells + res$processed_cells
    blocks_done <<- blocks_done + 1
  }

  if (PARALLEL_ENABLE && PARALLEL_CORES > 1) {
    cat(sprintf("  启用并行块处理: %d cores\n", PARALLEL_CORES))
    cl <- makeCluster(PARALLEL_CORES)
    clusterEvalQ(cl, library(raster))
    clusterExport(
      cl,
      c("bs", "mask_r", "fixed_window_length_r", "stacks", "na_values", "years",
        "min_years", "DETREND_PIXEL_ENABLE", "coef_names", "r2_names",
        "sanitize_values", "calc_r2", "regress_beta_p", "calc_vif_vector",
        "detrend_series", "scale_vec", "zero_coef", "zero_r2", "process_block",
        "delta_p_two", "delta_p_three", "delta_p_ratio", "safe_mediation_ratio",
        "MEDIATION_DENOM_EPS",
        "filter_sem_outputs", "FILTER_SEM_OUTLIERS", "SEM_COEF_ABS_MAX",
        "SEM_P_MIN", "SEM_P_MAX", "SEM_R2_MIN", "SEM_R2_MAX",
        "zero_filter", "NODATA_OUT", "NODATA_ABS_MAX"),
      envir = environment()
    )

    chunk_size <- max(1, PARALLEL_CHUNK_SIZE)
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

  coef_samples_all <- NULL
  if (length(coef_samples_list) > 0) {
    coef_samples_all <- do.call(rbind, coef_samples_list)
  }

  boot_ci_low <- rep(NA_real_, length(coef_names))
  boot_ci_high <- rep(NA_real_, length(coef_names))
  boot_sig <- rep(NA, length(coef_names))

  if (PIXELWISE_BOOTSTRAP_ENABLE && !is.null(coef_samples_all)) {
    cat(sprintf("\n  像元均值 bootstrap CI: %d 次\n", PIXELWISE_BOOTSTRAP_N))
    for (k in seq_along(coef_names)) {
      ci <- bootstrap_ci_mean(
        coef_samples_all[, k],
        PIXELWISE_BOOTSTRAP_N,
        seed = PIXELWISE_BOOTSTRAP_SEED + k
      )
      boot_ci_low[k] <- ci["ci_low"]
      boot_ci_high[k] <- ci["ci_high"]
      if (is.finite(boot_ci_low[k]) && is.finite(boot_ci_high[k])) {
        boot_sig[k] <- !(boot_ci_low[k] <= 0 && boot_ci_high[k] >= 0)
      }
    }
  }

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
  safe_write_csv(pixel_filter_df,
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
    sd_sig = sd_sig,
    boot_ci_low = boot_ci_low,
    boot_ci_high = boot_ci_high,
    boot_sig = boot_sig
  )

  safe_write_csv(summary_df, file.path(PIXELWISE_DIR, "SEM_dual_timescale_parameters.csv"),
                 row.names = FALSE)

  boot_df <- data.frame(
    path = coef_names,
    boot_ci_low = boot_ci_low,
    boot_ci_high = boot_ci_high,
    boot_sig = boot_sig
  )
  safe_write_csv(boot_df,
                 file.path(PIXELWISE_DIR, "SEM_dual_timescale_pixelwise_bootstrap_ci.csv"),
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

  safe_write_csv(r2_summary, file.path(PIXELWISE_DIR, "SEM_dual_timescale_R2_detail.csv"),
                 row.names = FALSE)
  safe_write_csv(r2_summary[, c("variable", "R2")],
                 file.path(PIXELWISE_DIR, "SEM_dual_timescale_R2.csv"),
                 row.names = FALSE)

  summary_path <- file.path(PIXELWISE_DIR, "SEM_dual_timescale_summary.txt")
  if (should_write(summary_path)) {
    sink(summary_path)
    cat("像元时间序列SEM汇总\n\n")
    print(summary_df)
    if (FILTER_SEM_OUTLIERS) {
      cat("\n输出异常值过滤统计（像元级）:\n")
      print(pixel_filter_df)
    }
    if (PIXELWISE_BOOTSTRAP_ENABLE) {
      cat("\n像元均值Bootstrap CI:\n")
      print(boot_df)
    }
    cat("\nR2汇总:\n")
    print(r2_summary)
    sink()
  } else {
    cat(sprintf("  [skip] %s\n", summary_path))
  }

  cat("\n像元时间序列SEM摘要:\n")
  print(summary_df)

  summary_df
}

# ==================== SEM分析（全局合并）====================
# ===【修改7】添加POSav气候态计算===
run_dual_timescale_sem <- function(years, mask_r) {
  cat("\n=== 计算物候气候态（像元多年平均：SOSav, POSav, EOSav） ===\n")

  # 读取所有年份的SOS、POS、EOS数据并计算多年平均
  sos_stack_list <- list()
  pos_stack_list <- list()
  eos_stack_list <- list()
  n_loaded_sos <- 0
  n_loaded_pos <- 0
  n_loaded_eos <- 0

  for (year in years) {
    sos_file <- file.path(PHENO_DIR, sprintf("sos_gpp_%d.tif", year))
    pos_file <- file.path(PHENO_DIR, sprintf("pos_doy_gpp_%d.tif", year))
    eos_file <- file.path(PHENO_DIR, sprintf("eos_gpp_%d.tif", year))

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

    if (file.exists(eos_file)) {
      eos_r <- raster(eos_file)
      eos_r <- set_nodata_if_missing(eos_r)
      eos_r <- mask_raster(eos_r, mask_r)
      eos_stack_list[[length(eos_stack_list) + 1]] <- eos_r
      n_loaded_eos <- n_loaded_eos + 1
    }

    if ((n_loaded_sos %% 10 == 0) && (n_loaded_sos > 0)) {
      cat(sprintf("  已读取 %d/%d 年份的物候数据...\n", n_loaded_sos, length(years)))
    }
  }

  cat(sprintf("  ✓ 成功读取 %d/%d 年份的SOS数据\n", n_loaded_sos, length(years)))
  cat(sprintf("  ✓ 成功读取 %d/%d 年份的POS数据\n", n_loaded_pos, length(years)))
  cat(sprintf("  ✓ 成功读取 %d/%d 年份的EOS数据\n", n_loaded_eos, length(years)))

  if (n_loaded_sos == 0 || n_loaded_pos == 0 || n_loaded_eos == 0) {
    stop("未找到足够的SOS/POS/EOS数据文件")
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

  cat("  正在计算EOS多年平均（EOSav）...\n")
  eos_stack <- stack(eos_stack_list)
  eos_climatology_r <- calc(eos_stack, fun = function(x) {
    mean(x, na.rm = TRUE)
  })
  eos_climatology_r <- mask_raster(eos_climatology_r, mask_r)

  # ===【修改：先统计掩膜内有效像元，再验证POS > SOS、EOS > POS】===
  mask_vals <- getValues(mask_r)
  in_mask <- !is.na(mask_vals) & (mask_vals > 0)

  sos_clim_vals_pre <- getValues(sos_climatology_r)
  sos_clim_vals_pre <- sanitize_values(sos_clim_vals_pre, NAvalue(sos_climatology_r), FALSE)
  pos_clim_vals_pre <- getValues(pos_climatology_r)
  pos_clim_vals_pre <- sanitize_values(pos_clim_vals_pre, NAvalue(pos_climatology_r), FALSE)
  eos_clim_vals_pre <- getValues(eos_climatology_r)
  eos_clim_vals_pre <- sanitize_values(eos_clim_vals_pre, NAvalue(eos_climatology_r), FALSE)

  valid_sos_pre <- in_mask & !is.na(sos_clim_vals_pre)
  valid_pos_pre <- in_mask & !is.na(pos_clim_vals_pre)
  valid_eos_pre <- in_mask & !is.na(eos_clim_vals_pre)
  valid_overlap_pre <- valid_sos_pre & valid_pos_pre & valid_eos_pre

  cat(sprintf("  掩膜内有效像元(预过滤): SOSav=%d, POSav=%d, EOSav=%d, 重叠=%d\n",
              sum(valid_sos_pre), sum(valid_pos_pre), sum(valid_eos_pre), sum(valid_overlap_pre)))

  cat("  验证POSav > SOSav, EOSav > POSav...\n")
  valid_window <- valid_overlap_pre &
    (pos_clim_vals_pre > sos_clim_vals_pre) &
    (eos_clim_vals_pre > pos_clim_vals_pre)
  n_invalid <- sum(valid_overlap_pre & !valid_window)
  cat(sprintf("  POSav/EOSav窗口过滤掉的像元数: %d\n", n_invalid))
  if (n_invalid > 0) {
    cat(sprintf("  ⚠️ 过滤了 %d 个无效像元\n", n_invalid))
  }

  valid_window_mask <- setValues(raster(mask_r), as.integer(valid_window))
  sos_climatology_r <- mask(sos_climatology_r, valid_window_mask, maskvalue = 0)
  pos_climatology_r <- mask(pos_climatology_r, valid_window_mask, maskvalue = 0)
  eos_climatology_r <- mask(eos_climatology_r, valid_window_mask, maskvalue = 0)

  # 统计气候态信息（过滤后）
  sos_clim_vals <- getValues(sos_climatology_r)
  sos_clim_vals <- sanitize_values(sos_clim_vals, NAvalue(sos_climatology_r), FALSE)
  pos_clim_vals <- getValues(pos_climatology_r)
  pos_clim_vals <- sanitize_values(pos_clim_vals, NAvalue(pos_climatology_r), FALSE)
  eos_clim_vals <- getValues(eos_climatology_r)
  eos_clim_vals <- sanitize_values(eos_clim_vals, NAvalue(eos_climatology_r), FALSE)

  n_valid_sos <- sum(!is.na(sos_clim_vals))
  n_valid_pos <- sum(!is.na(pos_clim_vals))
  n_valid_eos <- sum(!is.na(eos_clim_vals))
  sos_min <- min(sos_clim_vals, na.rm = TRUE)
  sos_max <- max(sos_clim_vals, na.rm = TRUE)
  sos_mean <- mean(sos_clim_vals, na.rm = TRUE)
  pos_min <- min(pos_clim_vals, na.rm = TRUE)
  pos_max <- max(pos_clim_vals, na.rm = TRUE)
  pos_mean <- mean(pos_clim_vals, na.rm = TRUE)
  eos_min <- min(eos_clim_vals, na.rm = TRUE)
  eos_max <- max(eos_clim_vals, na.rm = TRUE)
  eos_mean <- mean(eos_clim_vals, na.rm = TRUE)

  cat(sprintf("  ✓ 物候气候态计算完成:\n"))
  cat(sprintf("    SOSav - 有效像元: %d, 范围: %.1f - %.1f DOY, 平均: %.1f DOY\n",
              n_valid_sos, sos_min, sos_max, sos_mean))
  cat(sprintf("    POSav - 有效像元: %d, 范围: %.1f - %.1f DOY, 平均: %.1f DOY\n",
              n_valid_pos, pos_min, pos_max, pos_mean))
  cat(sprintf("    EOSav - 有效像元: %d, 范围: %.1f - %.1f DOY, 平均: %.1f DOY\n",
              n_valid_eos, eos_min, eos_max, eos_mean))

  # 保存物候气候态
  sos_clim_file <- file.path(DERIVED_DIR, "SOS_climatology.tif")
  pos_clim_file <- file.path(DERIVED_DIR, "POS_climatology.tif")
  eos_clim_file <- file.path(DERIVED_DIR, "EOS_climatology.tif")
  safe_write_raster(sos_climatology_r, sos_clim_file, datatype = "FLT4S")
  safe_write_raster(pos_climatology_r, pos_clim_file, datatype = "FLT4S")
  safe_write_raster(eos_climatology_r, eos_clim_file, datatype = "FLT4S")
  cat(sprintf("  ✓ 物候气候态已保存: %s, %s, %s\n",
              basename(sos_clim_file), basename(pos_clim_file), basename(eos_clim_file)))

  # 读取Fixed_Window_Length（用于计算Fixed_Trate）
  cat("\n=== 读取Fixed_Window_Length ===\n")
  fixed_len_file <- file.path(DECOMP_DIR, "Fixed_Window_Length.tif")
  if (!file.exists(fixed_len_file)) {
    stop(sprintf("Fixed_Window_Length文件不存在: %s", fixed_len_file))
  }
  fixed_window_length_r <- raster(fixed_len_file)
  fixed_window_length_r <- set_nodata_if_missing(fixed_window_length_r)
  fixed_window_length_r <- mask_raster(fixed_window_length_r, mask_r)
  fixed_window_length_r[fixed_window_length_r <= 0] <- NA
  cat(sprintf("  ✓ Fixed_Window_Length已加载: %s\n", basename(fixed_len_file)))

  cat("\n=== 准备SEM数据 ===\n")

  sem_data <- data.frame()

  # ===【修改8】调用时传入POSav和fixed_window_length_r参数，数据框列名改为Fixed_Trate===
  if (SEM_SAMPLE_MODE == "annual_mean" && PARALLEL_ENABLE && PARALLEL_CORES > 1 && length(years) > 1) {
    cat(sprintf("  使用 %d 个核心并行准备全域均值数据...\n", PARALLEL_CORES))
    cl <- makeCluster(PARALLEL_CORES)
    clusterEvalQ(cl, library(raster))
    clusterExport(
      cl,
      c("years", "sos_climatology_r", "pos_climatology_r", "eos_climatology_r",
        "fixed_window_length_r", "mask_r",
        "prepare_dual_timescale_data", "calc_window_climate_fixed",
        "build_daily_path", "sanitize_values", "mask_raster", "safe_write_raster",
        "should_write", "OVERWRITE", "set_nodata_if_missing", "PARALLEL_CORES",
        "NODATA_OUT", "NODATA_ABS_MAX", "MIN_VALID_FRAC", "is_leap_year",
        "doy_to_date_noleap", "DECOMP_DIR", "PHENO_DIR", "DERIVED_DIR",
        "GPP_DAILY_DIR", "GPP_DAILY_PATTERN",
        "PRECIP_DAILY_DIR", "PRECIP_DAILY_PATTERN",
        "TA_DAILY_DIR", "TA_DAILY_PATTERN",
        "SW_DAILY_DIR", "SW_DAILY_PATTERN"),
      envir = environment()
    )
    res_list <- parLapply(cl, years, function(year) {
      tryCatch({
        rasters <- prepare_dual_timescale_data(year, sos_climatology_r, pos_climatology_r,
                                               eos_climatology_r,
                                               fixed_window_length_r, mask_r,
                                               parallel_inner = FALSE)
        if (is.null(rasters)) {
          return(list(year = year, ok = FALSE, msg = "输出为空", row = NULL))
        }
        year_row <- data.frame(
          year = year,
          Fixed_Trate = cellStats(rasters$fixed_trate, mean, na.rm = TRUE),
          EOS = cellStats(rasters$eos, mean, na.rm = TRUE),
          Fixed_GPPrate = cellStats(rasters$fixed_gpprate, mean, na.rm = TRUE),
          P_pre = cellStats(rasters$p_pre, mean, na.rm = TRUE),
          T_pre = cellStats(rasters$t_pre, mean, na.rm = TRUE),
          SW_pre = cellStats(rasters$sw_pre, mean, na.rm = TRUE),
          P_season = cellStats(rasters$p_season, mean, na.rm = TRUE),
          T_season = cellStats(rasters$t_season, mean, na.rm = TRUE),
          SW_season = cellStats(rasters$sw_season, mean, na.rm = TRUE)
        )
        list(year = year, ok = TRUE, msg = NULL, row = year_row)
      }, error = function(e) {
        list(year = year, ok = FALSE, msg = conditionMessage(e), row = NULL)
      })
    })
    stopCluster(cl)

    for (res in res_list) {
      if (isTRUE(res$ok)) {
        sem_data <- rbind(sem_data, res$row)
        cat(sprintf("  年份 %d: 全域均值已添加\n", res$year))
      } else {
        cat(sprintf("  年份 %d: 跳过（%s）\n", res$year, res$msg))
      }
    }
  } else {
    for (year in years) {
      rasters <- prepare_dual_timescale_data(year, sos_climatology_r, pos_climatology_r,
                                             eos_climatology_r,
                                             fixed_window_length_r, mask_r)
      if (is.null(rasters)) next

      if (SEM_SAMPLE_MODE == "annual_mean") {
        year_row <- data.frame(
          year = year,
          Fixed_Trate = cellStats(rasters$fixed_trate, mean, na.rm = TRUE),
          EOS = cellStats(rasters$eos, mean, na.rm = TRUE),
          Fixed_GPPrate = cellStats(rasters$fixed_gpprate, mean, na.rm = TRUE),
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
        vals <- data.frame(
          year = year,
          Fixed_Trate = sanitize_values(getValues(rasters$fixed_trate), NAvalue(rasters$fixed_trate), TRUE),
          EOS = sanitize_values(getValues(rasters$eos), NAvalue(rasters$eos), FALSE),
          Fixed_GPPrate = sanitize_values(getValues(rasters$fixed_gpprate), NAvalue(rasters$fixed_gpprate), FALSE),
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
  }

  # 数据质量检查
  cat("\n数据质量检查:\n")
  cat(sprintf("  总样本数: %d\n", nrow(sem_data)))
  cat(sprintf("  完整案例数: %d\n", sum(complete.cases(sem_data))))

  # 修改：TRc → Fixed_Trate
  for (col in c("Fixed_Trate", "EOS", "Fixed_GPPrate", "P_pre", "T_pre", "SW_pre",
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
      safe_write_csv(sem_data_for_sem, file.path(DATA_DIR, "sem_dual_timescale_detrended.csv"), row.names = FALSE)
      cat("✓ 去趋势完成: sem_dual_timescale_detrended.csv\n")
    }
  }

  # Z-score标准化
  sem_data_std <- sem_data_for_sem
  sem_vars <- setdiff(names(sem_data_std), "year")
  sem_data_std[sem_vars] <- scale(sem_data_std[sem_vars])

  # 保存原始和标准化数据
  safe_write_csv(sem_data, file.path(DATA_DIR, "sem_dual_timescale_raw.csv"), row.names = FALSE)
  safe_write_csv(sem_data_std, file.path(DATA_DIR, "sem_dual_timescale_standardized.csv"), row.names = FALSE)

  cat("✓ SEM数据准备完成\n")

  # ===== 构建双时间尺度SEM模型 =====
  # ===【修改9】SEM模型：Fixed_Trate替代TRc===
  cat("\n=== 构建双时间尺度SEM模型（固定窗口方法） ===\n")
  cat("✓ 气候因子: 降水(P)、气温(T)、短波辐射(SW)\n")
  cat("✓ 时间窗口: 早季气候（SOSav-POSav）+ 晚季气候（POSav-EOSav）\n")
  cat("✓ 因变量: Fixed_Trate（固定窗口蒸腾速率）\n")

  # 完整路径模型（降水 + 温度 + 短波辐射）
  sem_model <- '
    # 第一层：早/晚季气候 → EOS（物候响应）
    EOS ~ a1*P_pre + a2*T_pre + a3*SW_pre +
          a4*P_season + a5*T_season + a6*SW_season

    # 第二层：EOS + 早/晚季气候 → Fixed_GPPrate（碳固定）
    Fixed_GPPrate ~ b*EOS + f1*P_pre + f2*T_pre + f3*SW_pre +
                 c1*P_season + c2*T_season + c3*SW_season

    # 第三层：EOS + Fixed_GPPrate + 早/晚季气候 → Fixed_Trate（固定窗口蒸腾速率）
    Fixed_Trate ~ g*EOS + d*Fixed_GPPrate +
          h1*P_pre + h2*T_pre + h3*SW_pre +
          e1*P_season + e2*T_season + e3*SW_season

    # === 间接效应分解 ===

    # 早季气候的间接效应路径
    P_pre_via_EOS_GPP  := a1 * b * d
    T_pre_via_EOS_GPP  := a2 * b * d
    SW_pre_via_EOS_GPP := a3 * b * d

    P_pre_via_EOS  := a1 * g
    T_pre_via_EOS  := a2 * g
    SW_pre_via_EOS := a3 * g

    P_pre_via_GPP  := f1 * d
    T_pre_via_GPP  := f2 * d
    SW_pre_via_GPP := f3 * d

    # 早季气候的总间接效应
    P_pre_indirect  := a1*b*d + a1*g + f1*d
    T_pre_indirect  := a2*b*d + a2*g + f2*d
    SW_pre_indirect := a3*b*d + a3*g + f3*d

    # 晚季气候经由EOS的间接效应
    P_season_via_EOS_GPP  := a4 * b * d
    T_season_via_EOS_GPP  := a5 * b * d
    SW_season_via_EOS_GPP := a6 * b * d

    P_season_via_EOS  := a4 * g
    T_season_via_EOS  := a5 * g
    SW_season_via_EOS := a6 * g

    # 晚季气候通过GPP的间接效应
    P_season_via_GPP  := c1 * d
    T_season_via_GPP  := c2 * d
    SW_season_via_GPP := c3 * d

    # 晚季气候的总间接效应
    P_season_indirect  := a4*b*d + a4*g + c1*d
    T_season_indirect  := a5*b*d + a5*g + c2*d
    SW_season_indirect := a6*b*d + a6*g + c3*d

    # EOS通过GPP的间接效应
    EOS_via_GPP := b * d

    # GPP的中介比例
    P_GPP_mediation := (c1*d) / (e1 + c1*d)
    T_GPP_mediation := (c2*d) / (e2 + c2*d)
    SW_GPP_mediation := (c3*d) / (e3 + c3*d)
  '

  cat("\nSEM模型结构:\n")
  cat(sem_model)
  cat("\n")

  # ===== 拟合SEM模型 =====
  cat("\n=== 拟合SEM模型 ===\n")

  fit <- sem(sem_model, data = sem_data_std, estimator = "MLR")

  fit_summary <- summary(fit, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

  # 保存结果
  params <- parameterEstimates(fit, standardized = TRUE)
  param_filter <- filter_sem_param_table(params)
  params <- param_filter$params
  if (FILTER_SEM_OUTLIERS) {
    cat(sprintf("\n  输出异常值过滤（全域SEM）：系数极端=%d, p值异常=%d\n",
                param_filter$info["coef_extreme"], param_filter$info["p_invalid"]))
  }

  # ===【Bug Fix - Problem 1】中介比例分母保护===
  # 修复：当分母 e + c*d ≈ 0 时，mediation ratio会爆炸为±Inf
  # 解决方案：使用safe_mediation_ratio()重新计算，分母过小时置NA
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
    if (length(idx) == 0) {
      idx <- which(params_df$lhs == label & params_df$op == ":=")
    }
    if (length(idx) > 0) {
      params_df$est[idx] <- value
      params_df$std.all[idx] <- value
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

  update_summary_ratio <- function(sum_obj, label, value) {
    idx <- which(sum_obj$pe$label == label)
    if (length(idx) == 0) {
      idx <- which(sum_obj$pe$lhs == label & sum_obj$pe$op == ":=")
    }
    if (length(idx) > 0) {
      sum_obj$pe$est[idx] <- value
      sum_obj$pe$std.all[idx] <- value
    }
    sum_obj
  }
  fit_summary <- update_summary_ratio(fit_summary, "P_GPP_mediation", p_ratio)
  fit_summary <- update_summary_ratio(fit_summary, "T_GPP_mediation", t_ratio)
  fit_summary <- update_summary_ratio(fit_summary, "SW_GPP_mediation", sw_ratio)

  get_label_stats <- function(params_df, label) {
    idx <- which(params_df$label == label)
    if (length(idx) == 0) {
      return(list(std = NA_real_, se_std = NA_real_))
    }
    est <- params_df$est[idx][1]
    std <- params_df$std.all[idx][1]
    se <- params_df$se[idx][1]
    se_std <- if (is.finite(est) && est != 0 && is.finite(std) && is.finite(se)) {
      abs(std / est) * se
    } else {
      NA_real_
    }
    list(std = std, se_std = se_std)
  }

  update_defined_param <- function(params_df, name, std_val, p_val) {
    idx <- which(params_df$lhs == name & params_df$op == ":=")
    if (length(idx) > 0) {
      params_df$std.all[idx] <- std_val
      params_df$pvalue[idx] <- p_val
    }
    params_df
  }

  a1_stats <- get_label_stats(params, "a1")
  a2_stats <- get_label_stats(params, "a2")
  a3_stats <- get_label_stats(params, "a3")
  a4_stats <- get_label_stats(params, "a4")
  a5_stats <- get_label_stats(params, "a5")
  a6_stats <- get_label_stats(params, "a6")
  b_stats <- get_label_stats(params, "b")
  f1_stats <- get_label_stats(params, "f1")
  f2_stats <- get_label_stats(params, "f2")
  f3_stats <- get_label_stats(params, "f3")
  c1_stats <- get_label_stats(params, "c1")
  c2_stats <- get_label_stats(params, "c2")
  c3_stats <- get_label_stats(params, "c3")
  g_stats <- get_label_stats(params, "g")
  d_stats <- get_label_stats(params, "d")
  e1_stats <- get_label_stats(params, "e1")
  e2_stats <- get_label_stats(params, "e2")
  e3_stats <- get_label_stats(params, "e3")

  p_ratio_std <- safe_mediation_ratio(c1_stats$std, d_stats$std, e1_stats$std)
  t_ratio_std <- safe_mediation_ratio(c2_stats$std, d_stats$std, e2_stats$std)
  sw_ratio_std <- safe_mediation_ratio(c3_stats$std, d_stats$std, e3_stats$std)

  params <- update_defined_param(params, "P_pre_via_EOS_GPP",
                                 a1_stats$std * b_stats$std * d_stats$std,
                                 delta_p_three(a1_stats$std, b_stats$std, d_stats$std,
                                               a1_stats$se_std, b_stats$se_std, d_stats$se_std))
  params <- update_defined_param(params, "T_pre_via_EOS_GPP",
                                 a2_stats$std * b_stats$std * d_stats$std,
                                 delta_p_three(a2_stats$std, b_stats$std, d_stats$std,
                                               a2_stats$se_std, b_stats$se_std, d_stats$se_std))
  params <- update_defined_param(params, "SW_pre_via_EOS_GPP",
                                 a3_stats$std * b_stats$std * d_stats$std,
                                 delta_p_three(a3_stats$std, b_stats$std, d_stats$std,
                                               a3_stats$se_std, b_stats$se_std, d_stats$se_std))
  params <- update_defined_param(params, "P_pre_via_EOS",
                                 a1_stats$std * g_stats$std,
                                 delta_p_two(a1_stats$std, g_stats$std,
                                             a1_stats$se_std, g_stats$se_std))
  params <- update_defined_param(params, "T_pre_via_EOS",
                                 a2_stats$std * g_stats$std,
                                 delta_p_two(a2_stats$std, g_stats$std,
                                             a2_stats$se_std, g_stats$se_std))
  params <- update_defined_param(params, "SW_pre_via_EOS",
                                 a3_stats$std * g_stats$std,
                                 delta_p_two(a3_stats$std, g_stats$std,
                                             a3_stats$se_std, g_stats$se_std))
  params <- update_defined_param(params, "P_pre_via_GPP",
                                 f1_stats$std * d_stats$std,
                                 delta_p_two(f1_stats$std, d_stats$std,
                                             f1_stats$se_std, d_stats$se_std))
  params <- update_defined_param(params, "T_pre_via_GPP",
                                 f2_stats$std * d_stats$std,
                                 delta_p_two(f2_stats$std, d_stats$std,
                                             f2_stats$se_std, d_stats$se_std))
  params <- update_defined_param(params, "SW_pre_via_GPP",
                                 f3_stats$std * d_stats$std,
                                 delta_p_two(f3_stats$std, d_stats$std,
                                             f3_stats$se_std, d_stats$se_std))
  params <- update_defined_param(params, "P_season_via_EOS_GPP",
                                 a4_stats$std * b_stats$std * d_stats$std,
                                 delta_p_three(a4_stats$std, b_stats$std, d_stats$std,
                                               a4_stats$se_std, b_stats$se_std, d_stats$se_std))
  params <- update_defined_param(params, "T_season_via_EOS_GPP",
                                 a5_stats$std * b_stats$std * d_stats$std,
                                 delta_p_three(a5_stats$std, b_stats$std, d_stats$std,
                                               a5_stats$se_std, b_stats$se_std, d_stats$se_std))
  params <- update_defined_param(params, "SW_season_via_EOS_GPP",
                                 a6_stats$std * b_stats$std * d_stats$std,
                                 delta_p_three(a6_stats$std, b_stats$std, d_stats$std,
                                               a6_stats$se_std, b_stats$se_std, d_stats$se_std))
  params <- update_defined_param(params, "P_season_via_EOS",
                                 a4_stats$std * g_stats$std,
                                 delta_p_two(a4_stats$std, g_stats$std,
                                             a4_stats$se_std, g_stats$se_std))
  params <- update_defined_param(params, "T_season_via_EOS",
                                 a5_stats$std * g_stats$std,
                                 delta_p_two(a5_stats$std, g_stats$std,
                                             a5_stats$se_std, g_stats$se_std))
  params <- update_defined_param(params, "SW_season_via_EOS",
                                 a6_stats$std * g_stats$std,
                                 delta_p_two(a6_stats$std, g_stats$std,
                                             a6_stats$se_std, g_stats$se_std))
  params <- update_defined_param(params, "P_season_via_GPP",
                                 c1_stats$std * d_stats$std,
                                 delta_p_two(c1_stats$std, d_stats$std,
                                             c1_stats$se_std, d_stats$se_std))
  params <- update_defined_param(params, "T_season_via_GPP",
                                 c2_stats$std * d_stats$std,
                                 delta_p_two(c2_stats$std, d_stats$std,
                                             c2_stats$se_std, d_stats$se_std))
  params <- update_defined_param(params, "SW_season_via_GPP",
                                 c3_stats$std * d_stats$std,
                                 delta_p_two(c3_stats$std, d_stats$std,
                                             c3_stats$se_std, d_stats$se_std))
  params <- update_defined_param(params, "EOS_via_GPP",
                                 b_stats$std * d_stats$std,
                                 delta_p_two(b_stats$std, d_stats$std,
                                             b_stats$se_std, d_stats$se_std))
  params <- update_defined_param(params, "P_GPP_mediation",
                                 p_ratio_std,
                                 delta_p_ratio(c1_stats$std, d_stats$std, e1_stats$std,
                                               c1_stats$se_std, d_stats$se_std, e1_stats$se_std))
  params <- update_defined_param(params, "T_GPP_mediation",
                                 t_ratio_std,
                                 delta_p_ratio(c2_stats$std, d_stats$std, e2_stats$std,
                                               c2_stats$se_std, d_stats$se_std, e2_stats$se_std))
  params <- update_defined_param(params, "SW_GPP_mediation",
                                 sw_ratio_std,
                                 delta_p_ratio(c3_stats$std, d_stats$std, e3_stats$std,
                                               c3_stats$se_std, d_stats$se_std, e3_stats$se_std))

  fit_summary$pe <- update_defined_param(fit_summary$pe, "P_pre_via_EOS_GPP",
                                         a1_stats$std * b_stats$std * d_stats$std,
                                         delta_p_three(a1_stats$std, b_stats$std, d_stats$std,
                                                       a1_stats$se_std, b_stats$se_std, d_stats$se_std))
  fit_summary$pe <- update_defined_param(fit_summary$pe, "T_pre_via_EOS_GPP",
                                         a2_stats$std * b_stats$std * d_stats$std,
                                         delta_p_three(a2_stats$std, b_stats$std, d_stats$std,
                                                       a2_stats$se_std, b_stats$se_std, d_stats$se_std))
  fit_summary$pe <- update_defined_param(fit_summary$pe, "SW_pre_via_EOS_GPP",
                                         a3_stats$std * b_stats$std * d_stats$std,
                                         delta_p_three(a3_stats$std, b_stats$std, d_stats$std,
                                                       a3_stats$se_std, b_stats$se_std, d_stats$se_std))
  fit_summary$pe <- update_defined_param(fit_summary$pe, "P_pre_via_EOS",
                                         a1_stats$std * g_stats$std,
                                         delta_p_two(a1_stats$std, g_stats$std,
                                                     a1_stats$se_std, g_stats$se_std))
  fit_summary$pe <- update_defined_param(fit_summary$pe, "T_pre_via_EOS",
                                         a2_stats$std * g_stats$std,
                                         delta_p_two(a2_stats$std, g_stats$std,
                                                     a2_stats$se_std, g_stats$se_std))
  fit_summary$pe <- update_defined_param(fit_summary$pe, "SW_pre_via_EOS",
                                         a3_stats$std * g_stats$std,
                                         delta_p_two(a3_stats$std, g_stats$std,
                                                     a3_stats$se_std, g_stats$se_std))
  fit_summary$pe <- update_defined_param(fit_summary$pe, "P_pre_via_GPP",
                                         f1_stats$std * d_stats$std,
                                         delta_p_two(f1_stats$std, d_stats$std,
                                                     f1_stats$se_std, d_stats$se_std))
  fit_summary$pe <- update_defined_param(fit_summary$pe, "T_pre_via_GPP",
                                         f2_stats$std * d_stats$std,
                                         delta_p_two(f2_stats$std, d_stats$std,
                                                     f2_stats$se_std, d_stats$se_std))
  fit_summary$pe <- update_defined_param(fit_summary$pe, "SW_pre_via_GPP",
                                         f3_stats$std * d_stats$std,
                                         delta_p_two(f3_stats$std, d_stats$std,
                                                     f3_stats$se_std, d_stats$se_std))
  fit_summary$pe <- update_defined_param(fit_summary$pe, "P_season_via_EOS_GPP",
                                         a4_stats$std * b_stats$std * d_stats$std,
                                         delta_p_three(a4_stats$std, b_stats$std, d_stats$std,
                                                       a4_stats$se_std, b_stats$se_std, d_stats$se_std))
  fit_summary$pe <- update_defined_param(fit_summary$pe, "T_season_via_EOS_GPP",
                                         a5_stats$std * b_stats$std * d_stats$std,
                                         delta_p_three(a5_stats$std, b_stats$std, d_stats$std,
                                                       a5_stats$se_std, b_stats$se_std, d_stats$se_std))
  fit_summary$pe <- update_defined_param(fit_summary$pe, "SW_season_via_EOS_GPP",
                                         a6_stats$std * b_stats$std * d_stats$std,
                                         delta_p_three(a6_stats$std, b_stats$std, d_stats$std,
                                                       a6_stats$se_std, b_stats$se_std, d_stats$se_std))
  fit_summary$pe <- update_defined_param(fit_summary$pe, "P_season_via_EOS",
                                         a4_stats$std * g_stats$std,
                                         delta_p_two(a4_stats$std, g_stats$std,
                                                     a4_stats$se_std, g_stats$se_std))
  fit_summary$pe <- update_defined_param(fit_summary$pe, "T_season_via_EOS",
                                         a5_stats$std * g_stats$std,
                                         delta_p_two(a5_stats$std, g_stats$std,
                                                     a5_stats$se_std, g_stats$se_std))
  fit_summary$pe <- update_defined_param(fit_summary$pe, "SW_season_via_EOS",
                                         a6_stats$std * g_stats$std,
                                         delta_p_two(a6_stats$std, g_stats$std,
                                                     a6_stats$se_std, g_stats$se_std))
  fit_summary$pe <- update_defined_param(fit_summary$pe, "P_season_via_GPP",
                                         c1_stats$std * d_stats$std,
                                         delta_p_two(c1_stats$std, d_stats$std,
                                                     c1_stats$se_std, d_stats$se_std))
  fit_summary$pe <- update_defined_param(fit_summary$pe, "T_season_via_GPP",
                                         c2_stats$std * d_stats$std,
                                         delta_p_two(c2_stats$std, d_stats$std,
                                                     c2_stats$se_std, d_stats$se_std))
  fit_summary$pe <- update_defined_param(fit_summary$pe, "SW_season_via_GPP",
                                         c3_stats$std * d_stats$std,
                                         delta_p_two(c3_stats$std, d_stats$std,
                                                     c3_stats$se_std, d_stats$se_std))
  fit_summary$pe <- update_defined_param(fit_summary$pe, "EOS_via_GPP",
                                         b_stats$std * d_stats$std,
                                         delta_p_two(b_stats$std, d_stats$std,
                                                     b_stats$se_std, d_stats$se_std))
  fit_summary$pe <- update_defined_param(fit_summary$pe, "P_GPP_mediation",
                                         p_ratio_std,
                                         delta_p_ratio(c1_stats$std, d_stats$std, e1_stats$std,
                                                       c1_stats$se_std, d_stats$se_std, e1_stats$se_std))
  fit_summary$pe <- update_defined_param(fit_summary$pe, "T_GPP_mediation",
                                         t_ratio_std,
                                         delta_p_ratio(c2_stats$std, d_stats$std, e2_stats$std,
                                                       c2_stats$se_std, d_stats$se_std, e2_stats$se_std))
  fit_summary$pe <- update_defined_param(fit_summary$pe, "SW_GPP_mediation",
                                         sw_ratio_std,
                                         delta_p_ratio(c3_stats$std, d_stats$std, e3_stats$std,
                                                       c3_stats$se_std, d_stats$se_std, e3_stats$se_std))

  cat("\nSEM拟合摘要:\n")
  print(fit_summary)

  safe_write_csv(params, file.path(OUTPUT_DIR, "SEM_dual_timescale_parameters.csv"), row.names = FALSE)

  fit_measures <- fitMeasures(
    fit,
    c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr")
  )
  safe_write_csv(t(as.data.frame(fit_measures)), file.path(OUTPUT_DIR, "SEM_dual_timescale_fitindices.csv"))

  r2_vals <- inspect(fit, "r2")
  r2_filter <- filter_sem_r2(r2_vals)
  r2_vals <- r2_filter$r2
  if (FILTER_SEM_OUTLIERS) {
    cat(sprintf("  输出异常值过滤（全域SEM）：R²异常=%d\n",
                r2_filter$info["r2_invalid"]))
  }
  safe_write_csv(as.data.frame(r2_vals), file.path(OUTPUT_DIR, "SEM_dual_timescale_R2.csv"))

  global_filter_df <- data.frame(
    metric = c("coef_extreme", "p_invalid", "r2_invalid"),
    count = as.integer(c(param_filter$info["coef_extreme"],
                         param_filter$info["p_invalid"],
                         r2_filter$info["r2_invalid"]))
  )
  safe_write_csv(global_filter_df,
                 file.path(OUTPUT_DIR, "SEM_dual_timescale_outlier_filtering.csv"),
                 row.names = FALSE)

  summary_path <- file.path(OUTPUT_DIR, "SEM_dual_timescale_summary.txt")
  if (should_write(summary_path)) {
    sink(summary_path)
    print(fit_summary)
    if (FILTER_SEM_OUTLIERS) {
      cat("\n输出异常值过滤统计（全域SEM）:\n")
      print(global_filter_df)
    }
    sink()
  } else {
    cat(sprintf("  [skip] %s\n", summary_path))
  }

  # 路径图
  pdf_path <- file.path(OUTPUT_DIR, "SEM_dual_timescale_pathdiagram.pdf")
  plot_open <- should_write(pdf_path)
  if (plot_open) {
    pdf(pdf_path, width = 12, height = 10)
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
  } else {
    cat(sprintf("  [skip] %s\n", pdf_path))
  }

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

  sem_mode_orig <- SEM_SAMPLE_MODE
  detrend_orig <- DETREND_ENABLE

  # 运行双时间尺度SEM（先全域均值，再像元级）
  SEM_SAMPLE_MODE <<- "annual_mean"  # Bug Fix 5补充: 必须用<<-修改全局变量
  run_annual_mean_once <- function(enable_detrend, label) {
    SEM_SAMPLE_MODE <<- "annual_mean"  # Bug Fix 5: 设置annual_mean模式
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
    if (RUN_MODE == "skip" && annual_outputs_ready("raw")) {
      cat("  ✓ annual_mean(raw) 输出齐全，已跳过\n")
    } else {
      run_annual_mean_once(FALSE, "raw")
    }
    if (RUN_MODE == "skip" && annual_outputs_ready("detrended")) {
      cat("  ✓ annual_mean(detrended) 输出齐全，已跳过\n")
    } else {
      run_annual_mean_once(TRUE, "detrended")
    }
  } else {
    label <- if (DETREND_ENABLE) "detrended" else "raw"
    if (RUN_MODE == "skip" && annual_outputs_ready(label)) {
      cat(sprintf("  ✓ annual_mean(%s) 输出齐全，已跳过\n", label))
    } else {
      run_annual_mean_once(DETREND_ENABLE, label)
    }
  }

  # ===【修改10】像元时间序列模式：添加POSav气候态计算===
  set_output_dirs("")
  DETREND_ENABLE <- detrend_orig
  SEM_SAMPLE_MODE <<- "pixel_time_series"  # 必须用<<-修改全局变量
  detrend_pixel_orig <- DETREND_PIXEL_ENABLE

  pixel_runs <- list()
  if (RUN_BOTH_DETREND) {
    pixel_runs <- list(
      list(enable = FALSE, label = "raw", suffix = ""),
      list(enable = TRUE, label = "detrended", suffix = "_detrended")
    )
  } else {
    pixel_runs <- list(
      list(
        enable = DETREND_PIXEL_ENABLE,
        label = ifelse(DETREND_PIXEL_ENABLE, "detrended", "raw"),
        suffix = ifelse(DETREND_PIXEL_ENABLE, "_detrended", "")
      )
    )
  }

  pixel_all_skipped <- FALSE
  if (RUN_MODE == "skip") {
    pixel_all_skipped <- all(sapply(pixel_runs, function(pr) pixel_outputs_ready(pr$suffix)))
  }

  cat("\n=== 像元级SEM（pixel_time_series） ===\n")
  if (pixel_all_skipped) {
    cat("  ✓ 像元级输出齐全，已跳过\n")
  } else {
    # 先计算物候气候态（SOSav, POSav, EOSav）
    cat("\n=== 计算物候气候态（像元多年平均：SOSav, POSav, EOSav） ===\n")
    sos_clim_file <- file.path(DERIVED_DIR, "SOS_climatology.tif")
    pos_clim_file <- file.path(DERIVED_DIR, "POS_climatology.tif")
    eos_clim_file <- file.path(DERIVED_DIR, "EOS_climatology.tif")

    if (USE_EXISTING_CLIMATOLOGY &&
        file.exists(sos_clim_file) && file.exists(pos_clim_file) && file.exists(eos_clim_file)) {
      cat("  ✓ 复用已存在的物候气候态缓存\n")
      sos_climatology_r <- raster(sos_clim_file)
      pos_climatology_r <- raster(pos_clim_file)
      eos_climatology_r <- raster(eos_clim_file)
      sos_climatology_r <- set_nodata_if_missing(sos_climatology_r)
      pos_climatology_r <- set_nodata_if_missing(pos_climatology_r)
      eos_climatology_r <- set_nodata_if_missing(eos_climatology_r)
      sos_climatology_r <- mask_raster(sos_climatology_r, mask_r)
      pos_climatology_r <- mask_raster(pos_climatology_r, mask_r)
      eos_climatology_r <- mask_raster(eos_climatology_r, mask_r)
    } else {
      sos_stack_list <- list()
      pos_stack_list <- list()
      eos_stack_list <- list()
      n_loaded_sos <- 0
      n_loaded_pos <- 0
      n_loaded_eos <- 0

      for (year in years) {
        sos_file <- file.path(PHENO_DIR, sprintf("sos_gpp_%d.tif", year))
        pos_file <- file.path(PHENO_DIR, sprintf("pos_doy_gpp_%d.tif", year))
        eos_file <- file.path(PHENO_DIR, sprintf("eos_gpp_%d.tif", year))

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

        if (file.exists(eos_file)) {
          eos_r <- raster(eos_file)
          eos_r <- set_nodata_if_missing(eos_r)
          eos_r <- mask_raster(eos_r, mask_r)
          eos_stack_list[[length(eos_stack_list) + 1]] <- eos_r
          n_loaded_eos <- n_loaded_eos + 1
        }

        if ((n_loaded_sos %% 10 == 0) && (n_loaded_sos > 0)) {
          cat(sprintf("  已读取 %d/%d 年份的物候数据...\n", n_loaded_sos, length(years)))
        }
      }

      cat(sprintf("  ✓ 成功读取 %d/%d 年份的SOS数据\n", n_loaded_sos, length(years)))
      cat(sprintf("  ✓ 成功读取 %d/%d 年份的POS数据\n", n_loaded_pos, length(years)))
      cat(sprintf("  ✓ 成功读取 %d/%d 年份的EOS数据\n", n_loaded_eos, length(years)))

      if (n_loaded_sos == 0 || n_loaded_pos == 0 || n_loaded_eos == 0) {
        stop("未找到足够的SOS/POS/EOS数据文件")
      }

      cat("  正在计算SOS多年平均（SOSav）...\n")
      sos_stack <- stack(sos_stack_list)
      sos_climatology_r <- calc(sos_stack, fun = function(x) { mean(x, na.rm = TRUE) })
      sos_climatology_r <- mask_raster(sos_climatology_r, mask_r)

      cat("  正在计算POS多年平均（POSav）...\n")
      pos_stack <- stack(pos_stack_list)
      pos_climatology_r <- calc(pos_stack, fun = function(x) { mean(x, na.rm = TRUE) })
      pos_climatology_r <- mask_raster(pos_climatology_r, mask_r)

      cat("  正在计算EOS多年平均（EOSav）...\n")
      eos_stack <- stack(eos_stack_list)
      eos_climatology_r <- calc(eos_stack, fun = function(x) { mean(x, na.rm = TRUE) })
      eos_climatology_r <- mask_raster(eos_climatology_r, mask_r)
    }

    # ===【修改：添加POS > SOS和EOS > POS验证】===
    cat("  验证POSav > SOSav, EOSav > POSav...\n")
    valid_window_mask <- (pos_climatology_r > sos_climatology_r) & (eos_climatology_r > pos_climatology_r)
    valid_window_mask[is.na(valid_window_mask)] <- FALSE
    sos_climatology_r <- mask(sos_climatology_r, valid_window_mask, maskvalue = 0)
    pos_climatology_r <- mask(pos_climatology_r, valid_window_mask, maskvalue = 0)
    eos_climatology_r <- mask(eos_climatology_r, valid_window_mask, maskvalue = 0)

    n_invalid <- sum(getValues(valid_window_mask) == 0, na.rm = TRUE)
    cat(sprintf("  过滤掉的像元数: %d\n", n_invalid))
    if (n_invalid > 0) {
      cat(sprintf("  ⚠️ 过滤了 %d 个无效像元\n", n_invalid))
    }

    # 统计气候态信息
    sos_clim_vals <- getValues(sos_climatology_r)
    sos_clim_vals <- sanitize_values(sos_clim_vals, NAvalue(sos_climatology_r), FALSE)
    pos_clim_vals <- getValues(pos_climatology_r)
    pos_clim_vals <- sanitize_values(pos_clim_vals, NAvalue(pos_climatology_r), FALSE)
    eos_clim_vals <- getValues(eos_climatology_r)
    eos_clim_vals <- sanitize_values(eos_clim_vals, NAvalue(eos_climatology_r), FALSE)

    n_valid_sos <- sum(!is.na(sos_clim_vals))
    n_valid_pos <- sum(!is.na(pos_clim_vals))
    n_valid_eos <- sum(!is.na(eos_clim_vals))
    sos_min <- min(sos_clim_vals, na.rm = TRUE)
    sos_max <- max(sos_clim_vals, na.rm = TRUE)
    sos_mean <- mean(sos_clim_vals, na.rm = TRUE)
    pos_min <- min(pos_clim_vals, na.rm = TRUE)
    pos_max <- max(pos_clim_vals, na.rm = TRUE)
    pos_mean <- mean(pos_clim_vals, na.rm = TRUE)
    eos_min <- min(eos_clim_vals, na.rm = TRUE)
    eos_max <- max(eos_clim_vals, na.rm = TRUE)
    eos_mean <- mean(eos_clim_vals, na.rm = TRUE)

    cat(sprintf("  ✓ 物候气候态计算完成:\n"))
    cat(sprintf("    SOSav - 有效像元: %d, 范围: %.1f - %.1f DOY, 平均: %.1f DOY\n",
                n_valid_sos, sos_min, sos_max, sos_mean))
    cat(sprintf("    POSav - 有效像元: %d, 范围: %.1f - %.1f DOY, 平均: %.1f DOY\n",
                n_valid_pos, pos_min, pos_max, pos_mean))
    cat(sprintf("    EOSav - 有效像元: %d, 范围: %.1f - %.1f DOY, 平均: %.1f DOY\n",
                n_valid_eos, eos_min, eos_max, eos_mean))

    # 保存物候气候态（仅在未复用缓存时写入）
    if (!(USE_EXISTING_CLIMATOLOGY &&
          file.exists(sos_clim_file) && file.exists(pos_clim_file) && file.exists(eos_clim_file))) {
      safe_write_raster(sos_climatology_r, sos_clim_file, datatype = "FLT4S")
      safe_write_raster(pos_climatology_r, pos_clim_file, datatype = "FLT4S")
      safe_write_raster(eos_climatology_r, eos_clim_file, datatype = "FLT4S")
      cat(sprintf("  ✓ 物候气候态已保存: %s, %s, %s\n",
                  basename(sos_clim_file), basename(pos_clim_file), basename(eos_clim_file)))
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
    fixed_window_length_r[fixed_window_length_r <= 0] <- NA
    cat(sprintf("  ✓ Fixed_Window_Length已加载: %s\n", basename(fixed_len_file)))

    # 准备缓存
    prepare_sem_caches(years, sos_climatology_r, pos_climatology_r,
                      eos_climatology_r, fixed_window_length_r, mask_r)

    # 运行像元时间序列SEM（raw + detrended）
    for (pr in pixel_runs) {
      if (RUN_MODE == "skip" && pixel_outputs_ready(pr$suffix)) {
        cat(sprintf("  ✓ pixel_time_series(%s) 输出齐全，已跳过\n", pr$label))
        next
      }
      set_output_dirs(pr$suffix)
      DETREND_PIXEL_ENABLE <<- pr$enable  # Bug Fix 5补充: 必须用<<-修改全局变量
      cat(sprintf("\n=== 像元级SEM（pixel_time_series）[%s] ===\n", pr$label))
      run_sem_pixel_time_series(years, sos_climatology_r, fixed_window_length_r, mask_r)
    }
  }

  DETREND_PIXEL_ENABLE <<- detrend_pixel_orig
  SEM_SAMPLE_MODE <<- sem_mode_orig  # 必须用<<-修改全局变量

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
  cat("  Rscript 05b_SEM_analysis_dual_timescale_EOS.R\n")
  cat("======================================================================\n\n")
}
