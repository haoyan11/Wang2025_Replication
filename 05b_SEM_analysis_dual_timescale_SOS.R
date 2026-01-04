#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Module 05b: 双时间尺度SEM分析 - 基于N04路径结构 + 04c固定窗口方法
#
# 设计思路：
# 1. 借鉴N04_Cal_SEM的双时间尺度路径结构
# 2. 使用04c代码的固定窗口方法（Fixed_Trate + 固定窗口[SOSav, POSav]）
# 3. 计算季前气候因子（SOSav前3个月）
# 4. 计算固定窗口生长季（SOSav-POSav）的GPP和气象因子
# 5. 构建完整路径：季前气候 → SOS → GPP_season → Fixed_Trate
# 6. 气候因子：降水(P)、气温(T)、短波辐射(SW)
#
# 路径结构（参考N04）：
#   季前气候(P,T,SW) → SOS → GPP_season → Fixed_Trate
#      ↓                ↓         ↑          ↑
#      └────────────────┴─────────┴──────────┘ (直接路径)
#   生长季气候(P,T,SW) ──────────┴──────────┘ (同期路径)
#
# ⚠️ 关键修改（对比原05b版本）：
# 1. TRproduct → Fixed_Trate（数据源：TR_fixed_window_{year}.tif，临时计算）
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

# 自动生成的路径
OUTPUT_ROOT <- file.path(ROOT, "Wang2025_Analysis")
PHENO_DIR <- file.path(ROOT, "Phenology_Output_1", "GPP_phenology_EPSG4326")
# ===【修改1】Fixed_Trate数据源（04c固定窗口方法）===
# 原路径: Decomposition/TRproduct_{year}.tif（03模块输出）
# 新路径: Decomposition/TR_fixed_window_{year}.tif（03c模块输出）
# 说明：Fixed_Trate需临时计算 = TR_fixed_window / Fixed_Window_Length
DECOMP_DIR <- file.path(OUTPUT_ROOT, "Decomposition")
MASK_FILE <- file.path(OUTPUT_ROOT, "masks", "combined_mask.tif")

# 日尺度数据路径
GPP_DAILY_DIR <- file.path(ROOT, "GLASS_GPP", "GLASS_GPP_daily_interpolated")
PRECIP_DAILY_DIR <- file.path(ROOT, "Meteorological Data", "ERA5_Land", "Pre", "Pre_Daily", "Pre_Daily_2")
TA_DAILY_DIR <- file.path(ROOT, "Meteorological Data", "ERA5_Land", "Tem", "Tem_Daily", "Tem_Daily_2")
SW_DAILY_DIR <- file.path(ROOT, "Meteorological Data", "ERA5_Land", "DSW", "DSW_Daily", "DSW_Daily_2")

# 输出目录
DATA_DIR <- file.path(OUTPUT_ROOT, "SEM_Dual_Timescale")
DERIVED_DIR <- file.path(DATA_DIR, "Derived")
OUTPUT_DIR <- file.path(OUTPUT_ROOT, "SEM_Dual_Timescale_Results")
PIXELWISE_DIR <- file.path(OUTPUT_DIR, "Pixelwise")

# 创建输出目录
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DERIVED_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PIXELWISE_DIR, showWarnings = FALSE, recursive = TRUE)

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
SEM_SAMPLE_MODE <- "annual_mean"  # 建议测试时用"annual_mean"，正式运行用"pixel_time_series"

USE_CACHE <- TRUE                  # 启用缓存以加快处理速度
MIN_VALID_FRAC <- 0.60            # SOS-POS窗口内有效数据最低比例
MIN_VALID_YEARS <- 15             # 像元时间序列最少有效年份
WRITE_PIXELWISE_RASTERS <- FALSE  # 是否输出像元级栅格（会生成大量文件）

# 日尺度文件命名模式
GPP_DAILY_PATTERN <- "GPP_{date}.tif"
PRECIP_DAILY_PATTERN <- "ERA5L_PrecipDaily_mm_{date}.tif"
TA_DAILY_PATTERN <- "ERA5L_T2mDaily_C_{date}.tif"
SW_DAILY_PATTERN <- "ERA5L_SWDaily_MJ_{date}.tif"

cat("\n======================================================================\n")
cat("双时间尺度SEM分析 - 基于N04路径结构 + 05数据源\n")
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
  "Decomposition directory" = DECOMP_DIR,  # 修改：TRproduct → Fixed_Trate（03c输出）
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

# 读取栅格块
read_block_values <- function(file_path, row, nrows, allow_negative = TRUE) {
  r <- raster(file_path)
  vals <- getValues(r, row = row, nrows = nrows)
  nodata <- NAvalue(r)
  sanitize_values(vals, nodata, allow_negative)
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

  # 逐像元计算（按窗口分组以提高效率）
  # 1. 将像元按SOS气候态分组
  sos_int <- round(sos_clim_valid)
  sos_int <- pmax(90, pmin(365, sos_int))  # 限制到合理范围

  unique_sos <- sort(unique(sos_int))
  n_groups <- length(unique_sos)

  cat(sprintf("      不同的SOS窗口数: %d\n", n_groups))

  # 2. 对每个SOS窗口，批量处理对应的像元
  processed_cells <- 0

  for (sos_doy in unique_sos) {
    # 找到该SOS值对应的像元
    cells_this_sos <- valid_cells[sos_int == sos_doy]
    n_cells <- length(cells_this_sos)

    # 确定季前窗口
    preseason_start_doy <- sos_doy - 90
    preseason_end_doy <- sos_doy - 1

    # 处理跨年情况
    if (preseason_start_doy < 1) {
      doys_prev_year <- (365 + preseason_start_doy):365
      doys_this_year <- 1:preseason_end_doy
      doys <- c(doys_prev_year, doys_this_year)
    } else {
      doys <- preseason_start_doy:preseason_end_doy
    }

    # 初始化这组像元的累加数组
    sum_vals_group <- rep(0, n_cells)
    cnt_vals_group <- rep(0L, n_cells)

    # 读取并累加该窗口内的日尺度数据
    for (doy in doys) {
      # 处理跨年
      if (doy > 300 && preseason_start_doy < 1) {
        year_to_use <- year - 1
      } else {
        year_to_use <- year
      }

      # 将DOY转换为日期
      date_obj <- as.Date(doy - 1, origin = paste0(year_to_use, "-01-01"))
      date_str <- strftime(date_obj, "%Y%m%d")

      file_path <- build_daily_path(daily_dir, pattern, date_str)

      if (!file.exists(file_path)) next

      tryCatch({
        daily_r <- raster(file_path)
        daily_vals_all <- getValues(daily_r)
        daily_vals_all <- sanitize_values(daily_vals_all, NAvalue(daily_r), allow_negative)

        # 只提取这组像元的值
        daily_vals_group <- daily_vals_all[cells_this_sos]

        valid_mask <- !is.na(daily_vals_group)

        sum_vals_group[valid_mask] <- sum_vals_group[valid_mask] + daily_vals_group[valid_mask]
        cnt_vals_group[valid_mask] <- cnt_vals_group[valid_mask] + 1L
      }, error = function(e) {
        # 静默跳过
      })
    }

    # 计算这组像元的均值
    good <- cnt_vals_group > 0
    if (any(good)) {
      mean_vals[cells_this_sos[good]] <- sum_vals_group[good] / cnt_vals_group[good]
    }

    processed_cells <- processed_cells + n_cells

    # 每处理20%的像元，打印一次进度
    if (processed_cells %% max(1, floor(n_valid / 5)) == 0 ||
        sos_doy == unique_sos[length(unique_sos)]) {
      cat(sprintf("      进度: %d/%d 像元 (%.1f%%), SOS=%d DOY, 窗口=%d-%d\n",
                  processed_cells, n_valid, 100 * processed_cells / n_valid,
                  sos_doy, preseason_start_doy, preseason_end_doy))
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

  valid_pheno <- !is.na(sos_int) & !is.na(pos_int) & (sos_int <= pos_int)

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

  # 生成日期序列
  year_start <- as.Date(paste0(year, "-01-01"))
  year_dates <- seq(year_start, by = "day", length.out = 365)

  # 只处理物候窗口内的日期
  window_dates <- year_dates[min_doy:max_doy]

  # 逐日累加
  for (i in seq_along(window_dates)) {
    date_obj <- window_dates[i]  # 显式提取Date对象
    doy <- as.integer(format(date_obj, "%j"))
    date_str <- strftime(date_obj, "%Y%m%d")

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

  cat(sprintf("      有效像元: %d (cnt > 0)\n", sum(cnt_vals > 0)))

  # 计算均值
  mean_vals <- rep(NA_real_, ncell(template))
  good <- cnt_vals > 0
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
  sos_r <- raster(sos_file)

  NAvalue(tr_fixed_r) <- NODATA_OUT
  NAvalue(sos_r) <- NODATA_OUT

  # 应用掩膜
  sos_r <- mask_raster(sos_r, mask_r)
  tr_fixed_r <- mask_raster(tr_fixed_r, mask_r)

  # 计算Fixed_Trate = TR_fixed_window / Fixed_Window_Length
  cat("  [计算Fixed_Trate] TR_fixed_window / Fixed_Window_Length...\n")
  fixed_trate_r <- tr_fixed_r / fixed_window_length_r
  fixed_trate_r <- mask_raster(fixed_trate_r, mask_r)

  # 缓存文件路径
  cache_prefix <- file.path(DATA_DIR, sprintf("%s_%d.tif", c(
    "P_pre", "T_pre", "SW_pre",
    "P_season", "T_season", "SW_season",
    "GPP_season"
  ), year))

  names(cache_prefix) <- c("P_pre", "T_pre", "SW_pre",
                           "P_season", "T_season", "SW_season",
                           "GPP_season")

  # ===== 计算季前气候因子（使用SOS气候态） =====
  cat("  [1/7] 季前气候因子（基于像元多年平均SOS）:\n")

  p_pre_r <- calc_preseason_climate(year, sos_climatology_r, PRECIP_DAILY_DIR, PRECIP_DAILY_PATTERN,
                                     mask_r, cache_prefix["P_pre"], "P_pre", allow_negative = FALSE)

  t_pre_r <- calc_preseason_climate(year, sos_climatology_r, TA_DAILY_DIR, TA_DAILY_PATTERN,
                                     mask_r, cache_prefix["T_pre"], "T_pre", allow_negative = TRUE)

  sw_pre_r <- calc_preseason_climate(year, sos_climatology_r, SW_DAILY_DIR, SW_DAILY_PATTERN,
                                      mask_r, cache_prefix["SW_pre"], "SW_pre", allow_negative = FALSE)

  # ===【修改3b】计算固定窗口生长季气候因子 =====
  cat("  [2/7] 生长季气候因子（固定窗口[SOSav, POSav]）:\n")

  p_season_r <- calc_season_climate_fixed(year, sos_climatology_r, pos_climatology_r, PRECIP_DAILY_DIR, PRECIP_DAILY_PATTERN,
                                           mask_r, cache_prefix["P_season"], "P_season", allow_negative = FALSE)

  t_season_r <- calc_season_climate_fixed(year, sos_climatology_r, pos_climatology_r, TA_DAILY_DIR, TA_DAILY_PATTERN,
                                           mask_r, cache_prefix["T_season"], "T_season", allow_negative = TRUE)

  sw_season_r <- calc_season_climate_fixed(year, sos_climatology_r, pos_climatology_r, SW_DAILY_DIR, SW_DAILY_PATTERN,
                                            mask_r, cache_prefix["SW_season"], "SW_season", allow_negative = FALSE)

  # ===【修改3c】计算固定窗口生长季GPP =====
  cat("  [3/7] 生长季GPP（固定窗口[SOSav, POSav]）:\n")
  gpp_season_r <- calc_season_climate_fixed(year, sos_climatology_r, pos_climatology_r, GPP_DAILY_DIR, GPP_DAILY_PATTERN,
                                             mask_r, cache_prefix["GPP_season"], "GPP_season", allow_negative = FALSE)

  cat("  ✓ 数据准备完成\n")

  # ===【修改3d】返回Fixed_Trate（替代TRc）===
  list(
    fixed_trate = fixed_trate_r,  # 修改：trc → fixed_trate
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
  lm_model <- lm(Fixed_Trate ~ SOS + GPP_season + P_pre + T_pre + SW_pre +
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
    TR_fixed_window = file.path(DECOMP_DIR, sprintf("TR_fixed_window_%d.tif", years)),  # 修改：TRc → TR_fixed_window
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

  if (length(years) < MIN_VALID_YEARS) {
    stop("像元时间序列SEM所需的完整输入年份不足")
  }

  template <- raster(files$TR_fixed_window[1])  # 修改：TRc → TR_fixed_window
  bs <- blockSize(template)

  # 提取Fixed_Window_Length值（只需读取一次）
  fixed_window_length_vals <- getValues(fixed_window_length_r)

  # 路径名称（双时间尺度）
  coef_names <- c(
    "SOS~P_pre",
    "SOS~T_pre",
    "SOS~SW_pre",
    "GPP_season~SOS",
    "GPP_season~P_pre",
    "GPP_season~T_pre",
    "GPP_season~SW_pre",
    "GPP_season~P_season",
    "GPP_season~T_season",
    "GPP_season~SW_season",
    "Fixed_Trate~SOS",  # 修改：TRc → Fixed_Trate
    "Fixed_Trate~GPP_season",
    "Fixed_Trate~P_pre",
    "Fixed_Trate~T_pre",
    "Fixed_Trate~SW_pre",
    "Fixed_Trate~P_season",
    "Fixed_Trate~T_season",
    "Fixed_Trate~SW_season",
    "P_pre_via_SOS_GPP",
    "T_pre_via_SOS_GPP",
    "SW_pre_via_SOS_GPP"
  )

  coef_sum <- setNames(rep(0, length(coef_names)), coef_names)
  coef_sumsq <- setNames(rep(0, length(coef_names)), coef_names)
  coef_count <- setNames(rep(0, length(coef_names)), coef_names)
  sig_sum <- setNames(rep(0, length(coef_names)), coef_names)
  sig_sumsq <- setNames(rep(0, length(coef_names)), coef_names)
  sig_count <- setNames(rep(0, length(coef_names)), coef_names)

  for (i in seq_len(bs$n)) {
    row <- bs$row[i]
    nrows <- bs$nrows[i]

    mask_vals <- getValues(mask_r, row = row, nrows = nrows)
    valid_cells <- which(!is.na(mask_vals))
    if (length(valid_cells) == 0) {
      next
    }

    n_cells <- length(mask_vals)
    n_years <- length(years)

    # 获取当前块的Fixed_Window_Length值
    fixed_len_block <- fixed_window_length_vals[((row-1)*ncol(template)+1):((row+nrows-1)*ncol(template))]
    if (length(fixed_len_block) != n_cells) {
      # 处理最后一个不完整的块
      fixed_len_block <- fixed_window_length_vals[((row-1)*ncol(template)+1):min(length(fixed_window_length_vals), (row+nrows-1)*ncol(template))]
    }

    mats <- list(
      Fixed_Trate = matrix(NA_real_, nrow = n_years, ncol = n_cells),  # 修改：TRc → Fixed_Trate
      SOS = matrix(NA_real_, nrow = n_years, ncol = n_cells),
      GPP_season = matrix(NA_real_, nrow = n_years, ncol = n_cells),
      P_pre = matrix(NA_real_, nrow = n_years, ncol = n_cells),
      T_pre = matrix(NA_real_, nrow = n_years, ncol = n_cells),
      SW_pre = matrix(NA_real_, nrow = n_years, ncol = n_cells),
      P_season = matrix(NA_real_, nrow = n_years, ncol = n_cells),
      T_season = matrix(NA_real_, nrow = n_years, ncol = n_cells),
      SW_season = matrix(NA_real_, nrow = n_years, ncol = n_cells)
    )

    for (y in seq_len(n_years)) {
      # 读取TR_fixed_window并计算Fixed_Trate
      tr_fixed_vals <- read_block_values(files$TR_fixed_window[y], row, nrows, allow_negative = TRUE)
      # Fixed_Trate = TR_fixed_window / Fixed_Window_Length
      mats$Fixed_Trate[y, ] <- tr_fixed_vals / fixed_len_block

      mats$SOS[y, ] <- read_block_values(files$SOS[y], row, nrows, allow_negative = FALSE)
      mats$GPP_season[y, ] <- read_block_values(files$GPP_season[y], row, nrows, allow_negative = FALSE)
      mats$P_pre[y, ] <- read_block_values(files$P_pre[y], row, nrows, allow_negative = FALSE)
      mats$T_pre[y, ] <- read_block_values(files$T_pre[y], row, nrows, allow_negative = TRUE)
      mats$SW_pre[y, ] <- read_block_values(files$SW_pre[y], row, nrows, allow_negative = FALSE)
      mats$P_season[y, ] <- read_block_values(files$P_season[y], row, nrows, allow_negative = FALSE)
      mats$T_season[y, ] <- read_block_values(files$T_season[y], row, nrows, allow_negative = TRUE)
      mats$SW_season[y, ] <- read_block_values(files$SW_season[y], row, nrows, allow_negative = FALSE)
    }

    for (idx in valid_cells) {
      df <- data.frame(
        Fixed_Trate = mats$Fixed_Trate[, idx],  # 修改：TRc → Fixed_Trate
        SOS = mats$SOS[, idx],
        GPP_season = mats$GPP_season[, idx],
        P_pre = mats$P_pre[, idx],
        T_pre = mats$T_pre[, idx],
        SW_pre = mats$SW_pre[, idx],
        P_season = mats$P_season[, idx],
        T_season = mats$T_season[, idx],
        SW_season = mats$SW_season[, idx]
      )

      df <- df[complete.cases(df), ]
      if (nrow(df) < MIN_VALID_YEARS) {
        next
      }

      # Z-score标准化
      fixed_trate_z <- scale_vec(df$Fixed_Trate)  # 修改：trc_z → fixed_trate_z
      sos_z <- scale_vec(df$SOS)
      gpp_z <- scale_vec(df$GPP_season)
      p_pre_z <- scale_vec(df$P_pre)
      t_pre_z <- scale_vec(df$T_pre)
      sw_pre_z <- scale_vec(df$SW_pre)
      p_season_z <- scale_vec(df$P_season)
      t_season_z <- scale_vec(df$T_season)
      sw_season_z <- scale_vec(df$SW_season)

      if (any(is.na(c(fixed_trate_z, sos_z, gpp_z, p_pre_z, t_pre_z, sw_pre_z,  # 修改：trc_z → fixed_trate_z
                      p_season_z, t_season_z, sw_season_z)))) {
        next
      }

      # SOS ~ P_pre + T_pre + SW_pre
      X_sos <- cbind(p_pre_z, t_pre_z, sw_pre_z)
      a_res <- regress_beta_p(X_sos, sos_z)
      a <- a_res$beta
      p_a <- a_res$p

      # GPP_season ~ SOS + P_pre + T_pre + SW_pre + P_season + T_season + SW_season
      X_gpp <- cbind(sos_z, p_pre_z, t_pre_z, sw_pre_z, p_season_z, t_season_z, sw_season_z)
      bc_res <- regress_beta_p(X_gpp, gpp_z)
      bc <- bc_res$beta
      p_bc <- bc_res$p
      b <- bc[1]
      c <- bc[2:4]
      f <- bc[5:7]
      p_b <- p_bc[1]
      p_c <- p_bc[2:4]
      p_f <- p_bc[5:7]

      # Fixed_Trate ~ SOS + GPP + P_pre + T_pre + SW_pre + P_season + T_season + SW_season
      vars_tr <- c("SOS", "GPP_season", "P_pre", "T_pre", "SW_pre",
                   "P_season", "T_season", "SW_season")
      X_tr_full <- cbind(sos_z, gpp_z, p_pre_z, t_pre_z, sw_pre_z,
                         p_season_z, t_season_z, sw_season_z)
      X_tr <- X_tr_full
      vars_current <- vars_tr

      # VIF过滤
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
        d_res <- regress_beta_p(X_tr, fixed_trate_z)  # 修改：trc_z → fixed_trate_z
        if (!any(is.na(d_res$beta)) && !any(is.na(d_res$p))) {
          for (j in seq_along(vars_current)) {
            idx_match <- match(vars_current[j], vars_tr)
            d[idx_match] <- d_res$beta[j]
            p_d[idx_match] <- d_res$p[j]
          }
        }
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

      for (k in seq_along(coef_names)) {
        if (is.finite(vals[k])) {
          coef_sum[coef_names[k]] <- coef_sum[coef_names[k]] + vals[k]
          coef_sumsq[coef_names[k]] <- coef_sumsq[coef_names[k]] + vals[k]^2
          coef_count[coef_names[k]] <- coef_count[coef_names[k]] + 1
          if (is.finite(p_vals[k]) && p_vals[k] < 0.05) {
            sig_sum[coef_names[k]] <- sig_sum[coef_names[k]] + vals[k]
            sig_sumsq[coef_names[k]] <- sig_sumsq[coef_names[k]] + vals[k]^2
            sig_count[coef_names[k]] <- sig_count[coef_names[k]] + 1
          }
        }
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
    pos_file <- file.path(PHENO_DIR, sprintf("pos_gpp_%d.tif", year))

    if (file.exists(sos_file)) {
      sos_r <- raster(sos_file)
      NAvalue(sos_r) <- NODATA_OUT
      sos_r <- mask_raster(sos_r, mask_r)
      sos_stack_list[[length(sos_stack_list) + 1]] <- sos_r
      n_loaded_sos <- n_loaded_sos + 1
    }

    if (file.exists(pos_file)) {
      pos_r <- raster(pos_file)
      NAvalue(pos_r) <- NODATA_OUT
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

  # 保存物候气候态
  sos_clim_file <- file.path(DATA_DIR, "SOS_climatology.tif")
  pos_clim_file <- file.path(DATA_DIR, "POS_climatology.tif")
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
  NAvalue(fixed_window_length_r) <- NODATA_OUT
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
        Fixed_Trate = cellStats(rasters$fixed_trate, mean, na.rm = TRUE),  # 修改：TRc → Fixed_Trate, trc → fixed_trate
        SOS = cellStats(rasters$sos, mean, na.rm = TRUE),
        GPP_season = cellStats(rasters$gpp_season, mean, na.rm = TRUE),
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
        Fixed_Trate = sanitize_values(getValues(rasters$fixed_trate), NAvalue(rasters$fixed_trate), TRUE),  # 修改：TRc → Fixed_Trate, trc → fixed_trate
        SOS = sanitize_values(getValues(rasters$sos), NAvalue(rasters$sos), FALSE),
        GPP_season = sanitize_values(getValues(rasters$gpp_season), NAvalue(rasters$gpp_season), FALSE),
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

  # 修改：TRc → Fixed_Trate
  for (col in c("Fixed_Trate", "SOS", "GPP_season", "P_pre", "T_pre", "SW_pre",
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

  # Z-score标准化
  sem_data_std <- sem_data
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
  cat("✓ 因变量: Fixed_Trate（固定窗口蒸腾速率）\n")

  # 完整路径模型（降水 + 温度 + 短波辐射）
  sem_model <- '
    # 第一层：季前气候 → SOS（物候响应）
    SOS ~ a1*P_pre + a2*T_pre + a3*SW_pre

    # 第二层：SOS + 季前气候 + 生长季气候 → GPP_season（碳固定）
    GPP_season ~ b*SOS + f1*P_pre + f2*T_pre + f3*SW_pre +
                 c1*P_season + c2*T_season + c3*SW_season

    # 第三层：SOS + GPP_season + 季前气候 + 生长季气候 → Fixed_Trate（固定窗口蒸腾速率）
    Fixed_Trate ~ g*SOS + d*GPP_season +
          h1*P_pre + h2*T_pre + h3*SW_pre +
          e1*P_season + e2*T_season + e3*SW_season

    # === 间接效应分解 ===

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

  cat("\nSEM模型结构:\n")
  cat(sem_model)
  cat("\n")

  # ===== 拟合SEM模型 =====
  cat("\n=== 拟合SEM模型 ===\n")

  fit <- sem(sem_model, data = sem_data_std, estimator = "MLR")

  cat("\nSEM拟合摘要:\n")
  summary(fit, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

  # 保存结果
  params <- parameterEstimates(fit, standardized = TRUE)
  write.csv(params, file.path(OUTPUT_DIR, "SEM_dual_timescale_parameters.csv"), row.names = FALSE)

  fit_measures <- fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))
  write.csv(t(as.data.frame(fit_measures)), file.path(OUTPUT_DIR, "SEM_dual_timescale_fitindices.csv"))

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

  # 测试读取一个实际的日尺度文件
  cat("\n[测试读取日尺度文件]\n")
  test_gpp_file <- file.path(GPP_DAILY_DIR, "GPP_19820301.tif")
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

  cat("----------------------------------------------------------------------\n")

  # 运行双时间尺度SEM
  # ===【修改10】像元时间序列模式：添加POSav气候态计算===
  if (SEM_SAMPLE_MODE == "pixel_time_series") {
    # 先计算物候气候态（SOSav, POSav）
    cat("\n=== 计算物候气候态（像元多年平均：SOSav, POSav） ===\n")
    sos_stack_list <- list()
    pos_stack_list <- list()
    n_loaded_sos <- 0
    n_loaded_pos <- 0

    for (year in years) {
      sos_file <- file.path(PHENO_DIR, sprintf("sos_gpp_%d.tif", year))
      pos_file <- file.path(PHENO_DIR, sprintf("pos_gpp_%d.tif", year))

      if (file.exists(sos_file)) {
        sos_r <- raster(sos_file)
        NAvalue(sos_r) <- NODATA_OUT
        sos_r <- mask_raster(sos_r, mask_r)
        sos_stack_list[[length(sos_stack_list) + 1]] <- sos_r
        n_loaded_sos <- n_loaded_sos + 1
      }

      if (file.exists(pos_file)) {
        pos_r <- raster(pos_file)
        NAvalue(pos_r) <- NODATA_OUT
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

    # 保存物候气候态
    sos_clim_file <- file.path(DATA_DIR, "SOS_climatology.tif")
    pos_clim_file <- file.path(DATA_DIR, "POS_climatology.tif")
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
    NAvalue(fixed_window_length_r) <- NODATA_OUT
    fixed_window_length_r <- mask_raster(fixed_window_length_r, mask_r)
    cat(sprintf("  ✓ Fixed_Window_Length已加载: %s\n", basename(fixed_len_file)))

    # 准备缓存
    prepare_sem_caches(years, sos_climatology_r, pos_climatology_r,
                      fixed_window_length_r, mask_r)

    # 运行像元时间序列SEM
    run_sem_pixel_time_series(years, sos_climatology_r, fixed_window_length_r, mask_r)
  } else {
    # 全局合并模式
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
  }

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
  cat("  Rscript 05b_SEM_analysis_dual_timescale_SOS.R\n")
  cat("======================================================================\n\n")
}
