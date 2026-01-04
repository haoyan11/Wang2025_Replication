#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Module 05: SEM analysis - Modified to use SOS instead of LSP
#
# 修改说明：
# - 原代码使用LSP（生长季长度 = POS - SOS + 1）
# - 本版本改用SOS（生长季起始日期），与04a代码保持一致
# - SOS更直接反映物候变化（提前/推迟）
#
# ⚠️ 重要发现：实际结果与Wang (2025)论文不同！
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 偏相关分析结果：
#   - SOS偏相关系数为负（r < 0）
#   - 含义：SOS提前（变小）→ TRproduct升高
#   - SEM预期：d5系数为负值
#
# Wang (2025)论文结果：
#   - TRproduct ~ +1.04 * ΔSOS
#   - 含义：SOS提前（ΔSOS < 0）→ TRproduct降低
#   - 论文预期：d5系数为正值
#
# 本代码采用实际数据结果，与论文相反！
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#
# SEM模型结构（基于Wang 2025 Figure 7）：
#   TRproduct (因变量)
#     ↑
#     ├── SOS (直接效应：更早的SOS → 更高的TRproduct)
#     ├── Ta (直接效应：温度)
#     ├── GPPspr (直接负效应 + 通过SMroot的间接效应)
#     ├── GPPsum (直接正效应 + 通过SMroot的间接效应)
#     └── SMroot (中介变量：土壤水分限制)
#
# Version: 2.1.0 (SOS-based, 基于实际数据结果)
# Author: Wang2025 Replication Project

suppressPackageStartupMessages({
  library(lavaan)
  library(semPlot)
  library(raster)
})

# ==================== 全局配置 ====================
# ⚠️ 重要：请根据您的数据路径修改以下配置！
# 提示：这些路径应与 _config.py 中的配置一致

# 根目录（请修改为您的实际路径）
# 自动检测操作系统并选择合适的路径格式
if (.Platform$OS.type == "windows") {
  # Windows环境（包括Windows上的RStudio）
  ROOT <- "I:/F/Data4"
} else {
  # Unix/Linux环境（WSL, Linux, macOS）
  if (dir.exists("/mnt/i/F/Data4")) {
    # WSL环境
    ROOT <- "/mnt/i/F/Data4"
  } else {
    # 其他Unix环境，尝试Windows路径
    ROOT <- "I:/F/Data4"
  }
}

# 自动生成的路径（通常不需要修改）
OUTPUT_ROOT <- file.path(ROOT, "Wang2025_Analysis")
PHENO_DIR <- file.path(ROOT, "Phenology_Output_1", "GPP_phenology_EPSG4326")
TRPRODUCT_DIR <- file.path(OUTPUT_ROOT, "Decomposition")
MASK_FILE <- file.path(OUTPUT_ROOT, "masks", "combined_mask.tif")

# 日尺度数据路径（与实际数据保持一致）
GPP_DAILY_DIR <- file.path(ROOT, "GLASS_GPP", "GLASS_GPP_daily_interpolated")
SMROOT_DAILY_DIR <- file.path(ROOT, "Meteorological Data", "GLEAM", "SMrz", "SMrz_Daily_1")  # 恢复原配置
TA_DAILY_DIR <- file.path(ROOT, "Meteorological Data", "ERA5_Land", "Tem", "Tem_Daily", "Tem_Daily_2")  # ✓ 与04a一致

# 输出目录
DATA_DIR <- file.path(OUTPUT_ROOT, "SEM_Data_SOS")
DERIVED_DIR <- file.path(DATA_DIR, "Derived")
OUTPUT_DIR <- file.path(OUTPUT_ROOT, "SEM_Results_SOS")
PIXELWISE_DIR <- file.path(OUTPUT_DIR, "Pixelwise")

# 创建输出目录
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DERIVED_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PIXELWISE_DIR, showWarnings = FALSE, recursive = TRUE)

# 年份范围（恢复完整范围以获得足够的SEM样本）
YEAR_START <- 1982
YEAR_END <- 2018

# 常量
NODATA_OUT <- -9999
NODATA_ABS_MAX <- 1e20

# ==================== 控制选项 ====================
# SEM分析模式：
#   "pixel_time_series" - 每个像元的时间序列单独分析（计算量大，推荐用于最终结果）
#   "annual_mean"       - 全域年度均值合并分析（快速，适合测试）
#   "pixel_year"        - 像元-年份对合并分析（中等，适合大规模）
SEM_SAMPLE_MODE <- "annual_mean"  # 建议测试时用"annual_mean"，正式运行用"pixel_time_series"

USE_CACHE <- TRUE                  # 启用缓存以加快处理速度（bug已修复）
MIN_VALID_FRAC <- 0.60            # LSP窗口内有效数据最低比例
MIN_VALID_YEARS <- 15             # 像元时间序列最少有效年份
WRITE_PIXELWISE_RASTERS <- FALSE  # 是否输出像元级栅格（会生成大量文件）

# 日尺度文件命名模式
GPP_DAILY_PATTERN <- "GPP_{date}.tif"
SMROOT_DAILY_PATTERN <- "SMrz_{date}.tif"
TA_DAILY_PATTERN <- "ERA5L_T2mDaily_C_{date}.tif"

# ==================== 路径验证 ====================
cat("\n======================================================================\n")
cat("05_SEM_analysis_SOS.R - 结构方程模型分析（基于SOS）\n")
cat("======================================================================\n")

cat("\n[环境检测]\n")
cat(sprintf("  操作系统类型: %s\n", .Platform$OS.type))
cat(sprintf("  R版本: %s\n", R.version.string))
cat(sprintf("  根目录: %s\n", ROOT))

cat("\n[路径检查]\n")
critical_paths <- list(
  "Root directory" = ROOT,
  "Phenology directory" = PHENO_DIR,
  "TRproduct directory" = TRPRODUCT_DIR,
  "Mask file" = MASK_FILE,
  "GPP daily directory" = GPP_DAILY_DIR,
  "SMroot daily directory" = SMROOT_DAILY_DIR,
  "Ta daily directory" = TA_DAILY_DIR
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

# ==================== 季节GPP计算（固定日历月份）====================
# 说明：
# - GPPspr: 3-5月GPP均值
# - GPPsum: 6-8月GPP均值
# - 使用固定日历月份，与论文方法一致
calc_seasonal_mean <- function(year, months, daily_dir, pattern, mask_r, cache_file, allow_negative = FALSE) {
  if (USE_CACHE && file.exists(cache_file)) {
    return(raster(cache_file))
  }

  dates <- date_seq_year(year)
  months_vec <- as.integer(strftime(dates, "%m"))
  keep_dates <- dates[months_vec %in% months]

  template <- mask_r
  sum_vals <- rep(0, ncell(template))
  cnt_vals <- rep(0L, ncell(template))

  files_found <- 0
  total_checked <- 0

  # 修复：使用索引循环而不是直接循环Date对象
  # 原因：for (d in keep_dates) 会将Date转换为numeric，导致日期错误（变成1970-01-01）
  for (i in seq_along(keep_dates)) {
    d <- keep_dates[i]  # 显式提取Date对象
    date_str <- strftime(d, "%Y%m%d")
    f <- build_daily_path(daily_dir, pattern, date_str)
    total_checked <- total_checked + 1
    if (!file.exists(f)) {
      next
    }
    files_found <- files_found + 1

    tryCatch({
      daily_r <- raster(f)
      daily_vals <- getValues(daily_r)
      daily_vals <- sanitize_values(daily_vals, NAvalue(daily_r), allow_negative)

      valid <- !is.na(daily_vals)
      sum_vals[valid] <- sum_vals[valid] + daily_vals[valid]
      cnt_vals[valid] <- cnt_vals[valid] + 1L
    }, error = function(e) {
      warning(sprintf("读取文件失败: %s", f))
    })
  }

  if (files_found == 0) {
    warning(sprintf("年份 %d 月份 %s 没有找到任何日尺度数据 (检查了 %d 个日期)",
                    year, paste(months, collapse=","), total_checked))
  } else {
    cat(sprintf("      找到 %d/%d 个文件\n", files_found, total_checked))
  }

  mean_vals <- sum_vals / cnt_vals
  mean_vals[cnt_vals == 0] <- NA

  out_r <- setValues(template, mean_vals)
  out_r <- mask_raster(out_r, mask_r)

  if (USE_CACHE) {
    tryCatch({
      writeRaster(out_r, cache_file, overwrite = TRUE, datatype = "FLT4S")
    }, error = function(e) {
      warning(sprintf("写入缓存文件失败: %s\n错误: %s", cache_file, conditionMessage(e)))
    })
  }
  out_r
}

# ==================== LSP期间气象变量均值（变化窗口）====================
# 说明：
# - 使用每年实际的[SOS, POS]窗口
# - 变量：Ta, SMroot（LSP期间均值）
# - 与论文方法一致
calc_lsp_mean <- function(year, sos_r, pos_r, daily_dir, pattern, mask_r, cache_file, allow_negative = TRUE) {
  if (USE_CACHE && file.exists(cache_file)) {
    return(raster(cache_file))
  }

  template <- mask_r
  sos_vals <- getValues(sos_r)
  pos_vals <- getValues(pos_r)

  sos_vals <- sanitize_values(sos_vals, NAvalue(sos_r), allow_negative = FALSE)
  pos_vals <- sanitize_values(pos_vals, NAvalue(pos_r), allow_negative = FALSE)

  sos_int_raw <- as.integer(round(sos_vals))
  pos_int_raw <- as.integer(round(pos_vals))

  # 处理NA值
  sos_int_raw[is.na(sos_int_raw)] <- 0L
  pos_int_raw[is.na(pos_int_raw)] <- 0L

  valid_pheno <- sos_int_raw > 0 & pos_int_raw > 0 &
    pos_int_raw >= sos_int_raw &
    sos_int_raw <= 366 & pos_int_raw <= 366

  sos_int <- pmin(pmax(sos_int_raw, 1), 365)
  pos_int <- pmin(pmax(pos_int_raw, 1), 365)

  sum_vals <- rep(0, ncell(template))
  cnt_vals <- rep(0L, ncell(template))

  if (!any(valid_pheno)) {
    out_r <- setValues(template, rep(NA_real_, ncell(template)))
    out_r <- mask_raster(out_r, mask_r)
    if (USE_CACHE) {
      writeRaster(out_r, cache_file, overwrite = TRUE)
    }
    return(out_r)
  }

  doy_start <- max(1, min(sos_int[valid_pheno], na.rm = TRUE))
  doy_end <- min(365, max(pos_int[valid_pheno], na.rm = TRUE))

  dates <- date_seq_year(year)
  files_found <- 0
  total_in_window <- 0

  # 调试：显示第一个文件路径（bug已修复，减少输出）
  debug_count <- 0
  debug_detail_count <- 0

  # 修复：使用索引循环而不是直接循环Date对象
  # 原因：for (d in dates) 会将Date转换为numeric，导致日期错误（变成1970-01-01）
  for (i in seq_along(dates)) {
    d <- dates[i]  # 显式提取Date对象
    doy <- as.integer(strftime(d, "%j"))
    if (doy < doy_start || doy > doy_end) {
      next
    }
    total_in_window <- total_in_window + 1

    date_str <- strftime(d, "%Y%m%d")
    f <- build_daily_path(daily_dir, pattern, date_str)

    # 调试：显示第1个文件（验证日期格式正确）
    if (debug_count < 1) {
      cat(sprintf("      [LSP调试] DOY %d: %s -> 存在=%s\n",
                  doy, basename(f), ifelse(file.exists(f), "是", "否")))
      debug_count <- debug_count + 1
    }

    if (!file.exists(f)) {
      next
    }
    files_found <- files_found + 1

    # 修复：将数据处理和累加分离，避免tryCatch作用域问题
    daily_result <- tryCatch({
      daily_r <- raster(f)
      daily_vals <- getValues(daily_r)
      daily_vals <- sanitize_values(daily_vals, NAvalue(daily_r), allow_negative)

      valid_daily <- !is.na(daily_vals)
      in_window <- valid_pheno & (sos_int <= doy) & (pos_int >= doy)
      use_mask <- in_window & valid_daily

      # 详细调试：显示第1次迭代的详细信息（验证累加逻辑）
      if (debug_detail_count < 1) {
        cat(sprintf("        [详细调试 DOY %d]\n", doy))
        cat(sprintf("          valid_daily 中TRUE数: %d / %d\n", sum(valid_daily), length(valid_daily)))
        cat(sprintf("          in_window 中TRUE数: %d / %d\n", sum(in_window), length(in_window)))
        cat(sprintf("          use_mask 中TRUE数: %d / %d\n", sum(use_mask), length(use_mask)))
        if (sum(valid_daily) > 0) {
          cat(sprintf("          数据值范围: [%.3f, %.3f]\n",
                      min(daily_vals[valid_daily], na.rm=TRUE),
                      max(daily_vals[valid_daily], na.rm=TRUE)))
        }
        debug_detail_count <<- debug_detail_count + 1
      }

      # 返回数据供外部累加
      list(daily_vals = daily_vals, use_mask = use_mask)
    }, error = function(e) {
      warning(sprintf("读取日尺度文件失败 (DOY=%d): %s", doy, f))
      NULL
    })

    # 在tryCatch外部进行累加，彻底避免作用域问题
    if (!is.null(daily_result) && any(daily_result$use_mask)) {
      sum_vals[daily_result$use_mask] <- sum_vals[daily_result$use_mask] + daily_result$daily_vals[daily_result$use_mask]
      cnt_vals[daily_result$use_mask] <- cnt_vals[daily_result$use_mask] + 1L

      # 调试：验证累加是否成功（仅显示第一次）
      if (debug_detail_count <= 1) {
        cat(sprintf("          ✓ 累加成功: %d 个像元，cnt_vals范围 [%d, %d]\n",
                    sum(daily_result$use_mask),
                    min(cnt_vals[daily_result$use_mask]),
                    max(cnt_vals[daily_result$use_mask])))
      }
    }
  }

  if (files_found == 0) {
    warning(sprintf("年份 %d LSP窗口内没有找到任何日尺度数据 (DOY %d-%d, 检查了 %d 个日期)",
                    year, doy_start, doy_end, total_in_window))
  } else {
    cat(sprintf("      找到 %d/%d 个文件 (DOY %d-%d)\n", files_found, total_in_window, doy_start, doy_end))
  }

  window_len <- pos_int - sos_int + 1
  good <- valid_pheno & (cnt_vals >= MIN_VALID_FRAC * window_len)

  # 调试：显示过滤统计
  cat(sprintf("      [过滤诊断]\n"))
  cat(sprintf("        有效物候像元: %d\n", sum(valid_pheno)))
  cat(sprintf("        通过数据量阈值(%.0f%%): %d\n", MIN_VALID_FRAC * 100, sum(good)))
  if (sum(valid_pheno) > 0) {
    cnt_sample <- cnt_vals[valid_pheno][1:min(5, sum(valid_pheno))]
    window_sample <- window_len[valid_pheno][1:min(5, sum(valid_pheno))]
    required_sample <- MIN_VALID_FRAC * window_sample
    cat(sprintf("        示例像元统计 (前5个):\n"))
    for (i in seq_along(cnt_sample)) {
      cat(sprintf("          像元%d: 有效天数=%d, 窗口长度=%d, 需要>=%d\n",
                  i, cnt_sample[i], window_sample[i], ceiling(required_sample[i])))
    }
  }

  mean_vals <- rep(NA_real_, ncell(template))
  mean_vals[good] <- sum_vals[good] / cnt_vals[good]

  out_r <- setValues(template, mean_vals)
  out_r <- mask_raster(out_r, mask_r)

  if (USE_CACHE) {
    tryCatch({
      writeRaster(out_r, cache_file, overwrite = TRUE, datatype = "FLT4S")
    }, error = function(e) {
      warning(sprintf("写入缓存文件失败: %s\n错误: %s", cache_file, conditionMessage(e)))
    })
  }
  out_r
}

# ==================== 数据准备 ====================
prepare_year_rasters <- function(year, mask_r) {
  trproduct_file <- file.path(TRPRODUCT_DIR, sprintf("TRproduct_%d.tif", year))
  sos_file <- file.path(PHENO_DIR, sprintf("sos_gpp_%d.tif", year))
  pos_file <- file.path(PHENO_DIR, sprintf("pos_doy_gpp_%d.tif", year))

  if (!file.exists(trproduct_file)) {
    cat(sprintf("    ⚠ 缺少TRproduct文件: %s\n", trproduct_file))
    return(NULL)
  }
  if (!file.exists(sos_file)) {
    cat(sprintf("    ⚠ 缺少SOS文件: %s\n", sos_file))
    return(NULL)
  }
  if (!file.exists(pos_file)) {
    cat(sprintf("    ⚠ 缺少POS文件: %s\n", pos_file))
    return(NULL)
  }

  trproduct_r <- tryCatch(raster(trproduct_file), error = function(e) {
    cat(sprintf("    ✗ 读取TRproduct失败: %s\n", conditionMessage(e)))
    return(NULL)
  })
  if (is.null(trproduct_r)) return(NULL)

  sos_r <- tryCatch(raster(sos_file), error = function(e) {
    cat(sprintf("    ✗ 读取SOS失败: %s\n", conditionMessage(e)))
    return(NULL)
  })
  if (is.null(sos_r)) return(NULL)

  pos_r <- tryCatch(raster(pos_file), error = function(e) {
    cat(sprintf("    ✗ 读取POS失败: %s\n", conditionMessage(e)))
    return(NULL)
  })
  if (is.null(pos_r)) return(NULL)

  NAvalue(trproduct_r) <- NODATA_OUT
  NAvalue(sos_r) <- NODATA_OUT
  NAvalue(pos_r) <- NODATA_OUT

  sos_r <- mask_raster(sos_r, mask_r)
  pos_r <- mask_raster(pos_r, mask_r)

  # 修改：直接使用SOS，不计算LSP
  # SOS: 生长季起始日期（DOY）
  sos_r <- mask_raster(sos_r, mask_r)

  gpp_spr_cache <- file.path(DERIVED_DIR, sprintf("GPPspr_%d.tif", year))
  gpp_sum_cache <- file.path(DERIVED_DIR, sprintf("GPPsum_%d.tif", year))
  smroot_cache <- file.path(DERIVED_DIR, sprintf("SMroot_LSP_%d.tif", year))
  ta_cache <- file.path(DERIVED_DIR, sprintf("Ta_LSP_%d.tif", year))

  cat(sprintf("    计算GPPspr (3-5月GPP均值)...\n"))
  gpp_spr_r <- calc_seasonal_mean(year, c(3, 4, 5), GPP_DAILY_DIR, GPP_DAILY_PATTERN,
                                 mask_r, gpp_spr_cache, allow_negative = FALSE)

  cat(sprintf("    计算GPPsum (6-8月GPP均值)...\n"))
  gpp_sum_r <- calc_seasonal_mean(year, c(6, 7, 8), GPP_DAILY_DIR, GPP_DAILY_PATTERN,
                                 mask_r, gpp_sum_cache, allow_negative = FALSE)

  cat(sprintf("    计算SMroot (LSP期间均值)...\n"))
  smroot_r <- calc_lsp_mean(year, sos_r, pos_r, SMROOT_DAILY_DIR, SMROOT_DAILY_PATTERN,
                           mask_r, smroot_cache, allow_negative = FALSE)

  cat(sprintf("    计算Ta (LSP期间均值)...\n"))
  ta_r <- calc_lsp_mean(year, sos_r, pos_r, TA_DAILY_DIR, TA_DAILY_PATTERN,
                       mask_r, ta_cache, allow_negative = TRUE)

  trproduct_r <- mask_raster(trproduct_r, mask_r)

  list(
    trproduct = trproduct_r,
    sos = sos_r,          # 修改：改为SOS
    gppspr = gpp_spr_r,
    gppsum = gpp_sum_r,
    smroot = smroot_r,
    ta = ta_r
  )
}

prepare_sem_caches <- function(years, mask_r) {
  cat("\n=== 准备SEM缓存数据 ===\n")
  for (year in years) {
    cat(sprintf("  年份: %d\n", year))
    tryCatch({
      prepare_year_rasters(year, mask_r)
      cat(sprintf("    ✓ 完成\n"))
    }, error = function(e) {
      cat(sprintf("    ✗ 错误: %s\n", conditionMessage(e)))
      cat(sprintf("    提示: 检查年份 %d 的输入数据是否完整\n", year))
    })
  }
}

prepare_sem_data <- function(years, mask_r) {
  cat("\n=== 准备SEM数据 ===\n")

  sem_data <- data.frame()

  for (year in years) {
    cat(sprintf("  年份: %d\n", year))

    rasters <- prepare_year_rasters(year, mask_r)
    if (is.null(rasters)) {
      cat(sprintf("    跳过 %d (缺少TRproduct或SOS/POS)\n", year))
      next
    }

    if (SEM_SAMPLE_MODE == "annual_mean") {
      year_row <- data.frame(
        year = year,
        TRproduct = cellStats(rasters$trproduct, mean, na.rm = TRUE),
        SOS = cellStats(rasters$sos, mean, na.rm = TRUE),  # 修改：改为SOS
        GPPspr = cellStats(rasters$gppspr, mean, na.rm = TRUE),
        GPPsum = cellStats(rasters$gppsum, mean, na.rm = TRUE),
        SMroot = cellStats(rasters$smroot, mean, na.rm = TRUE),
        Ta = cellStats(rasters$ta, mean, na.rm = TRUE)
      )
      sem_data <- rbind(sem_data, year_row)
    } else {
      vals <- data.frame(
        year = year,
        TRproduct = sanitize_values(getValues(rasters$trproduct), NAvalue(rasters$trproduct), allow_negative = TRUE),
        SOS = sanitize_values(getValues(rasters$sos), NAvalue(rasters$sos), allow_negative = FALSE),  # 修改：改为SOS
        GPPspr = sanitize_values(getValues(rasters$gppspr), NAvalue(rasters$gppspr), allow_negative = FALSE),
        GPPsum = sanitize_values(getValues(rasters$gppsum), NAvalue(rasters$gppsum), allow_negative = FALSE),
        SMroot = sanitize_values(getValues(rasters$smroot), NAvalue(rasters$smroot), allow_negative = FALSE),
        Ta = sanitize_values(getValues(rasters$ta), NAvalue(rasters$ta), allow_negative = TRUE)
      )
      vals <- vals[complete.cases(vals), ]
      sem_data <- rbind(sem_data, vals)
    }
  }

  # 数据质量检查
  cat("\n数据质量检查:\n")
  cat(sprintf("  总样本数: %d\n", nrow(sem_data)))
  cat(sprintf("  完整案例数: %d\n", sum(complete.cases(sem_data))))

  for (col in c("TRproduct", "SOS", "GPPspr", "GPPsum", "SMroot", "Ta")) {
    valid_n <- sum(is.finite(sem_data[[col]]))
    mean_val <- mean(sem_data[[col]], na.rm = TRUE)
    cat(sprintf("  %s: %d个有效值, 均值=%.3f\n", col, valid_n, mean_val))
  }

  # 检查是否有足够的有效数据
  if (sum(complete.cases(sem_data)) < 3) {
    cat("\n✗ 错误：完整案例数不足，无法进行SEM分析\n")
    cat("可能原因:\n")
    cat("  1. 日尺度数据文件路径不正确\n")
    cat("  2. 日尺度数据缺失\n")
    cat("  3. 掩膜过于严格导致无有效像元\n")
    cat("\n请检查以下路径:\n")
    cat(sprintf("  GPP日尺度: %s\n", GPP_DAILY_DIR))
    cat(sprintf("  SMroot日尺度: %s\n", SMROOT_DAILY_DIR))
    cat(sprintf("  Ta日尺度: %s\n", TA_DAILY_DIR))
    stop("数据准备失败：完整案例数不足")
  }

  # Z-score标准化
  sem_data_std <- sem_data
  sem_vars <- setdiff(names(sem_data_std), "year")
  sem_data_std[sem_vars] <- scale(sem_data_std[sem_vars])

  write.csv(sem_data, file.path(DATA_DIR, "sem_data_raw.csv"), row.names = FALSE)
  write.csv(sem_data_std, file.path(DATA_DIR, "sem_data_standardized.csv"), row.names = FALSE)

  cat("✓ SEM数据准备完成\n")
  sem_data_std
}

# ==================== VIF诊断 ====================
calculate_vif <- function(data) {
  # 修改：改为SOS
  lm_model <- lm(TRproduct ~ GPPspr + GPPsum + SMroot + Ta + SOS, data = data)
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

# ==================== SEM模型（全局合并）====================
run_sem_pooled <- function(data) {
  cat("\n=== SEM分析（全局合并）===\n")

  # 修改：将LSP改为SOS
  # ⚠️ 重要：SOS符号解释（基于实际偏相关结果）
  #   - 实际偏相关结果：SOS偏相关系数为负
  #   - 含义：SOS提前（变小）→ TRproduct升高
  #   - 因此：模型中 d5 预期为负值
  #   - 注意：这与Wang (2025)论文结果相反！论文是SOS提前→TRproduct降低
  model <- '
    # 携带效应（Carryover）
    GPPsum ~ a1*GPPspr

    # 土壤水分中介
    SMroot ~ b1*GPPspr + b2*GPPsum

    # 气候 -> 物候（温度影响SOS）
    # 注意：Ta升高 → SOS提前（变小），预期负相关
    SOS ~ c1*Ta

    # 蒸腾生产力组分
    # 修改：LSP改为SOS
    # ⚠️ 关键：基于偏相关结果，SOS提前（变小）→ TRproduct升高
    #    因此 d5 应为负值（SOS越小，TRproduct越大）
    TRproduct ~ d1*GPPspr + d2*GPPsum + d3*SMroot + d4*Ta + d5*SOS

    # 通过SMroot的间接效应
    indirect_GPPspr_SM_TRproduct := b1 * d3
    indirect_GPPsum_SM_TRproduct := b2 * d3
  '

  fit <- sem(model, data = data, estimator = "MLR")

  cat("\nSEM拟合摘要:\n")
  summary(fit, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

  params <- parameterEstimates(fit, standardized = TRUE)
  write.csv(params, file.path(OUTPUT_DIR, "SEM_parameters.csv"), row.names = FALSE)

  fit_measures <- fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))
  write.csv(t(as.data.frame(fit_measures)), file.path(OUTPUT_DIR, "SEM_fitindices.csv"))

  # 路径图
  pdf(file.path(OUTPUT_DIR, "SEM_pathdiagram.pdf"), width = 10, height = 8)
  semPaths(
    fit,
    what = "std",
    edge.label.cex = 1.2,
    curvePivot = TRUE,
    layout = "tree2",
    style = "lisrel",
    edge.color = "black",
    nodeLabels = c("TRproduct", "SOS", "Ta", "SMroot", "GPPspr", "GPPsum"),  # 修改：LSP改为SOS
    sizeMan = 10,
    residuals = FALSE,
    exoCov = FALSE
  )
  title("SEM (Wang 2025 - SOS版本)", line = 3)
  dev.off()

  cat("✓ SEM全局分析完成\n")
  cat(sprintf("  输出目录: %s\n", OUTPUT_DIR))
  fit
}

# ==================== SEM模型（像元时间序列）====================
run_sem_pixel_time_series <- function(years, mask_r) {
  cat("\n=== SEM分析（像元时间序列）===\n")
  cat("⚠️ 注意：此模式计算量极大，可能需要数小时！\n")
  cat("  建议：先用annual_mean模式测试\n\n")

  # 修改：将LSP改为SOS
  files <- list(
    TRproduct = file.path(TRPRODUCT_DIR, sprintf("TRproduct_%d.tif", years)),
    SOS = file.path(PHENO_DIR, sprintf("sos_gpp_%d.tif", years)),  # 修改：LSP改为SOS
    GPPspr = file.path(DERIVED_DIR, sprintf("GPPspr_%d.tif", years)),
    GPPsum = file.path(DERIVED_DIR, sprintf("GPPsum_%d.tif", years)),
    SMroot = file.path(DERIVED_DIR, sprintf("SMroot_LSP_%d.tif", years)),
    Ta = file.path(DERIVED_DIR, sprintf("Ta_LSP_%d.tif", years))
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

  template <- raster(files$TRproduct[1])
  bs <- blockSize(template)

  # 修改：路径名称
  coef_names <- c(
    "GPPsum~GPPspr",
    "SMroot~GPPspr",
    "SMroot~GPPsum",
    "SOS~Ta",                   # 修改：LSP改为SOS
    "TRproduct~GPPspr",
    "TRproduct~GPPsum",
    "TRproduct~SMroot",
    "TRproduct~Ta",
    "TRproduct~SOS",            # 修改：LSP改为SOS
    "indirect_GPPspr_SM_TRproduct",
    "indirect_GPPsum_SM_TRproduct"
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

    mats <- list(
      TRproduct = matrix(NA_real_, nrow = n_years, ncol = n_cells),
      SOS = matrix(NA_real_, nrow = n_years, ncol = n_cells),  # 修改：LSP改为SOS
      GPPspr = matrix(NA_real_, nrow = n_years, ncol = n_cells),
      GPPsum = matrix(NA_real_, nrow = n_years, ncol = n_cells),
      SMroot = matrix(NA_real_, nrow = n_years, ncol = n_cells),
      Ta = matrix(NA_real_, nrow = n_years, ncol = n_cells)
    )

    for (y in seq_len(n_years)) {
      mats$TRproduct[y, ] <- read_block_values(files$TRproduct[y], row, nrows, allow_negative = TRUE)
      mats$SOS[y, ] <- read_block_values(files$SOS[y], row, nrows, allow_negative = FALSE)  # 修改：LSP改为SOS
      mats$GPPspr[y, ] <- read_block_values(files$GPPspr[y], row, nrows, allow_negative = FALSE)
      mats$GPPsum[y, ] <- read_block_values(files$GPPsum[y], row, nrows, allow_negative = FALSE)
      mats$SMroot[y, ] <- read_block_values(files$SMroot[y], row, nrows, allow_negative = FALSE)
      mats$Ta[y, ] <- read_block_values(files$Ta[y], row, nrows, allow_negative = TRUE)
    }

    for (idx in valid_cells) {
      df <- data.frame(
        TRproduct = mats$TRproduct[, idx],
        SOS = mats$SOS[, idx],        # 修改：LSP改为SOS
        GPPspr = mats$GPPspr[, idx],
        GPPsum = mats$GPPsum[, idx],
        SMroot = mats$SMroot[, idx],
        Ta = mats$Ta[, idx]
      )

      df <- df[complete.cases(df), ]
      if (nrow(df) < MIN_VALID_YEARS) {
        next
      }

      # Z-score标准化
      trp_z <- scale_vec(df$TRproduct)
      sos_z <- scale_vec(df$SOS)        # 修改：LSP改为SOS
      gppspr_z <- scale_vec(df$GPPspr)
      gppsum_z <- scale_vec(df$GPPsum)
      sm_z <- scale_vec(df$SMroot)
      ta_z <- scale_vec(df$Ta)

      if (any(is.na(trp_z)) || any(is.na(sos_z)) || any(is.na(gppspr_z)) ||
          any(is.na(gppsum_z)) || any(is.na(sm_z)) || any(is.na(ta_z))) {
        next
      }

      # GPPsum ~ GPPspr
      a1_res <- corr_with_p(gppsum_z, gppspr_z)
      a1 <- a1_res$r
      p_a1 <- a1_res$p

      # SMroot ~ GPPspr + GPPsum
      X_sm <- cbind(gppspr_z, gppsum_z)
      b_res <- regress_beta_p(X_sm, sm_z)
      b <- b_res$beta
      p_b <- b_res$p

      # SOS ~ Ta（修改：LSP改为SOS）
      c1_res <- corr_with_p(sos_z, ta_z)
      c1 <- c1_res$r
      p_c1 <- c1_res$p

      # TRproduct ~ GPPspr + GPPsum + SMroot + Ta + SOS（修改：LSP改为SOS）
      vars_tr <- c("GPPspr", "GPPsum", "SMroot", "Ta", "SOS")  # 修改：LSP改为SOS
      X_tr_full <- cbind(gppspr_z, gppsum_z, sm_z, ta_z, sos_z)  # 修改：LSP改为SOS
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
        d_res <- regress_beta_p(X_tr, trp_z)
        if (!any(is.na(d_res$beta)) && !any(is.na(d_res$p))) {
          for (j in seq_along(vars_current)) {
            idx <- match(vars_current[j], vars_tr)
            d[idx] <- d_res$beta[j]
            p_d[idx] <- d_res$p[j]
          }
        }
      }

      if (any(is.na(c(a1, b, c1))) || any(is.na(c(p_a1, p_b, p_c1))) || all(is.na(d))) {
        next
      }

      vals <- c(
        a1,
        b[1], b[2],
        c1,
        d[1], d[2], d[3], d[4], d[5],
        b[1] * d[3],
        b[2] * d[3]
      )

      p_vals <- c(
        p_a1,
        p_b[1], p_b[2],
        p_c1,
        p_d[1], p_d[2], p_d[3], p_d[4], p_d[5],
        NA_real_,
        NA_real_
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

  write.csv(summary_df, file.path(PIXELWISE_DIR, "SEM_pixel_time_series_summary.csv"), row.names = FALSE)

  cat("\n像元时间序列SEM摘要:\n")
  print(summary_df)

  summary_df
}

# ==================== 主函数 ====================
main <- function() {
  cat("\n======================================================================\n")
  cat("SEM分析（基于SOS） - 与Wang (2025)对齐\n")
  cat("======================================================================\n")
  cat(sprintf("分析模式: %s\n", SEM_SAMPLE_MODE))
  cat(sprintf("年份范围: %d-%d\n", YEAR_START, YEAR_END))
  cat("----------------------------------------------------------------------\n")

  years <- YEAR_START:YEAR_END

  # 诊断日尺度数据文件（使用第一年作为测试）
  cat("\n[日尺度文件诊断 - 测试年份: 1982]\n")
  diagnose_daily_files(1982, GPP_DAILY_DIR, GPP_DAILY_PATTERN, sample_dates = 3)
  diagnose_daily_files(1982, SMROOT_DAILY_DIR, SMROOT_DAILY_PATTERN, sample_dates = 3)
  diagnose_daily_files(1982, TA_DAILY_DIR, TA_DAILY_PATTERN, sample_dates = 3)

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

  prepare_sem_caches(years, mask_r)

  if (SEM_SAMPLE_MODE == "pixel_time_series") {
    run_sem_pixel_time_series(years, mask_r)
  } else {
    data <- prepare_sem_data(years, mask_r)
    vif_vals <- calculate_vif(data)
    if (any(vif_vals > 10)) {
      cat("\n⚠️ 警告：检测到VIF > 10，建议移除高VIF变量\n")
    }
    run_sem_pooled(data)
  }

  cat("\n======================================================================\n")
  cat("✓ SEM分析完成\n")
  cat(sprintf("输出目录: %s\n", OUTPUT_DIR))
  cat("======================================================================\n")
}

# 自动运行逻辑
# - 非交互模式（Rscript）：自动运行
# - 交互模式（RStudio）：提示用户手动调用
if (!interactive()) {
  main()
} else {
  cat("\n======================================================================\n")
  cat("⚠️ 交互模式检测\n")
  cat("======================================================================\n")
  cat("请手动执行以下命令以运行SEM分析:\n\n")
  cat("  main()\n\n")
  cat("或者使用以下命令在非交互模式运行:\n\n")
  cat("  Rscript 05_SEM_analysis_SOS.R\n")
  cat("======================================================================\n\n")
}
