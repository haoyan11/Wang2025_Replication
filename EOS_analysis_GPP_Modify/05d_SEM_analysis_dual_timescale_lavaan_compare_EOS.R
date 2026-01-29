#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Pixel-time-series SEM (lavaan) with dual screening comparison
#
# Purpose:
# - Keep the 05b model structure (EOS -> GPP -> Fixed_Trate).
# - Use per-pixel year series (39-41 years) and lavaan::sem().
# - Keep 05b missing-data handling (valid-year fraction, no interpolation).
# - Compare two screening/significance strategies:
#     (A) "ours": 05b-style outlier filtering + p-value significance summary
#     (B) "other": GFI-based pixel filtering + p-value masking (non-sig -> NA)
#
# Output (aligned with 05b/05c pixelwise format):
# - Pixelwise summary files (ours vs other) with 05b-style columns
# - Pixelwise R2 summary files (ours vs other)
# - Comparison table (mean + sig_frac deltas)
#
# This script is standalone and does NOT modify 05b.

suppressPackageStartupMessages({
  library(raster)
  library(lavaan)
  library(parallel)
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

cat(sprintf("\n[并行化配置] 使用 %d 个CPU核心进行并行计算\n", PARALLEL_CORES))

# === Run mode ===
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

# ==================== Paths ====================
if (.Platform$OS.type == "windows") {
  ROOT <- "I:/F/Data4"
} else {
  if (dir.exists("/mnt/i/F/Data4")) {
    ROOT <- "/mnt/i/F/Data4"
  } else {
    ROOT <- "I:/F/Data4"
  }
}

OUTPUT_ROOT <- file.path(ROOT, "Wang2025_Analysis_EOS_GPP_Modify")
PHENO_DIR <- file.path(ROOT, "Phenology_Output_1", "GPP_phenology")
DECOMP_DIR <- file.path(OUTPUT_ROOT, "Decomposition_FixedWindow")
MASK_FILE <- file.path(OUTPUT_ROOT, "masks", "combined_mask.tif")
TEMPLATE_FILE <- file.path(OUTPUT_ROOT, "masks", "template_grid.tif")
DERIVED_DIR <- file.path(OUTPUT_ROOT, "SEM_Data_Dual_Fixed", "Derived")

OUTPUT_DIR_BASE <- file.path(OUTPUT_ROOT, "SEM_Results_Dual_Fixed_Lavaan_Compare")
OUTPUT_DIR <- OUTPUT_DIR_BASE
PIXELWISE_DIR <- file.path(OUTPUT_DIR, "Pixelwise")
COMPARE_DIR <- file.path(OUTPUT_DIR_BASE, "Compare")

set_output_dirs <- function(suffix = "", method_tag = "Ours") {
  compare_dir <- file.path(OUTPUT_DIR_BASE, paste0("Compare", suffix))
  output_dir <- paste0(OUTPUT_DIR_BASE, "_", method_tag, suffix)
  pixelwise_dir <- file.path(output_dir, "Pixelwise")

  OUTPUT_DIR <<- output_dir
  PIXELWISE_DIR <<- pixelwise_dir
  COMPARE_DIR <<- compare_dir

  dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
  dir.create(PIXELWISE_DIR, showWarnings = FALSE, recursive = TRUE)
  dir.create(COMPARE_DIR, showWarnings = FALSE, recursive = TRUE)

  list(output_dir = output_dir, pixelwise_dir = pixelwise_dir, compare_dir = compare_dir)
}

outputs_ready <- function(suffix = "") {
  ours_dir <- paste0(OUTPUT_DIR_BASE, "_Ours", suffix)
  other_dir <- paste0(OUTPUT_DIR_BASE, "_Other", suffix)
  compare_dir <- file.path(OUTPUT_DIR_BASE, paste0("Compare", suffix))

  required_pixelwise <- c(
    "SEM_dual_timescale_parameters.csv",
    "SEM_dual_timescale_pixelwise_bootstrap_ci.csv",
    "SEM_dual_timescale_R2_detail.csv",
    "SEM_dual_timescale_R2.csv",
    "SEM_dual_timescale_pixelwise_outlier_filtering.csv",
    "SEM_dual_timescale_summary.txt"
  )

  required_compare <- c(
    "SEM_dual_timescale_pixel_time_series_compare.csv",
    "SEM_dual_timescale_pixel_time_series_filtering.csv"
  )

  ours_ok <- all(file.exists(file.path(ours_dir, "Pixelwise", required_pixelwise)))
  other_ok <- all(file.exists(file.path(other_dir, "Pixelwise", required_pixelwise)))
  compare_ok <- all(file.exists(file.path(compare_dir, required_compare)))

  ours_ok && other_ok && compare_ok
}

# ==================== Options ====================
YEAR_START <- 1982
YEAR_END <- 2018

MIN_VALID_YEAR_FRAC <- 0.60  # "Ours"方法的最小年份比例
OTHER_MIN_VALID_YEAR_FRAC <- 1.00  # "Other"方法要求100%数据完整（对应N04）
DETREND_PIXEL_ENABLE <- TRUE
RUN_BOTH_DETREND <- TRUE
MEDIATION_DENOM_EPS <- 1e-6
PIXELWISE_BOOTSTRAP_ENABLE <- TRUE
PIXELWISE_BOOTSTRAP_N <- 800
PIXELWISE_BOOTSTRAP_SEED <- 202501


# 05b-style outlier filtering
FILTER_SEM_OUTLIERS <- TRUE
SEM_COEF_ABS_MAX <- 5
SEM_P_MIN <- 0
SEM_P_MAX <- 1
SEM_R2_MIN <- 0
SEM_R2_MAX <- 1

# "Other" screening
OTHER_GFI_MIN <- 0.90
OTHER_P_THRESHOLD <- 0.05

# ==================== Helpers ====================
check_file <- function(file_path, description) {
  if (!file.exists(file_path)) {
    stop(sprintf("✗ 文件不存在: %s\n  路径: %s", description, file_path))
  }
  cat(sprintf("  [OK] %s\n    %s\n", description, file_path))
}

check_dir <- function(dir_path, description) {
  if (!dir.exists(dir_path)) {
    stop(sprintf("✗ 目录不存在: %s\n  路径: %s", description, dir_path))
  }
  cat(sprintf("  [OK] %s\n    %s\n", description, dir_path))
}

mask_raster <- function(r, mask_r) {
  mask(r, mask_r)
}

sanitize_values <- function(vals, nodata, allow_negative = TRUE) {
  vals[vals == nodata] <- NA_real_
  vals[abs(vals) > 1e20] <- NA_real_  # 过滤极端异常值
  if (!allow_negative) {
    vals[vals < 0] <- NA_real_
  }
  vals[is.infinite(vals)] <- NA_real_
  vals
}

safe_mediation_ratio <- function(c_val, d_val, e_val, eps = MEDIATION_DENOM_EPS) {
  denom <- e_val + c_val * d_val
  if (!is.finite(denom) || abs(denom) < eps) {
    return(NA_real_)
  }
  (c_val * d_val) / denom
}

delta_var_two <- function(a, b, se_a, se_b) {
  if (!is.finite(a) || !is.finite(b) || !is.finite(se_a) || !is.finite(se_b)) {
    return(NA_real_)
  }
  (b^2) * (se_a^2) + (a^2) * (se_b^2)
}

delta_var_three <- function(a, b, c, se_a, se_b, se_c) {
  if (!is.finite(a) || !is.finite(b) || !is.finite(c) ||
      !is.finite(se_a) || !is.finite(se_b) || !is.finite(se_c)) {
    return(NA_real_)
  }
  (b^2 * c^2) * (se_a^2) +
    (a^2 * c^2) * (se_b^2) +
    (a^2 * b^2) * (se_c^2)
}

delta_p_two <- function(a, b, se_a, se_b) {
  var_ab <- delta_var_two(a, b, se_a, se_b)
  if (!is.finite(var_ab) || var_ab <= 0) {
    return(NA_real_)
  }
  z <- (a * b) / sqrt(var_ab)
  2 * (1 - stats::pnorm(abs(z)))
}

delta_p_three <- function(a, b, c, se_a, se_b, se_c) {
  var_abc <- delta_var_three(a, b, c, se_a, se_b, se_c)
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

detrend_series <- function(x, t) {
  if (length(x) < 3 || all(is.na(x))) return(x)
  ok <- is.finite(x) & is.finite(t)
  if (sum(ok) < 3) return(x)
  fit <- lm(x[ok] ~ t[ok])
  x[ok] <- x[ok] - fitted(fit)
  x
}

scale_vec <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  mu <- mean(x, na.rm = TRUE)
  sdv <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) return(rep(NA_real_, length(x)))
  (x - mu) / sdv
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

# ==================== SEM model (same as 05b) ====================
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

  # === 间接效应定义（lavaan自动计算系数和p值）===
  P_pre_via_EOS_GPP  := a1 * b * d
  T_pre_via_EOS_GPP  := a2 * b * d
  SW_pre_via_EOS_GPP := a3 * b * d
  P_pre_via_EOS  := a1 * g
  T_pre_via_EOS  := a2 * g
  SW_pre_via_EOS := a3 * g
  P_pre_via_GPP  := f1 * d
  T_pre_via_GPP  := f2 * d
  SW_pre_via_GPP := f3 * d
  P_pre_indirect  := a1*b*d + a1*g + f1*d
  T_pre_indirect  := a2*b*d + a2*g + f2*d
  SW_pre_indirect := a3*b*d + a3*g + f3*d
  P_season_via_EOS_GPP  := a4 * b * d
  T_season_via_EOS_GPP  := a5 * b * d
  SW_season_via_EOS_GPP := a6 * b * d
  P_season_via_EOS  := a4 * g
  T_season_via_EOS  := a5 * g
  SW_season_via_EOS := a6 * g
  P_season_via_GPP  := c1 * d
  T_season_via_GPP  := c2 * d
  SW_season_via_GPP := c3 * d
  P_season_indirect  := a4*b*d + a4*g + c1*d
  T_season_indirect  := a5*b*d + a5*g + c2*d
  SW_season_indirect := a6*b*d + a6*g + c3*d
  EOS_via_GPP := b * d

  # GPP中介比例
  P_GPP_mediation := (c1*d) / (e1 + c1*d)
  T_GPP_mediation := (c2*d) / (e2 + c2*d)
  SW_GPP_mediation := (c3*d) / (e3 + c3*d)
'

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
  "P_pre_indirect",
  "T_pre_indirect",
  "SW_pre_indirect",
  "P_season_via_EOS_GPP",
  "T_season_via_EOS_GPP",
  "SW_season_via_EOS_GPP",
  "P_season_via_EOS",
  "T_season_via_EOS",
  "SW_season_via_EOS",
  "P_season_via_GPP",
  "T_season_via_GPP",
  "SW_season_via_GPP",
  "P_season_indirect",
  "T_season_indirect",
  "SW_season_indirect",
  "EOS_via_GPP",
  "P_GPP_mediation",
  "T_GPP_mediation",
  "SW_GPP_mediation"
)

r2_names <- c("EOS", "Fixed_GPPrate", "Fixed_Trate")

extract_lavaan_results <- function(df_std) {
  fit <- tryCatch(
    lavaan::sem(sem_model, data = df_std, estimator = "MLR"),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(NULL)
  }

  params <- lavaan::parameterEstimates(fit, standardized = TRUE)

  get_label_std <- function(label) {
    idx <- which(params$label == label)
    if (length(idx) == 0) return(NA_real_)
    params$std.all[idx][1]
  }

  get_label_p <- function(label) {
    idx <- which(params$label == label)
    if (length(idx) == 0) return(NA_real_)
    params$pvalue[idx][1]
  }

  get_label_stats <- function(label) {
    idx <- which(params$label == label)
    if (length(idx) == 0) {
      return(list(std = NA_real_, se_std = NA_real_))
    }
    est <- params$est[idx][1]
    std <- params$std.all[idx][1]
    se <- params$se[idx][1]
    se_std <- if (is.finite(est) && est != 0 && is.finite(std) && is.finite(se)) {
      abs(std / est) * se
    } else {
      NA_real_
    }
    list(std = std, se_std = se_std)
  }

  a1 <- get_label_std("a1")
  a2 <- get_label_std("a2")
  a3 <- get_label_std("a3")
  a4 <- get_label_std("a4")
  a5 <- get_label_std("a5")
  a6 <- get_label_std("a6")
  b <- get_label_std("b")
  f1 <- get_label_std("f1")
  f2 <- get_label_std("f2")
  f3 <- get_label_std("f3")
  c1 <- get_label_std("c1")
  c2 <- get_label_std("c2")
  c3 <- get_label_std("c3")
  g <- get_label_std("g")
  d <- get_label_std("d")
  h1 <- get_label_std("h1")
  h2 <- get_label_std("h2")
  h3 <- get_label_std("h3")
  e1 <- get_label_std("e1")
  e2 <- get_label_std("e2")
  e3 <- get_label_std("e3")

  p_a1 <- get_label_p("a1")
  p_a2 <- get_label_p("a2")
  p_a3 <- get_label_p("a3")
  p_a4 <- get_label_p("a4")
  p_a5 <- get_label_p("a5")
  p_a6 <- get_label_p("a6")
  p_b <- get_label_p("b")
  p_f1 <- get_label_p("f1")
  p_f2 <- get_label_p("f2")
  p_f3 <- get_label_p("f3")
  p_c1 <- get_label_p("c1")
  p_c2 <- get_label_p("c2")
  p_c3 <- get_label_p("c3")
  p_g <- get_label_p("g")
  p_d <- get_label_p("d")
  p_h1 <- get_label_p("h1")
  p_h2 <- get_label_p("h2")
  p_h3 <- get_label_p("h3")
  p_e1 <- get_label_p("e1")
  p_e2 <- get_label_p("e2")
  p_e3 <- get_label_p("e3")

  # 间接效应（标准化系数乘积，保持与直接效应尺度一致）
  p_via_eos_gpp <- a1 * b * d
  t_via_eos_gpp <- a2 * b * d
  sw_via_eos_gpp <- a3 * b * d
  p_via_eos <- a1 * g
  t_via_eos <- a2 * g
  sw_via_eos <- a3 * g
  p_via_gpp <- f1 * d
  t_via_gpp <- f2 * d
  sw_via_gpp <- f3 * d
  p_pre_indirect <- a1 * b * d + a1 * g + f1 * d
  t_pre_indirect <- a2 * b * d + a2 * g + f2 * d
  sw_pre_indirect <- a3 * b * d + a3 * g + f3 * d
  p_season_via_eos_gpp <- a4 * b * d
  t_season_via_eos_gpp <- a5 * b * d
  sw_season_via_eos_gpp <- a6 * b * d
  p_season_via_eos <- a4 * g
  t_season_via_eos <- a5 * g
  sw_season_via_eos <- a6 * g
  p_season_via_gpp <- c1 * d
  t_season_via_gpp <- c2 * d
  sw_season_via_gpp <- c3 * d
  p_season_indirect <- a4 * b * d + a4 * g + c1 * d
  t_season_indirect <- a5 * b * d + a5 * g + c2 * d
  sw_season_indirect <- a6 * b * d + a6 * g + c3 * d
  eos_via_gpp <- b * d
  p_ratio <- safe_mediation_ratio(c1, d, e1)
  t_ratio <- safe_mediation_ratio(c2, d, e2)
  sw_ratio <- safe_mediation_ratio(c3, d, e3)

  a1_stats <- get_label_stats("a1")
  a2_stats <- get_label_stats("a2")
  a3_stats <- get_label_stats("a3")
  a4_stats <- get_label_stats("a4")
  a5_stats <- get_label_stats("a5")
  a6_stats <- get_label_stats("a6")
  b_stats <- get_label_stats("b")
  f1_stats <- get_label_stats("f1")
  f2_stats <- get_label_stats("f2")
  f3_stats <- get_label_stats("f3")
  c1_stats <- get_label_stats("c1")
  c2_stats <- get_label_stats("c2")
  c3_stats <- get_label_stats("c3")
  g_stats <- get_label_stats("g")
  d_stats <- get_label_stats("d")
  e1_stats <- get_label_stats("e1")
  e2_stats <- get_label_stats("e2")
  e3_stats <- get_label_stats("e3")

  p_p_via_eos_gpp <- delta_p_three(a1_stats$std, b_stats$std, d_stats$std,
                                   a1_stats$se_std, b_stats$se_std, d_stats$se_std)
  p_t_via_eos_gpp <- delta_p_three(a2_stats$std, b_stats$std, d_stats$std,
                                   a2_stats$se_std, b_stats$se_std, d_stats$se_std)
  p_sw_via_eos_gpp <- delta_p_three(a3_stats$std, b_stats$std, d_stats$std,
                                    a3_stats$se_std, b_stats$se_std, d_stats$se_std)
  p_p_via_eos <- delta_p_two(a1_stats$std, g_stats$std, a1_stats$se_std, g_stats$se_std)
  p_t_via_eos <- delta_p_two(a2_stats$std, g_stats$std, a2_stats$se_std, g_stats$se_std)
  p_sw_via_eos <- delta_p_two(a3_stats$std, g_stats$std, a3_stats$se_std, g_stats$se_std)
  p_p_via_gpp <- delta_p_two(f1_stats$std, d_stats$std, f1_stats$se_std, d_stats$se_std)
  p_t_via_gpp <- delta_p_two(f2_stats$std, d_stats$std, f2_stats$se_std, d_stats$se_std)
  p_sw_via_gpp <- delta_p_two(f3_stats$std, d_stats$std, f3_stats$se_std, d_stats$se_std)
  p_p_season_via_eos_gpp <- delta_p_three(a4_stats$std, b_stats$std, d_stats$std,
                                          a4_stats$se_std, b_stats$se_std, d_stats$se_std)
  p_t_season_via_eos_gpp <- delta_p_three(a5_stats$std, b_stats$std, d_stats$std,
                                          a5_stats$se_std, b_stats$se_std, d_stats$se_std)
  p_sw_season_via_eos_gpp <- delta_p_three(a6_stats$std, b_stats$std, d_stats$std,
                                           a6_stats$se_std, b_stats$se_std, d_stats$se_std)
  p_p_season_via_eos <- delta_p_two(a4_stats$std, g_stats$std, a4_stats$se_std, g_stats$se_std)
  p_t_season_via_eos <- delta_p_two(a5_stats$std, g_stats$std, a5_stats$se_std, g_stats$se_std)
  p_sw_season_via_eos <- delta_p_two(a6_stats$std, g_stats$std, a6_stats$se_std, g_stats$se_std)
  p_p_season_via_gpp <- delta_p_two(c1_stats$std, d_stats$std, c1_stats$se_std, d_stats$se_std)
  p_t_season_via_gpp <- delta_p_two(c2_stats$std, d_stats$std, c2_stats$se_std, d_stats$se_std)
  p_sw_season_via_gpp <- delta_p_two(c3_stats$std, d_stats$std, c3_stats$se_std, d_stats$se_std)
  p_eos_via_gpp <- delta_p_two(b_stats$std, d_stats$std, b_stats$se_std, d_stats$se_std)
  p_p_ratio <- delta_p_ratio(c1_stats$std, d_stats$std, e1_stats$std,
                             c1_stats$se_std, d_stats$se_std, e1_stats$se_std)
  p_t_ratio <- delta_p_ratio(c2_stats$std, d_stats$std, e2_stats$std,
                             c2_stats$se_std, d_stats$se_std, e2_stats$se_std)
  p_sw_ratio <- delta_p_ratio(c3_stats$std, d_stats$std, e3_stats$std,
                              c3_stats$se_std, d_stats$se_std, e3_stats$se_std)

  p_pre_indirect_p <- {
    var_sum <- sum(c(
      delta_var_three(a1_stats$std, b_stats$std, d_stats$std, a1_stats$se_std, b_stats$se_std, d_stats$se_std),
      delta_var_two(a1_stats$std, g_stats$std, a1_stats$se_std, g_stats$se_std),
      delta_var_two(f1_stats$std, d_stats$std, f1_stats$se_std, d_stats$se_std)
    ), na.rm = TRUE)
    if (is.finite(var_sum) && var_sum > 0) {
      2 * (1 - stats::pnorm(abs(p_pre_indirect / sqrt(var_sum))))
    } else {
      NA_real_
    }
  }

  t_pre_indirect_p <- {
    var_sum <- sum(c(
      delta_var_three(a2_stats$std, b_stats$std, d_stats$std, a2_stats$se_std, b_stats$se_std, d_stats$se_std),
      delta_var_two(a2_stats$std, g_stats$std, a2_stats$se_std, g_stats$se_std),
      delta_var_two(f2_stats$std, d_stats$std, f2_stats$se_std, d_stats$se_std)
    ), na.rm = TRUE)
    if (is.finite(var_sum) && var_sum > 0) {
      2 * (1 - stats::pnorm(abs(t_pre_indirect / sqrt(var_sum))))
    } else {
      NA_real_
    }
  }

  sw_pre_indirect_p <- {
    var_sum <- sum(c(
      delta_var_three(a3_stats$std, b_stats$std, d_stats$std, a3_stats$se_std, b_stats$se_std, d_stats$se_std),
      delta_var_two(a3_stats$std, g_stats$std, a3_stats$se_std, g_stats$se_std),
      delta_var_two(f3_stats$std, d_stats$std, f3_stats$se_std, d_stats$se_std)
    ), na.rm = TRUE)
    if (is.finite(var_sum) && var_sum > 0) {
      2 * (1 - stats::pnorm(abs(sw_pre_indirect / sqrt(var_sum))))
    } else {
      NA_real_
    }
  }

  p_season_indirect_p <- {
    var_sum <- sum(c(
      delta_var_three(a4_stats$std, b_stats$std, d_stats$std, a4_stats$se_std, b_stats$se_std, d_stats$se_std),
      delta_var_two(a4_stats$std, g_stats$std, a4_stats$se_std, g_stats$se_std),
      delta_var_two(c1_stats$std, d_stats$std, c1_stats$se_std, d_stats$se_std)
    ), na.rm = TRUE)
    if (is.finite(var_sum) && var_sum > 0) {
      2 * (1 - stats::pnorm(abs(p_season_indirect / sqrt(var_sum))))
    } else {
      NA_real_
    }
  }

  t_season_indirect_p <- {
    var_sum <- sum(c(
      delta_var_three(a5_stats$std, b_stats$std, d_stats$std, a5_stats$se_std, b_stats$se_std, d_stats$se_std),
      delta_var_two(a5_stats$std, g_stats$std, a5_stats$se_std, g_stats$se_std),
      delta_var_two(c2_stats$std, d_stats$std, c2_stats$se_std, d_stats$se_std)
    ), na.rm = TRUE)
    if (is.finite(var_sum) && var_sum > 0) {
      2 * (1 - stats::pnorm(abs(t_season_indirect / sqrt(var_sum))))
    } else {
      NA_real_
    }
  }

  sw_season_indirect_p <- {
    var_sum <- sum(c(
      delta_var_three(a6_stats$std, b_stats$std, d_stats$std, a6_stats$se_std, b_stats$se_std, d_stats$se_std),
      delta_var_two(a6_stats$std, g_stats$std, a6_stats$se_std, g_stats$se_std),
      delta_var_two(c3_stats$std, d_stats$std, c3_stats$se_std, d_stats$se_std)
    ), na.rm = TRUE)
    if (is.finite(var_sum) && var_sum > 0) {
      2 * (1 - stats::pnorm(abs(sw_season_indirect / sqrt(var_sum))))
    } else {
      NA_real_
    }
  }

  vals <- c(
    a1, a2, a3, a4, a5, a6,
    b, f1, f2, f3, c1, c2, c3,
    g, d, h1, h2, h3, e1, e2, e3,
    p_via_eos_gpp, t_via_eos_gpp, sw_via_eos_gpp,
    p_via_eos, t_via_eos, sw_via_eos,
    p_via_gpp, t_via_gpp, sw_via_gpp,
    p_pre_indirect, t_pre_indirect, sw_pre_indirect,
    p_season_via_eos_gpp, t_season_via_eos_gpp, sw_season_via_eos_gpp,
    p_season_via_eos, t_season_via_eos, sw_season_via_eos,
    p_season_via_gpp, t_season_via_gpp, sw_season_via_gpp,
    p_season_indirect, t_season_indirect, sw_season_indirect,
    eos_via_gpp,
    p_ratio, t_ratio, sw_ratio
  )

  p_vals <- c(
    p_a1, p_a2, p_a3, p_a4, p_a5, p_a6,
    p_b, p_f1, p_f2, p_f3, p_c1, p_c2, p_c3,
    p_g, p_d, p_h1, p_h2, p_h3, p_e1, p_e2, p_e3,
    p_p_via_eos_gpp, p_t_via_eos_gpp, p_sw_via_eos_gpp,
    p_p_via_eos, p_t_via_eos, p_sw_via_eos,
    p_p_via_gpp, p_t_via_gpp, p_sw_via_gpp,
    p_pre_indirect_p, t_pre_indirect_p, sw_pre_indirect_p,
    p_p_season_via_eos_gpp, p_t_season_via_eos_gpp, p_sw_season_via_eos_gpp,
    p_p_season_via_eos, p_t_season_via_eos, p_sw_season_via_eos,
    p_p_season_via_gpp, p_t_season_via_gpp, p_sw_season_via_gpp,
    p_season_indirect_p, t_season_indirect_p, sw_season_indirect_p,
    p_eos_via_gpp,
    p_p_ratio, p_t_ratio, p_sw_ratio
  )

  r2_vals <- rep(NA_real_, length(r2_names))
  r2_fit <- tryCatch(lavaan::inspect(fit, "r2"), error = function(e) NULL)
  if (!is.null(r2_fit)) {
    for (nm in r2_names) {
      if (nm %in% names(r2_fit)) {
        r2_vals[match(nm, r2_names)] <- r2_fit[[nm]]
      }
    }
  }

  fit_measures <- tryCatch(
    lavaan::fitMeasures(fit, c("gfi", "rmsea", "cfi", "chisq", "df", "pvalue")),
    error = function(e) NULL
  )

  list(vals = vals, p_vals = p_vals, r2_vals = r2_vals, fit_measures = fit_measures)
}

build_r2_summary <- function(r2_sum, r2_sumsq, r2_count) {
  r2_mean <- r2_sum / r2_count
  r2_sd <- sqrt(pmax(0, r2_sumsq / r2_count - r2_mean^2))
  r2_mean[!is.finite(r2_mean)] <- NA_real_
  r2_sd[!is.finite(r2_sd)] <- NA_real_
  data.frame(
    variable = r2_names,
    R2 = r2_mean,
    R2_sd = r2_sd,
    n = r2_count
  )
}

write_pixelwise_outputs <- function(pixelwise_dir, summary_df, r2_summary, filter_df, label) {
  dir.create(pixelwise_dir, showWarnings = FALSE, recursive = TRUE)

  safe_write_csv(summary_df,
                 file.path(pixelwise_dir, "SEM_dual_timescale_parameters.csv"),
                 row.names = FALSE)

  boot_df <- summary_df[, c("path", "boot_ci_low", "boot_ci_high", "boot_sig")]
  safe_write_csv(boot_df,
                 file.path(pixelwise_dir, "SEM_dual_timescale_pixelwise_bootstrap_ci.csv"),
                 row.names = FALSE)

  safe_write_csv(r2_summary,
                 file.path(pixelwise_dir, "SEM_dual_timescale_R2_detail.csv"),
                 row.names = FALSE)
  safe_write_csv(r2_summary[, c("variable", "R2")],
                 file.path(pixelwise_dir, "SEM_dual_timescale_R2.csv"),
                 row.names = FALSE)

  safe_write_csv(filter_df,
                 file.path(pixelwise_dir, "SEM_dual_timescale_pixelwise_outlier_filtering.csv"),
                 row.names = FALSE)

  summary_path <- file.path(pixelwise_dir, "SEM_dual_timescale_summary.txt")
  if (should_write(summary_path)) {
    sink(summary_path)
    cat(sprintf("像元时间序列SEM汇总 (%s)\n\n", label))
    print(summary_df)
    cat("\n输出异常值过滤统计（像元级）:\n")
    print(filter_df)
    sink()
  } else {
    cat(sprintf("  [skip] %s\n", summary_path))
  }
}

run_pixel_sem_lavaan_compare <- function(years, mask_r, fixed_window_length_r,
                                         pixelwise_dir_ours, pixelwise_dir_other, compare_dir) {
  cat("\n=== SEM分析（像元时间序列，lavaan）===\n")
  cat("⚠️ 注意：此模式计算量极大，可能需要数小时！\n\n")

  files <- list(
    TR_fixed_window = file.path(DECOMP_DIR, sprintf("TR_fixed_window_%d.tif", years)),
    EOS = file.path(PHENO_DIR, sprintf("eos_gpp_%d.tif", years)),
    Fixed_GPPrate = file.path(DECOMP_DIR, sprintf("Fixed_GPPrate_%d.tif", years)),
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
  other_min_years <- max(1, ceiling(OTHER_MIN_VALID_YEAR_FRAC * length(years)))
  if (length(years) < min_years) {
    stop("像元时间序列SEM所需的完整输入年份不足")
  }
  cat(sprintf("  'Ours'方法最小年份: %d (%.0f%%)\n", min_years, MIN_VALID_YEAR_FRAC * 100))
  cat(sprintf("  'Other'方法最小年份: %d (%.0f%%)\n", other_min_years, OTHER_MIN_VALID_YEAR_FRAC * 100))

  stacks <- lapply(files, stack)
  na_values <- lapply(files, function(paths) NAvalue(raster(paths[1])))

  template <- raster(files$TR_fixed_window[1])
  bs <- blockSize(template)
  mask_vals_all <- getValues(mask_r)
  total_valid_cells <- sum(!is.na(mask_vals_all))
  processed_cells <- 0

  zero_coef <- setNames(rep(0, length(coef_names)), coef_names)
  zero_r2 <- setNames(rep(0, length(r2_names)), r2_names)
  zero_filter <- c(coef_extreme = 0, p_invalid = 0, r2_invalid = 0)

  # "ours" stats
  coef_sum <- zero_coef
  coef_sumsq <- zero_coef
  coef_count <- zero_coef
  sig_sum <- zero_coef
  sig_sumsq <- zero_coef
  sig_count <- zero_coef
  r2_sum <- zero_r2
  r2_sumsq <- zero_r2
  r2_count <- zero_r2
  filter_info <- zero_filter
  coef_samples_list <- list()

  # "other" stats (GFI + p-value masking + 100% data requirement)
  other_sig_sum <- zero_coef
  other_sig_sumsq <- zero_coef
  other_sig_count <- zero_coef
  other_total_count <- zero_coef
  other_data_incomplete <- 0L  # 数据完整性<100%的像元数
  other_gfi_fail <- 0L
  other_fit_fail <- 0L
  other_r2_sum <- zero_r2
  other_r2_sumsq <- zero_r2
  other_r2_count <- zero_r2
  other_coef_samples_list <- list()

  process_block <- function(i) {
    row <- bs$row[i]
    nrows <- bs$nrows[i]

    mask_vals <- getValues(mask_r, row = row, nrows = nrows)
    valid_mask_cells <- which(!is.na(mask_vals))
    if (length(valid_mask_cells) == 0) {
      return(list(
        coef_sum = zero_coef, coef_sumsq = zero_coef, coef_count = zero_coef,
        sig_sum = zero_coef, sig_sumsq = zero_coef, sig_count = zero_coef,
        r2_sum = zero_r2, r2_sumsq = zero_r2, r2_count = zero_r2,
        filter_info = zero_filter,
        other_sig_sum = zero_coef, other_sig_sumsq = zero_coef,
        other_sig_count = zero_coef, other_total_count = zero_coef,
        other_data_incomplete = 0L, other_gfi_fail = 0L, other_fit_fail = 0L,
        other_r2_sum = zero_r2, other_r2_sumsq = zero_r2, other_r2_count = zero_r2,
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
        coef_sum = zero_coef, coef_sumsq = zero_coef, coef_count = zero_coef,
        sig_sum = zero_coef, sig_sumsq = zero_coef, sig_count = zero_coef,
        r2_sum = zero_r2, r2_sumsq = zero_r2, r2_count = zero_r2,
        filter_info = zero_filter,
        other_sig_sum = zero_coef, other_sig_sumsq = zero_coef,
        other_sig_count = zero_coef, other_total_count = zero_coef,
        other_data_incomplete = 0L, other_gfi_fail = 0L, other_fit_fail = 0L,
        other_r2_sum = zero_r2, other_r2_sumsq = zero_r2, other_r2_count = zero_r2,
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

    other_sig_sum_block <- zero_coef
    other_sig_sumsq_block <- zero_coef
    other_sig_count_block <- zero_coef
    other_total_count_block <- zero_coef
    other_data_incomplete_block <- 0L
    other_gfi_fail_block <- 0L
    other_fit_fail_block <- 0L
    other_r2_sum_block <- zero_r2
    other_r2_sumsq_block <- zero_r2
    other_r2_count_block <- zero_r2
    other_coef_samples_block <- matrix(NA_real_, nrow = length(valid_cells), ncol = length(coef_names))
    other_coef_sample_rows <- 0L

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

      df_std <- data.frame(
        Fixed_Trate = scale_vec(df$Fixed_Trate),
        EOS = scale_vec(df$EOS),
        Fixed_GPPrate = scale_vec(df$Fixed_GPPrate),
        P_pre = scale_vec(df$P_pre),
        T_pre = scale_vec(df$T_pre),
        SW_pre = scale_vec(df$SW_pre),
        P_season = scale_vec(df$P_season),
        T_season = scale_vec(df$T_season),
        SW_season = scale_vec(df$SW_season)
      )

      if (any(!is.finite(as.matrix(df_std)))) {
        next
      }

      fit_res <- extract_lavaan_results(df_std)
      if (is.null(fit_res)) {
        other_fit_fail_block <- other_fit_fail_block + 1L
        next
      }

      vals <- fit_res$vals
      p_vals <- fit_res$p_vals
      r2_vals <- fit_res$r2_vals

      # === Our screening ===
      filter_res <- filter_sem_outputs(vals, p_vals, r2_vals)
      vals_f <- filter_res$vals
      p_vals_f <- filter_res$p_vals
      r2_vals_f <- filter_res$r2_vals
      filter_info_block <- filter_info_block + filter_res$info

      coef_sample_rows <- coef_sample_rows + 1
      coef_samples_block[coef_sample_rows, ] <- vals_f

      for (k in seq_along(coef_names)) {
        if (is.finite(vals_f[k])) {
          coef_sum_block[coef_names[k]] <- coef_sum_block[coef_names[k]] + vals_f[k]
          coef_sumsq_block[coef_names[k]] <- coef_sumsq_block[coef_names[k]] + vals_f[k]^2
          coef_count_block[coef_names[k]] <- coef_count_block[coef_names[k]] + 1
          if (is.finite(p_vals_f[k]) && p_vals_f[k] < 0.05) {
            sig_sum_block[coef_names[k]] <- sig_sum_block[coef_names[k]] + vals_f[k]
            sig_sumsq_block[coef_names[k]] <- sig_sumsq_block[coef_names[k]] + vals_f[k]^2
            sig_count_block[coef_names[k]] <- sig_count_block[coef_names[k]] + 1
          }
        }
      }

      for (k in seq_along(r2_names)) {
        if (is.finite(r2_vals_f[k])) {
          r2_sum_block[r2_names[k]] <- r2_sum_block[r2_names[k]] + r2_vals_f[k]
          r2_sumsq_block[r2_names[k]] <- r2_sumsq_block[r2_names[k]] + r2_vals_f[k]^2
          r2_count_block[r2_names[k]] <- r2_count_block[r2_names[k]] + 1
        }
      }

      # === Other screening (100% data + GFI + p-mask) ===
      # Step 1: 检查数据完整性（Other方法要求100%）
      if (sum(valid_years) < other_min_years) {
        other_data_incomplete_block <- other_data_incomplete_block + 1L
        next
      }

      # Step 2: 检查GFI
      gfi_val <- NA_real_
      if (!is.null(fit_res$fit_measures) && "gfi" %in% names(fit_res$fit_measures)) {
        gfi_val <- fit_res$fit_measures[["gfi"]]
      }
      if (!is.finite(gfi_val) || gfi_val < OTHER_GFI_MIN) {
        other_gfi_fail_block <- other_gfi_fail_block + 1L
        next
      }

      other_vals <- vals
      other_vals[!(is.finite(p_vals) & p_vals < OTHER_P_THRESHOLD)] <- NA_real_
      other_coef_sample_rows <- other_coef_sample_rows + 1L
      other_coef_samples_block[other_coef_sample_rows, ] <- other_vals

      for (k in seq_along(r2_names)) {
        if (is.finite(r2_vals[k])) {
          other_r2_sum_block[r2_names[k]] <- other_r2_sum_block[r2_names[k]] + r2_vals[k]
          other_r2_sumsq_block[r2_names[k]] <- other_r2_sumsq_block[r2_names[k]] + r2_vals[k]^2
          other_r2_count_block[r2_names[k]] <- other_r2_count_block[r2_names[k]] + 1
        }
      }

      for (k in seq_along(coef_names)) {
        if (is.finite(vals[k])) {
          other_total_count_block[coef_names[k]] <- other_total_count_block[coef_names[k]] + 1
          if (is.finite(p_vals[k]) && p_vals[k] < OTHER_P_THRESHOLD) {
            other_sig_sum_block[coef_names[k]] <- other_sig_sum_block[coef_names[k]] + vals[k]
            other_sig_sumsq_block[coef_names[k]] <- other_sig_sumsq_block[coef_names[k]] + vals[k]^2
            other_sig_count_block[coef_names[k]] <- other_sig_count_block[coef_names[k]] + 1
          }
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
      other_sig_sum = other_sig_sum_block,
      other_sig_sumsq = other_sig_sumsq_block,
      other_sig_count = other_sig_count_block,
      other_total_count = other_total_count_block,
      other_data_incomplete = other_data_incomplete_block,
      other_gfi_fail = other_gfi_fail_block,
      other_fit_fail = other_fit_fail_block,
      other_r2_sum = other_r2_sum_block,
      other_r2_sumsq = other_r2_sumsq_block,
      other_r2_count = other_r2_count_block,
      coef_samples = if (coef_sample_rows > 0) {
        coef_samples_block[seq_len(coef_sample_rows), , drop = FALSE]
      } else {
        NULL
      },
      other_coef_samples = if (other_coef_sample_rows > 0) {
        other_coef_samples_block[seq_len(other_coef_sample_rows), , drop = FALSE]
      } else {
        NULL
      },
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

    other_sig_sum <<- other_sig_sum + res$other_sig_sum
    other_sig_sumsq <<- other_sig_sumsq + res$other_sig_sumsq
    other_sig_count <<- other_sig_count + res$other_sig_count
    other_total_count <<- other_total_count + res$other_total_count
    other_data_incomplete <<- other_data_incomplete + res$other_data_incomplete
    other_gfi_fail <<- other_gfi_fail + res$other_gfi_fail
    other_fit_fail <<- other_fit_fail + res$other_fit_fail
    other_r2_sum <<- other_r2_sum + res$other_r2_sum
    other_r2_sumsq <<- other_r2_sumsq + res$other_r2_sumsq
    other_r2_count <<- other_r2_count + res$other_r2_count
    if (!is.null(res$coef_samples)) {
      coef_samples_list[[length(coef_samples_list) + 1]] <<- res$coef_samples
    }
    if (!is.null(res$other_coef_samples)) {
      other_coef_samples_list[[length(other_coef_samples_list) + 1]] <<- res$other_coef_samples
    }

    processed_cells <<- processed_cells + res$processed_cells
    blocks_done <<- blocks_done + 1
  }

  if (PARALLEL_ENABLE && PARALLEL_CORES > 1) {
    cat(sprintf("  启用并行块处理: %d cores\n", PARALLEL_CORES))
    cl <- makeCluster(PARALLEL_CORES)
    clusterEvalQ(cl, library(raster))
    clusterEvalQ(cl, library(lavaan))
    clusterExport(
      cl,
      c("bs", "mask_r", "fixed_window_length_r", "stacks", "na_values", "years",
        "min_years", "other_min_years", "DETREND_PIXEL_ENABLE", "coef_names", "r2_names",
        "sanitize_values", "detrend_series", "scale_vec", "filter_sem_outputs",
        "FILTER_SEM_OUTLIERS", "SEM_COEF_ABS_MAX", "SEM_P_MIN", "SEM_P_MAX",
        "SEM_R2_MIN", "SEM_R2_MAX", "zero_coef", "zero_r2", "process_block",
        "extract_lavaan_results", "sem_model", "OTHER_GFI_MIN",
        "OTHER_P_THRESHOLD", "zero_filter",
        "safe_mediation_ratio", "MEDIATION_DENOM_EPS",
        "delta_var_two", "delta_var_three", "delta_p_two", "delta_p_three",
        "delta_p_ratio"),
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

  # === Summaries ===
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

  coef_samples_all <- NULL
  if (length(coef_samples_list) > 0) {
    coef_samples_all <- do.call(rbind, coef_samples_list)
  }

  other_coef_samples_all <- NULL
  if (length(other_coef_samples_list) > 0) {
    other_coef_samples_all <- do.call(rbind, other_coef_samples_list)
  }

  boot_ci_low <- rep(NA_real_, length(coef_names))
  boot_ci_high <- rep(NA_real_, length(coef_names))
  boot_sig <- rep(NA, length(coef_names))
  boot_ci_low_other <- rep(NA_real_, length(coef_names))
  boot_ci_high_other <- rep(NA_real_, length(coef_names))
  boot_sig_other <- rep(NA, length(coef_names))

  if (PIXELWISE_BOOTSTRAP_ENABLE && !is.null(coef_samples_all)) {
    cat(sprintf("\n  像元均值 bootstrap CI (Ours): %d 次\n", PIXELWISE_BOOTSTRAP_N))
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

  if (PIXELWISE_BOOTSTRAP_ENABLE && !is.null(other_coef_samples_all)) {
    cat(sprintf("\n  像元均值 bootstrap CI (Other): %d 次\n", PIXELWISE_BOOTSTRAP_N))
    for (k in seq_along(coef_names)) {
      ci <- bootstrap_ci_mean(
        other_coef_samples_all[, k],
        PIXELWISE_BOOTSTRAP_N,
        seed = PIXELWISE_BOOTSTRAP_SEED + 10000 + k
      )
      boot_ci_low_other[k] <- ci["ci_low"]
      boot_ci_high_other[k] <- ci["ci_high"]
      if (is.finite(boot_ci_low_other[k]) && is.finite(boot_ci_high_other[k])) {
        boot_sig_other[k] <- !(boot_ci_low_other[k] <= 0 && boot_ci_high_other[k] >= 0)
      }
    }
  }

  summary_ours <- data.frame(
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

  other_mean <- other_sig_sum / other_sig_count
  other_sd <- sqrt(pmax(0, other_sig_sumsq / other_sig_count - other_mean^2))
  other_sig_frac <- other_sig_count / other_total_count

  other_mean[!is.finite(other_mean)] <- NA_real_
  other_sd[!is.finite(other_sd)] <- NA_real_
  other_sig_frac[!is.finite(other_sig_frac)] <- NA_real_

  summary_other <- data.frame(
    path = coef_names,
    mean = other_mean,
    sd = other_sd,
    n = other_total_count,
    sig_n = other_sig_count,
    sig_frac = other_sig_frac,
    mean_sig = other_mean,
    sd_sig = other_sd,
    boot_ci_low = boot_ci_low_other,
    boot_ci_high = boot_ci_high_other,
    boot_sig = boot_sig_other
  )

  compare_df <- data.frame(
    path = coef_names,
    mean_ours = mean_vals,
    mean_other_sig = other_mean,
    delta_mean = other_mean - mean_vals,
    sig_frac_ours = sig_frac,
    sig_frac_other = other_sig_frac,
    n_ours = coef_count,
    n_other_total = other_total_count
  )

  r2_summary_ours <- build_r2_summary(r2_sum, r2_sumsq, r2_count)
  r2_summary_other <- build_r2_summary(other_r2_sum, other_r2_sumsq, other_r2_count)

  filter_df_ours <- data.frame(
    metric = c("coef_extreme", "p_invalid", "r2_invalid", "data_incomplete", "gfi_fail", "fit_fail"),
    count = c(filter_info["coef_extreme"], filter_info["p_invalid"], filter_info["r2_invalid"], 0, 0, 0)
  )
  filter_df_other <- data.frame(
    metric = c("coef_extreme", "p_invalid", "r2_invalid", "data_incomplete", "gfi_fail", "fit_fail"),
    count = c(0, 0, 0, other_data_incomplete, other_gfi_fail, other_fit_fail)
  )

  # === Write outputs ===
  write_pixelwise_outputs(pixelwise_dir_ours, summary_ours, r2_summary_ours, filter_df_ours, "ours")
  write_pixelwise_outputs(pixelwise_dir_other, summary_other, r2_summary_other, filter_df_other, "other")

  safe_write_csv(compare_df,
                 file.path(compare_dir, "SEM_dual_timescale_pixel_time_series_compare.csv"),
                 row.names = FALSE)

  filter_df <- data.frame(
    metric = c("coef_extreme", "p_invalid", "r2_invalid", "data_incomplete", "gfi_fail", "fit_fail"),
    count = c(filter_info["coef_extreme"], filter_info["p_invalid"], filter_info["r2_invalid"],
              other_data_incomplete, other_gfi_fail, other_fit_fail)
  )
  safe_write_csv(filter_df,
                 file.path(compare_dir, "SEM_dual_timescale_pixel_time_series_filtering.csv"),
                 row.names = FALSE)

  # === Print summary to console ===
  cat("\n=== 筛选统计汇总 ===\n")
  cat(sprintf("'Ours'方法:\n"))
  cat(sprintf("  - 总处理像元数: %d\n", processed_cells))
  cat(sprintf("  - 异常值过滤: %d (系数极值: %d, p值无效: %d, R2无效: %d)\n",
              sum(filter_info), filter_info["coef_extreme"],
              filter_info["p_invalid"], filter_info["r2_invalid"]))
  cat(sprintf("  - 有效像元数: %d (%.1f%%)\n",
              as.integer(mean(coef_count)),
              100 * mean(coef_count) / processed_cells))
  cat(sprintf("  - 显著像元数: %d (%.1f%%)\n",
              as.integer(mean(sig_count)),
              100 * mean(sig_count) / mean(coef_count)))

  cat(sprintf("\n'Other'方法 (100%% data + GFI≥%.2f + p<%.2f):\n",
              OTHER_GFI_MIN, OTHER_P_THRESHOLD))
  cat(sprintf("  - 总处理像元数: %d\n", processed_cells))
  cat(sprintf("  - 数据不完整过滤: %d (要求100%%完整)\n", other_data_incomplete))
  cat(sprintf("  - GFI筛选过滤: %d (GFI<%.2f)\n", other_gfi_fail, OTHER_GFI_MIN))
  cat(sprintf("  - 拟合失败过滤: %d\n", other_fit_fail))
  cat(sprintf("  - 通过筛选像元数: %d (%.1f%%)\n",
              as.integer(mean(other_total_count)),
              100 * mean(other_total_count) / processed_cells))
  cat(sprintf("  - 显著像元数: %d (%.1f%%)\n",
              as.integer(mean(other_sig_count)),
              100 * mean(other_sig_count) / mean(other_total_count)))

  cat(sprintf("\n像元数对比:\n"))
  cat(sprintf("  - 'Ours'有效像元: %d\n", as.integer(mean(coef_count))))
  cat(sprintf("  - 'Other'通过筛选: %d\n", as.integer(mean(other_total_count))))
  cat(sprintf("  - 差异: %d (%.1f%%)\n",
              as.integer(mean(coef_count) - mean(other_total_count)),
              100 * (mean(coef_count) - mean(other_total_count)) / mean(coef_count)))

  list(summary_ours = summary_ours, summary_other = summary_other, compare = compare_df)
}

main <- function() {
  cat("\n======================================================================\n")
  cat("双时间尺度SEM分析 - Lavaan像元级对比（Fixed Window）\n")
  cat("======================================================================\n")

  cat("[环境检测]\n")
  cat(sprintf("  操作系统类型: %s\n", .Platform$OS.type))
  cat(sprintf("  R版本: %s\n", R.version.string))
  cat(sprintf("  根目录: %s\n", ROOT))

  cat("\n[路径检查]\n")
  check_dir(ROOT, "Root directory")
  check_dir(PHENO_DIR, "Phenology directory")
  check_dir(DECOMP_DIR, "Decomposition FixedWindow directory")
  check_dir(DERIVED_DIR, "SEM derived cache directory")
  check_file(MASK_FILE, "Mask file")
  check_file(file.path(DECOMP_DIR, "Fixed_Window_Length.tif"), "Fixed_Window_Length file")

  cat("\n[基础信息]\n")
  cat(sprintf("  年份范围: %d-%d\n", YEAR_START, YEAR_END))
  cat("  分析模式: pixel_time_series\n")

  years <- YEAR_START:YEAR_END

  cat("\n[缓存文件检查]\n")
  missing_files <- 0L
  for (year in years) {
    files_to_check <- c(
      file.path(DECOMP_DIR, sprintf("TR_fixed_window_%d.tif", year)),
      file.path(PHENO_DIR, sprintf("eos_gpp_%d.tif", year)),
      file.path(DECOMP_DIR, sprintf("Fixed_GPPrate_%d.tif", year)),
      file.path(DERIVED_DIR, sprintf("P_pre_%d.tif", year)),
      file.path(DERIVED_DIR, sprintf("T_pre_%d.tif", year)),
      file.path(DERIVED_DIR, sprintf("SW_pre_%d.tif", year)),
      file.path(DERIVED_DIR, sprintf("P_season_%d.tif", year)),
      file.path(DERIVED_DIR, sprintf("T_season_%d.tif", year)),
      file.path(DERIVED_DIR, sprintf("SW_season_%d.tif", year))
    )
    for (f in files_to_check) {
      if (!file.exists(f)) {
        cat(sprintf("  ✗ 缺失: %s\n", basename(f)))
        missing_files <- missing_files + 1L
      }
    }
  }
  if (missing_files > 0) {
    stop(sprintf("缓存文件缺失: %d 个，请先运行05b生成缓存。", missing_files))
  }
  cat(sprintf("  ✓ 所有 %d 年的缓存文件验证通过\n", length(years)))

  mask_r <- raster(MASK_FILE)
  fixed_window_length_r <- raster(file.path(DECOMP_DIR, "Fixed_Window_Length.tif"))
  fixed_window_length_r <- mask_raster(fixed_window_length_r, mask_r)
  fixed_window_length_r[fixed_window_length_r <= 0] <- NA

  cat("\n[掩膜诊断]\n")
  mask_vals <- getValues(mask_r)
  n_total <- length(mask_vals)
  n_valid <- sum(is.finite(mask_vals) & mask_vals > 0)
  cat(sprintf("  栅格尺寸: %d x %d\n", nrow(mask_r), ncol(mask_r)))
  cat(sprintf("  总像元数: %d\n", n_total))
  cat(sprintf("  有效像元数 (>0): %d (%.2f%%)\n", n_valid, 100 * n_valid / n_total))

  cat("\n[网格一致性检查]\n")
  template_path <- if (file.exists(TEMPLATE_FILE)) {
    TEMPLATE_FILE
  } else {
    file.path(DECOMP_DIR, sprintf("TR_fixed_window_%d.tif", YEAR_START))
  }
  template_r <- raster(template_path)
  sample_rasters <- list(
    fixed_window_length_r,
    raster(file.path(DECOMP_DIR, sprintf("TR_fixed_window_%d.tif", YEAR_START))),
    raster(file.path(PHENO_DIR, sprintf("eos_gpp_%d.tif", YEAR_START))),
    raster(file.path(DECOMP_DIR, sprintf("Fixed_GPPrate_%d.tif", YEAR_START))),
    raster(file.path(DERIVED_DIR, sprintf("P_pre_%d.tif", YEAR_START)))
  )
  tryCatch({
    do.call(compareRaster, c(list(template_r), sample_rasters, list(stopiffalse = TRUE)))
    cat("  ✓ 栅格对齐检查通过\n")
  }, error = function(e) {
    stop(sprintf("网格一致性检查失败: %s", conditionMessage(e)))
  })

  if (RUN_BOTH_DETREND) {
    runs <- list(
      list(enable = FALSE, label = "raw", suffix = ""),
      list(enable = TRUE, label = "detrended", suffix = "_detrended")
    )
  } else {
    runs <- list(
      list(
        enable = DETREND_PIXEL_ENABLE,
        label = ifelse(DETREND_PIXEL_ENABLE, "detrended", "raw"),
        suffix = ifelse(DETREND_PIXEL_ENABLE, "_detrended", "")
      )
    )
  }

  for (pr in runs) {
    if (RUN_MODE == "skip" && outputs_ready(pr$suffix)) {
      cat(sprintf("\n=== 像元级SEM（lavaan）[%s] ===\n", pr$label))
      cat("  ✓ 输出齐全，已跳过\n")
      next
    }
    dirs_ours <- set_output_dirs(pr$suffix, "Ours")
    dirs_other <- set_output_dirs(pr$suffix, "Other")
    DETREND_PIXEL_ENABLE <<- pr$enable
    cat(sprintf("\n=== 像元级SEM（lavaan）[%s] ===\n", pr$label))
    run_pixel_sem_lavaan_compare(
      years,
      mask_r,
      fixed_window_length_r,
      dirs_ours$pixelwise_dir,
      dirs_other$pixelwise_dir,
      dirs_ours$compare_dir
    )
  }
}

if (!interactive()) {
  main()
} else {
  cat("\n请手动运行 main() 以执行 lavaan 像元级SEM对比分析。\n")
}
