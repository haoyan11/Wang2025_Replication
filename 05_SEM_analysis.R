#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Module 06: 结构方程模型（SEM）分析
# 使用lavaan包实现Wang (2025)的SEM分析

library(lavaan)
library(semPlot)
library(raster)
library(tidyverse)

# ==================== 全局配置 ====================
ROOT <- "I:/F/Data4/Wang2025_Analysis"
DATA_DIR <- file.path(ROOT, "SEM_Data")
OUTPUT_DIR <- file.path(ROOT, "SEM_Results")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

YEAR_START <- 1982
YEAR_END <- 2018

# ==================== 数据准备 ====================
prepare_sem_data <- function() {
  cat("\n=== 准备SEM数据 ===\n")

  years <- YEAR_START:YEAR_END

  # 初始化数据框
  sem_data <- data.frame()

  # 逐年读取栅格数据并提取均值
  for (year in years) {
    cat(sprintf("  处理年份: %d\n", year))

    # 读取各变量的年度数据
    TRc <- raster(file.path(ROOT, "TRc_annual", sprintf("TRc_%d.tif", year)))
    TRpheno <- raster(file.path(ROOT, "Decomposition", sprintf("TRpheno_%d.tif", year)))
    TRproduct <- raster(file.path(ROOT, "Decomposition", sprintf("TRproduct_%d.tif", year)))

    SOS <- raster(file.path(ROOT, "../Phenology_Output_1/SIF_phenology", sprintf("SOS_%d.tif", year)))
    POS <- raster(file.path(ROOT, "../Phenology_Output_1/SIF_phenology", sprintf("POS_%d.tif", year)))

    SIF <- raster(file.path(ROOT, "../SIF_Data/CSIF_annual", sprintf("SIF_annual_%d.tif", year)))
    SM <- raster(file.path(ROOT, "../Meteorological Data/GLEAM/Annual/SMrz", sprintf("SMrz_%d.tif", year)))

    # 计算LSP
    LSP <- POS - SOS

    # 读取掩膜
    mask <- raster(file.path(ROOT, "masks/combined_mask.tif"))

    # 应用掩膜并提取有效值
    TRc_masked <- mask(TRc, mask, maskvalue=0)
    TRpheno_masked <- mask(TRpheno, mask, maskvalue=0)
    TRproduct_masked <- mask(TRproduct, mask, maskvalue=0)
    SOS_masked <- mask(SOS, mask, maskvalue=0)
    LSP_masked <- mask(LSP, mask, maskvalue=0)
    SIF_masked <- mask(SIF, mask, maskvalue=0)
    SM_masked <- mask(SM, mask, maskvalue=0)

    # 计算空间平均值（或使用逐像元数据）
    year_data <- data.frame(
      year = year,
      TRc = cellStats(TRc_masked, mean, na.rm=TRUE),
      TRpheno = cellStats(TRpheno_masked, mean, na.rm=TRUE),
      TRproduct = cellStats(TRproduct_masked, mean, na.rm=TRUE),
      SOS = cellStats(SOS_masked, mean, na.rm=TRUE),
      LSP = cellStats(LSP_masked, mean, na.rm=TRUE),
      SIF = cellStats(SIF_masked, mean, na.rm=TRUE),
      SM = cellStats(SM_masked, mean, na.rm=TRUE)
    )

    sem_data <- rbind(sem_data, year_data)
  }

  # 标准化数据
  sem_data_std <- sem_data %>%
    mutate(across(-year, ~scale(.)[,1]))

  write.csv(sem_data, file.path(DATA_DIR, "sem_data_raw.csv"), row.names=FALSE)
  write.csv(sem_data_std, file.path(DATA_DIR, "sem_data_standardized.csv"), row.names=FALSE)

  cat("✓ 数据准备完成\n")
  return(sem_data_std)
}

# ==================== 原版SEM（Wang 2025） ====================
run_original_SEM <- function(data) {
  cat("\n=== 原版SEM分析（Wang 2025） ===\n")

  # 定义模型（原文路径）
  # 注意：原文存在问题（缺少SOS→SIF路径）
  model_original <- '
    # 直接效应
    TRproduct ~ b1*LSP + b2*SIF + b3*SM
    SIF ~ c1*SM

    # 间接效应（通过SIF）
    indirect_LSP_SIF := b1 * c1
    indirect_SM_SIF := b2 * c1

    # 总效应
    total_LSP := b1 + indirect_LSP_SIF
    total_SM := b3 + indirect_SM_SIF
  '

  # 拟合模型
  fit_original <- sem(model_original, data=data, estimator="MLR")

  # 输出结果
  cat("\n原版SEM拟合结果：\n")
  summary(fit_original, fit.measures=TRUE, standardized=TRUE, rsquare=TRUE)

  # 保存参数估计
  params <- parameterEstimates(fit_original, standardized=TRUE)
  write.csv(params, file.path(OUTPUT_DIR, "SEM_original_parameters.csv"), row.names=FALSE)

  # 拟合指标
  fit_measures <- fitMeasures(fit_original, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))
  write.csv(t(as.data.frame(fit_measures)), file.path(OUTPUT_DIR, "SEM_original_fitindices.csv"))

  # 绘制路径图
  pdf(file.path(OUTPUT_DIR, "SEM_original_pathdiagram.pdf"), width=10, height=8)
  semPaths(fit_original,
           what="std",
           edge.label.cex=1.2,
           curvePivot=TRUE,
           layout="tree2",
           style="lisrel",
           edge.color="black",
           nodeLabels=c("TRproduct", "LSP", "SIF", "SM"),
           sizeMan=10,
           residuals=FALSE,
           exoCov=FALSE)
  title("Original SEM (Wang 2025)", line=3)
  dev.off()

  cat("✓ 原版SEM分析完成\n")
  return(fit_original)
}

# ==================== 改进版SEM ====================
run_improved_SEM <- function(data) {
  cat("\n=== 改进版SEM分析（补充材料） ===\n")

  # 定义改进模型（添加缺失路径）
  model_improved <- '
    # 直接效应
    TRproduct ~ b1*LSP + b2*SIF + b3*SM
    SIF ~ c1*SM + c2*SOS  # 添加SOS→SIF路径
    LSP ~ d1*SOS          # 添加SOS→LSP路径

    # 间接效应
    indirect_SOS_LSP_TRproduct := d1 * b1
    indirect_SOS_SIF_TRproduct := c2 * b2
    indirect_SM_SIF_TRproduct := c1 * b2

    # 总效应
    total_SOS := indirect_SOS_LSP_TRproduct + indirect_SOS_SIF_TRproduct
    total_SM := b3 + indirect_SM_SIF_TRproduct
    total_LSP := b1
  '

  # 拟合模型
  fit_improved <- sem(model_improved, data=data, estimator="MLR")

  # 输出结果
  cat("\n改进版SEM拟合结果：\n")
  summary(fit_improved, fit.measures=TRUE, standardized=TRUE, rsquare=TRUE)

  # 保存参数估计
  params <- parameterEstimates(fit_improved, standardized=TRUE)
  write.csv(params, file.path(OUTPUT_DIR, "SEM_improved_parameters.csv"), row.names=FALSE)

  # 拟合指标
  fit_measures <- fitMeasures(fit_improved, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))
  write.csv(t(as.data.frame(fit_measures)), file.path(OUTPUT_DIR, "SEM_improved_fitindices.csv"))

  # 绘制路径图
  pdf(file.path(OUTPUT_DIR, "SEM_improved_pathdiagram.pdf"), width=10, height=8)
  semPaths(fit_improved,
           what="std",
           edge.label.cex=1.2,
           curvePivot=TRUE,
           layout="tree2",
           style="lisrel",
           edge.color="black",
           nodeLabels=c("TRproduct", "LSP", "SIF", "SM", "SOS"),
           sizeMan=10,
           residuals=FALSE,
           exoCov=FALSE)
  title("Improved SEM (with missing pathways)", line=3)
  dev.off()

  cat("✓ 改进版SEM分析完成\n")
  return(fit_improved)
}

# ==================== 模型比较 ====================
compare_models <- function(fit_original, fit_improved) {
  cat("\n=== 模型比较 ===\n")

  # AIC/BIC比较
  comparison <- data.frame(
    Model = c("Original", "Improved"),
    AIC = c(AIC(fit_original), AIC(fit_improved)),
    BIC = c(BIC(fit_original), BIC(fit_improved)),
    CFI = c(fitMeasures(fit_original, "cfi"), fitMeasures(fit_improved, "cfi")),
    RMSEA = c(fitMeasures(fit_original, "rmsea"), fitMeasures(fit_improved, "rmsea")),
    SRMR = c(fitMeasures(fit_original, "srmr"), fitMeasures(fit_improved, "srmr"))
  )

  print(comparison)
  write.csv(comparison, file.path(OUTPUT_DIR, "SEM_model_comparison.csv"), row.names=FALSE)

  # 卡方差异检验
  chi_diff <- anova(fit_original, fit_improved)
  print(chi_diff)
  write.csv(chi_diff, file.path(OUTPUT_DIR, "SEM_chisq_difference.csv"))

  cat("✓ 模型比较完成\n")
}

# ==================== VIF诊断 ====================
calculate_VIF <- function(data) {
  cat("\n=== VIF多重共线性诊断 ===\n")

  # 对TRproduct进行多元回归
  lm_model <- lm(TRproduct ~ LSP + SIF + SM + SOS, data=data)

  # 计算VIF
  library(car)
  vif_values <- vif(lm_model)

  cat("\nVIF值：\n")
  print(vif_values)

  vif_df <- data.frame(
    Variable = names(vif_values),
    VIF = vif_values
  )
  write.csv(vif_df, file.path(OUTPUT_DIR, "VIF_diagnostics.csv"), row.names=FALSE)

  # 判断标准：VIF > 10 表示严重多重共线性
  if (any(vif_values > 10)) {
    cat("\n⚠ 警告：存在严重多重共线性（VIF > 10）\n")
  } else if (any(vif_values > 5)) {
    cat("\n⚠ 注意：存在中度多重共线性（VIF > 5）\n")
  } else {
    cat("\n✓ VIF检验通过，无严重多重共线性\n")
  }

  return(vif_values)
}

# ==================== 主程序 ====================
main <- function() {
  cat("\n" paste0(rep("=", 70), collapse=""), "\n")
  cat("结构方程模型（SEM）分析\n")
  cat(paste0(rep("=", 70), collapse=""), "\n")

  # 1. 准备数据
  sem_data <- prepare_sem_data()

  # 2. 原版SEM
  fit_original <- run_original_SEM(sem_data)

  # 3. 改进版SEM
  fit_improved <- run_improved_SEM(sem_data)

  # 4. 模型比较
  compare_models(fit_original, fit_improved)

  # 5. VIF诊断
  vif_values <- calculate_VIF(sem_data)

  cat("\n", paste0(rep("=", 70), collapse=""), "\n")
  cat("✓ SEM分析完成！\n")
  cat(sprintf("输出目录: %s\n", OUTPUT_DIR))
  cat(paste0(rep("=", 70), collapse=""), "\n")
}

# 运行主程序
if (!interactive()) {
  main()
}
