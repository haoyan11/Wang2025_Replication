# 📚 文档索引

**最后更新**: 2025-12-24

---

## 🚀 快速入门

### 新手必读
1. **[项目说明](../README.md)** - 项目概述、背景和主要功能（根目录）
2. **[快速开始](00_README/START_HERE.md)** - 5分钟快速上手指南
3. **[依赖安装](00_README/requirements.txt)** - Python环境配置

---

## 📖 使用指南

### [01_guides/](01_guides/)
- **[数据配置说明.md](01_guides/数据配置说明.md)** - 数据路径配置、格式要求、准备步骤
- **[CHECKLIST.md](01_guides/CHECKLIST.md)** - 质量检查清单、验证步骤

---

## 📝 变更记录

### [02_changelogs/](02_changelogs/)
- **[CHANGELOG.md](02_changelogs/CHANGELOG.md)** - 完整版本历史（v1.0.0 → v2.0.0）
- **[今日修改总结_20251224.md](02_changelogs/今日修改总结_20251224.md)** ⭐ **最新**
  - 新增Timing/Shape分解方法
  - 创建T物候版本
  - 添加Trate分析
  - 16页详细记录
- **[技术修正记录.md](02_changelogs/技术修正记录.md)** ⭐ **整合文档**
  - TRc计算v1.3.1性能优化（18倍提速）
  - ΔSOS标准异常定义修正
  - 土壤湿度数据格式兼容性修复
  - 3个技术修正整合到一个文档

---

## 📊 分析报告

### [03_analysis/](03_analysis/)
- **[检查统计结果差异.md](03_analysis/检查统计结果差异.md)** - 回归结果与Wang 2025论文对比诊断

---

## 🎯 阶段总结

### [04_summaries/](04_summaries/)
- **[完成总结.md](04_summaries/完成总结.md)** - 项目整体完成情况、操作指南

---

## 🔧 核心脚本说明

### 数据准备（预处理）
| 脚本 | 功能 | 输出 |
|------|------|------|
| `00_data_preparation.py` | 土壤湿度格式转换、掩膜创建 | masks/combined_mask.tif |
| `01_phenology_extraction.py` | 从外部脚本读取物候指标 | GPP_phenology/, T_phenology/ |

### TRc累积蒸腾计算
| 脚本 | 物候源 | 输出 |
|------|--------|------|
| `02_TRc_calculation.py` | GPP物候 | TRc_annual/, Climatology/ |
| `02_TRc_calculation_T.py` ⭐ | T物候 | TRc_annual_T/, Climatology_T/ |

### TRc分解方法（三种）
| 脚本 | 方法 | 物候源 | 输出组分 |
|------|------|--------|---------|
| `03_decomposition_original.py` | Wang 2025原始 | GPP | TRpheno, TRproduct |
| `03_decomposition_original_T.py` ⭐ | Wang 2025原始 | T | TRpheno, TRproduct |
| `03_decomposition_timing_shape.py` ⭐⭐ | Timing/Shape新方法 | GPP | TRtiming, TRshape, TRsos, TRpos |
| `04_decomposition_improved.py` | 改进方法 | GPP | （开发中） |

### 统计分析（ΔSOS回归）
| 脚本 | 分解方法 | 物候源 | 响应变量 |
|------|---------|--------|---------|
| `05_statistical_analysis.py` | 原始分解 | GPP | TRc, TRpheno, TRproduct, Trate |
| `05_statistical_analysis_T.py` | 原始分解 | T | TRc, TRpheno, TRproduct, Trate |
| `05_statistical_analysis_timing_shape.py` ⭐⭐ | Timing/Shape | GPP | TRtiming, TRshape, TRsos, TRpos, Trate |

### 高级分析
| 脚本 | 功能 | 语言 |
|------|------|------|
| `06_SEM_analysis.R` | 结构方程模型 | R |
| `07_plotting_functions.py` | 可视化绘图 | Python |

### 工具脚本
| 脚本 | 功能 |
|------|------|
| `config.py` | 全局路径配置 |
| `verify_data.py` | 数据完整性验证 |
| `utils_climatology.py` | 气候态计算工具 |
| `utils_vegetation_stratification.py` | 植被分层工具 |

---

## 📂 数据输出目录

```
Wang2025_Analysis/
├── masks/                              掩膜文件
├── TRc_annual/                         TRc（GPP物候）
├── TRc_annual_T/ ⭐                     TRc（T物候）
├── Climatology/                        气候态（GPP物候）
├── Climatology_T/ ⭐                    气候态（T物候）
├── Decomposition/                      原始分解（GPP物候）
├── Decomposition_T/ ⭐                  原始分解（T物候）
├── Decomposition_TimingShape/ ⭐⭐       Timing/Shape分解
├── Statistical_Analysis/               统计分析（GPP-原始）
├── Statistical_Analysis_T/             统计分析（T-原始）
└── Statistical_Analysis_TimingShape/ ⭐⭐ 统计分析（Timing/Shape）
```

---

## 🎯 推荐阅读顺序

### 第一次使用
1. [项目说明](../README.md) - 从根目录README开始
2. [快速开始](00_README/START_HERE.md)
3. [数据配置](01_guides/数据配置说明.md)
4. [今日修改总结](02_changelogs/今日修改总结_20251224.md) - 了解最新功能

### 了解方法
1. [CHANGELOG](02_changelogs/CHANGELOG.md) - 完整版本历史
2. [技术修正记录](02_changelogs/技术修正记录.md) - 性能优化、符号修正、格式修复
3. [今日修改总结](02_changelogs/今日修改总结_20251224.md) - Timing/Shape新方法

### 结果分析
1. [检查统计结果差异](03_analysis/检查统计结果差异.md) - 诊断结果差异
2. [完成总结](04_summaries/完成总结.md) - 整体完成情况

---

## ⚠️ 重要提示

- ⭐ 标记表示2025-12-24新增/修改的文件
- ⭐⭐ 标记表示今日核心新增功能（Timing/Shape分解）
- 所有路径配置在 `config.py` 中统一管理
- 运行前请先阅读 [数据配置说明](01_guides/数据配置说明.md)

---

**文档生成**: 2025-12-24
**项目版本**: v2.0.0
