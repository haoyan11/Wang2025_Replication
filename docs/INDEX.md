# 📚 文档索引

**最后更新**: 2025-12-24

---

## 🚀 快速入门

### 新手必读
1. **[README.md](../README.md)** - 项目概述、快速开始、完整指南 ⭐ 从这里开始
2. **[CHECKLIST.md](01_guides/CHECKLIST.md)** - 快速检查清单（5分钟了解流程）
3. **[requirements.txt](../requirements.txt)** - Python环境依赖

---

## 📖 使用指南

### [01_guides/](01_guides/)
- **[数据配置说明.md](01_guides/数据配置说明.md)** - 数据路径配置、格式要求、准备步骤
- **[CHECKLIST.md](01_guides/CHECKLIST.md)** - 质量检查清单、验证步骤

---

## 📝 修改总结

### ⭐⭐ 推荐阅读
- **[项目修改总结汇编.md](项目修改总结汇编.md)** - **所有修改的完整汇总（35页）** ⭐ 2025-12-27更新
  - 项目完成总结
  - 今日修改总结 (2024-12-24)
  - 技术修正记录（3个重要修正）
  - 闰年DOY处理修复
  - TRpheno符号关系说明
  - **固定窗口方法与输出文件对比（第6章）** ⭐⭐ 新增

### [02_changelogs/](02_changelogs/)
- **[CHANGELOG.md](02_changelogs/CHANGELOG.md)** - 完整版本历史（v1.0.0 → v2.0.0）

---

## 📊 分析报告

### [03_analysis/](03_analysis/)
- **[检查统计结果差异.md](03_analysis/检查统计结果差异.md)** - 回归结果与Wang 2025论文对比诊断

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

### TRc分解方法（按方法分组）
| 脚本 | 方法 | 物候源 | 输出组分 |
|------|------|--------|---------|
| `03a_decomposition_wang2025.py` | 方法A: Wang 2025原始 | GPP | TRpheno, TRproduct |
| `03a_decomposition_wang2025_T.py` ⭐ | 方法A: Wang 2025原始 | T | TRpheno, TRproduct |
| `03b_decomposition_timing_shape.py` ⭐⭐ | 方法B: Timing/Shape新方法 | GPP | TRtiming, TRshape, TRsos, TRpos |

### 统计分析（对应各分解方法）
| 脚本 | 对应方法 | 物候源 | 响应变量 |
|------|---------|--------|---------|
| `04a_statistical_wang2025.py` | 方法A: Wang 2025 | GPP | TRc, TRpheno, TRproduct, Trate |
| `04a_statistical_wang2025_T.py` ⭐ | 方法A: Wang 2025 | T | TRc, TRpheno, TRproduct, Trate |
| `04b_statistical_timing_shape.py` ⭐⭐ | 方法B: Timing/Shape | GPP | TRtiming, TRshape, TRsos, TRpos, Trate |

### 高级分析
| 脚本 | 功能 | 语言 |
|------|------|------|
| `05_SEM_analysis.R` | 结构方程模型（SEM） | R |
| `06_plotting.py` | 可视化绘图 | Python |

### 工具脚本
| 脚本 | 功能 |
|------|------|
| `_config.py` | 全局路径配置 |
| `_verify_data.py` | 数据完整性验证 |
| `_utils_climatology.py` | 气候态计算工具 |
| `_utils_vegetation_stratification.py` | 植被分层工具 |

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
1. **[README.md](../README.md)** - 项目主页，从这里开始 ⭐
2. **[CHECKLIST.md](01_guides/CHECKLIST.md)** - 5分钟快速检查清单
3. **[数据配置说明.md](01_guides/数据配置说明.md)** - 数据准备和路径配置

### 了解所有修改 ⭐⭐
1. **[项目修改总结汇编](项目修改总结汇编.md)** - 一站式查看所有修改（推荐）
   - 包含完成总结、今日修改、技术修正、闰年处理、TRpheno符号说明、固定窗口方法
   - 35页完整汇总（2025-12-27更新）

### 查看详细历史
1. [CHANGELOG](02_changelogs/CHANGELOG.md) - 按版本组织的完整变更历史

### 结果诊断
1. [检查统计结果差异](03_analysis/检查统计结果差异.md) - 与Wang 2025论文对比

---

## ⚠️ 重要提示

- ⭐ 标记表示2025-12-24新增/修改的文件
- ⭐⭐ 标记表示今日核心新增功能（Timing/Shape分解）
- 所有路径配置在 `config.py` 中统一管理
- 运行前请先阅读 [数据配置说明](01_guides/数据配置说明.md)

---

**文档生成**: 2025-12-24
**最后更新**: 2025-12-27 (固定窗口方法章节合并)
**项目版本**: v2.2.0 (固定窗口方法集成版)
