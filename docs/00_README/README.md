# Wang (2025) 复现代码 - 已适配您的数据

> 分析春季物候对森林蒸腾的直接和间接影响
> 基于 Wang et al. (2025) 论文方法
> **已更新**: 适配ERA5-Land TR数据、GPP物候、SMrz深层土壤水分

---

## 📍 快速导航

| 我想... | 查看文档 | 运行命令 |
|--------|---------|---------|
| 🚀 **立即开始** | [START_HERE.md](START_HERE.md) | `python verify_data.py` |
| ✅ **运行前检查** | [CHECKLIST.md](CHECKLIST.md) | - |
| 📖 **详细操作指南** | [完成总结.md](完成总结.md) | - |
| 🔧 **修改配置** | [数据配置说明.md](数据配置说明.md) | 编辑 `config.py` |
| 📝 **查看修改记录** | [CHANGELOG.md](CHANGELOG.md) | - |

---

## ⚡ 快速开始（3步）

```bash
# 1. 验证数据（2-5分钟）
python verify_data.py

# 2. 准备掩膜（5-10分钟）
python 00_data_preparation.py

# 3. 运行分析（10-15小时）
python 00_master_pipeline.py
```

**详细步骤**: 请查看 [CHECKLIST.md](CHECKLIST.md)

---

## 🎯 当前配置

### 数据路径（已配置好）
```python
根目录: I:\F\Data4
TR数据: ERA5-Land格式 (I:\F\Data4\Meteorological Data\ERA5_Land\...)
物候数据: GPP物候 (I:\F\Data4\Phenology_Output_1\GPP_phenology)
土壤水分: SMrz深层 (推荐用于森林蒸腾研究)
```

### 分析范围
- **时间**: 1982-2018年（37年）
- **空间**: ≥30°N 所有陆地
- **策略**: 先全局分析，后植被分层（`USE_FOREST_MASK = False`）

### 数据格式
| 数据 | 格式 | 示例 |
|------|------|------|
| TR | `ERA5L_ET_transp_Daily_mm_YYYYMMDD.tif` | ERA5L_ET_transp_Daily_mm_19820101.tif |
| 物候 | `sos_gpp_{year}.tif`, `pos_doy_gpp_{year}.tif` | sos_gpp_1982.tif |
| 土壤水分 | `SMrz_YYYYMMDD.tif` ⭐ 日尺度 | SMrz_19801218.tif |

---

## 📂 核心文件说明

### 🔧 配置文件
- **[config.py](config.py)** - 所有路径和参数配置（已更新）

### 🐍 代码文件（按运行顺序）

| 模块 | 文件 | 功能 | 运行时间 | 必需？ |
|------|------|------|---------|-------|
| 00 | [00_data_preparation.py](00_data_preparation.py) | 创建掩膜、检查数据 | 5-10分钟 | ✅ 必需 |
| 01 | [01_phenology_extraction.py](01_phenology_extraction.py) | 提取物候 | ~2小时 | ⚠️ 您已有GPP物候，跳过 |
| 02 | [02_TRc_calculation.py](02_TRc_calculation.py) | 计算累积蒸腾 | 6-8小时 | ✅ 必需 |
| 03 | [03_decomposition_original.py](03_decomposition_original.py) | 原版分解方法 | ~1小时 | ✅ 必需 |
| 04 | [04_decomposition_improved.py](04_decomposition_improved.py) | 改进分解方法 | ~2小时 | ⚠️ 可选 |
| 05 | [05_statistical_analysis.py](05_statistical_analysis.py) | 统计分析 | 2-3小时 | ✅ 必需 |
| 06 | [06_SEM_analysis.R](06_SEM_analysis.R) | SEM路径分析 | ~10分钟 | ⚠️ 需R环境 |
| 07 | [07_plotting_functions.py](07_plotting_functions.py) | 绘图 | ~20分钟 | ✅ 必需 |

### 🛠️ 辅助工具
- **[verify_data.py](verify_data.py)** - 数据验证脚本（运行前必须用）
- **[utils_vegetation_stratification.py](utils_vegetation_stratification.py)** - 植被分层分析
- **[utils_climatology.py](utils_climatology.py)** - 气候平均态计算（改进方法需要）
- **[00_master_pipeline.py](00_master_pipeline.py)** - 一键运行所有模块

### 📄 文档（已整理）
- **[START_HERE.md](START_HERE.md)** - 导航页（从这里开始）
- **[CHECKLIST.md](CHECKLIST.md)** - 快速检查清单
- **[完成总结.md](完成总结.md)** - 详细操作指南
- **[数据配置说明.md](数据配置说明.md)** - 配置详解
- **[CHANGELOG.md](CHANGELOG.md)** - 修改记录
- **[README.md](README.md)** - 本文档

---

## 📊 输出结果

运行完成后，结果保存在 `I:\F\Data4\Wang2025_Analysis\`：

```
Wang2025_Analysis/
├── masks/                    # 掩膜文件
│   ├── lat_mask.tif
│   ├── forest_mask.tif
│   └── combined_mask.tif
├── TRc_annual/               # 年度TRc（1982-2018）
│   └── TRc_*.tif
├── Decomposition/            # 分解结果
│   ├── TRc_av.tif           # 多年平均
│   ├── LSP_av.tif
│   ├── TRpheno_*.tif        # 物候分量
│   └── TRintensity_*.tif    # 强度分量
├── Statistics/               # 统计分析
│   ├── Trends/              # 趋势分析
│   ├── Attribution/         # 归因分析
│   └── Moving_Window/       # 滑动窗口
└── Figures/                  # 生成的图表
```

---

## 🔧 修改配置

### 切换数据源

**GPP物候 ⇄ T物候**:
```python
# 编辑 config.py
PHENO_DIR = ROOT / "Phenology_Output_1" / "GPP_phenology"  # 当前
# PHENO_DIR = ROOT / "Phenology_Output_1" / "T_phenology"  # 切换
```

**深层土壤水分 ⇄ 表层土壤水分**:
```python
# 编辑 config.py
SM_DAILY_DIR = GLEAM_ROOT / "SMrz" / "SMrz_Daily"  # 深层（推荐）
# SM_DAILY_DIR = GLEAM_ROOT / "SMs" / "SMs_Daily"  # 表层（备选）
```

**全局分析 ⇄ 仅森林**:
```python
# 编辑 config.py
USE_FOREST_MASK = False  # 全局分析（当前）
# USE_FOREST_MASK = True  # 仅森林
```

详细说明请查看 [数据配置说明.md](数据配置说明.md)

---

## ⚠️ 重要注意事项

### 1. 物候数据格式
确认您的物候文件命名为：
- `sos_gpp_YYYY.tif` (小写)
- `pos_doy_gpp_YYYY.tif` (小写，DOY值)
- `eos_gpp_YYYY.tif` (小写)

### 2. 土壤水分数据优势 ⭐
- 您的土壤水分是**日尺度** (`SMrz_YYYYMMDD.tif`)
- 与TR和GPP的时间分辨率完全匹配
- **无需插值**，可直接用于逐日精确分析
- 比月尺度数据更能捕捉短期土壤水分变化

### 3. 土地覆盖数据（可选）
- **当前配置** (`USE_FOREST_MASK = False`): **不需要**土地覆盖数据
- 仅在切换到仅森林模式时才需要 `MCD12Q1_IGBP_2018.tif`
- 植被分层分析时，`utils_vegetation_stratification.py` 会使用土地覆盖数据

### 4. 质量控制（可选）
建议使用物候代码输出的质量标记筛选低质量像元：
```python
quality_flags = read_geotiff(PHENO_DIR / f"quality_flags_gpp_{year}.tif")
high_quality_mask = (quality_flags == 7)  # 7项检查都通过
```

---

## 🐛 故障排除

| 问题 | 解决方案 |
|------|---------|
| "找不到物候数据" | 检查 `I:\F\Data4\Phenology_Output_1\GPP_phenology` 是否存在，运行 `python verify_data.py` |
| "找不到TR数据" | 确认文件格式为 `ERA5L_ET_transp_Daily_mm_YYYYMMDD.tif` |
| 内存不足 | 编辑 `config.py`: `BLOCK_SIZE=64`, `MAX_WORKERS=1` |
| UnicodeEncodeError | Windows控制台编码问题，不影响功能，可忽略 |
| 验证脚本报错 | 查看详细输出，检查缺失的文件和目录 |

详细排错请查看 [完成总结.md](完成总结.md) 的"故障排除"章节

---

## 📚 环境要求

### Python依赖
```bash
pip install -r requirements.txt
```

主要包：numpy, scipy, pandas, rasterio, matplotlib, cartopy, tqdm

### R依赖（仅SEM分析需要）
```R
install.packages(c("lavaan", "semPlot", "raster", "tidyverse"))
```

---

## 🔄 更新记录

查看详细修改记录: [CHANGELOG.md](CHANGELOG.md)

**最新更新** (2025-12-23):
- ✅ 适配ERA5-Land TR数据格式
- ✅ 适配GPP物候（小写格式）
- ✅ 配置SMrz深层土壤水分
- ✅ 添加 `USE_FOREST_MASK` 掩膜策略配置
- ✅ 创建数据验证脚本
- ✅ 整理文档结构

---

## 📞 获取帮助

1. **运行前**: 查看 [CHECKLIST.md](CHECKLIST.md)
2. **修改配置**: 查看 [数据配置说明.md](数据配置说明.md)
3. **遇到错误**: 查看 [完成总结.md](完成总结.md) 的故障排除章节
4. **了解修改**: 查看 [CHANGELOG.md](CHANGELOG.md)

---

## 📖 引用

如使用本代码，请引用原文：

```
Wang, X. et al. (2025). Direct and indirect effects of spring phenology
on forest transpiration. [期刊名称], [卷期页码].
```

---

**开始分析**: 请先运行 `python verify_data.py` 验证数据 🚀
