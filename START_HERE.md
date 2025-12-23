# 📍 从这里开始

## 🎯 您需要看哪个文档？

```
┌─────────────────────────────────────────────────────────┐
│  根据您的需求，选择对应的文档：                           │
└─────────────────────────────────────────────────────────┘

❓ 我想知道...                          👉 看这个文档

1️⃣  如何快速开始运行                    → CHECKLIST.md（检查清单）
2️⃣  所有修改内容和详细步骤              → 完成总结.md（主文档）
3️⃣  如何修改配置或切换数据源            → 数据配置说明.md
4️⃣  查看所有修改记录                    → CHANGELOG.md（修改记录）
```

---

## 🚀 最快上手方式（3步）

```bash
# 第1步：验证数据（2分钟）
python verify_data.py

# 第2步：准备掩膜（5分钟）
python 00_data_preparation.py

# 第3步：运行分析（10-15小时）
python 00_master_pipeline.py
```

**详细步骤请看**: [CHECKLIST.md](CHECKLIST.md)

---

## 📂 核心文件说明

### 🔧 配置文件
- **[config.py](config.py)** - 所有数据路径和参数配置
  - 已更新为您的数据路径（ERA5-Land TR、GPP物候等）
  - `USE_FOREST_MASK = False` - 先全局后植被分层

### 🐍 核心代码（按运行顺序）
1. **[00_data_preparation.py](00_data_preparation.py)** - 准备掩膜
2. ~~01_phenology_extraction.py~~ - 物候提取（您已有GPP物候，跳过）
3. **[02_TRc_calculation.py](02_TRc_calculation.py)** - 计算TRc（6-8小时）
4. **[03_decomposition_original.py](03_decomposition_original.py)** - 分解分析
5. **[05_statistical_analysis.py](05_statistical_analysis.py)** - 统计分析
6. **[07_plotting_functions.py](07_plotting_functions.py)** - 绘图

### 🛠️ 辅助工具
- **[verify_data.py](verify_data.py)** - 验证数据完整性（运行前必须用）
- **[utils_vegetation_stratification.py](utils_vegetation_stratification.py)** - 植被分层分析（分析后使用）
- **[00_master_pipeline.py](00_master_pipeline.py)** - 一键运行所有模块

### 📄 文档
- **[START_HERE.md](START_HERE.md)** - 本文档（导航页）
- **[CHECKLIST.md](CHECKLIST.md)** - 快速检查清单
- **[完成总结.md](完成总结.md)** - 详细操作指南（主文档）
- **[数据配置说明.md](数据配置说明.md)** - 配置详解和数据源切换
- **[文件清单与清理建议.md](文件清单与清理建议.md)** - 文件说明和清理建议

---

## ⚡ 当前配置摘要

```python
# 数据路径（已配置）
ROOT = "I:\F\Data4"
TR数据 = ERA5-Land格式（ERA5L_ET_transp_Daily_mm_YYYYMMDD.tif）
物候数据 = GPP物候（sos_gpp_{year}.tif）
土壤水分 = SMrz深层（推荐）

# 分析范围
时间: 1982-2018年（37年）
空间: ≥30°N 所有陆地（先全局后植被分层）

# 掩膜策略
USE_FOREST_MASK = False  ✅ 先计算全局，后分层
```

---

## 🎯 核心输出

运行完成后，结果保存在 `I:\F\Data4\Wang2025_Analysis\`：

```
Wang2025_Analysis/
├── masks/              # 掩膜文件
├── TRc_annual/         # 年度TRc（1982-2018）
├── Decomposition/      # 分解结果
├── Statistics/         # 统计分析
└── Figures/            # 生成的图表
```

---

## ❓ 常见问题快速解答

### Q: 我应该先看哪个文档？
**A**: 先看 [CHECKLIST.md](CHECKLIST.md)，5分钟了解完整流程

### Q: 如何修改数据路径？
**A**: 编辑 [config.py](config.py)，参考 [数据配置说明.md](数据配置说明.md)

### Q: 如何验证数据是否准备好？
**A**: 运行 `python verify_data.py`

### Q: 如何切换到T物候或表层土壤水分？
**A**: 查看 [数据配置说明.md](数据配置说明.md) 第4节

### Q: 为什么有这么多文档？
**A**: 确实有些重复，请看 [文件清单与清理建议.md](文件清单与清理建议.md)

---

## 📞 获取帮助

1. **数据验证失败**？运行 `python verify_data.py` 查看详细错误
2. **配置问题**？查看 [数据配置说明.md](数据配置说明.md)
3. **运行错误**？查看 [完成总结.md](完成总结.md) 的"故障排除"章节
4. **文档太多**？按照 [文件清单与清理建议.md](文件清单与清理建议.md) 清理

---

**开始分析**: 运行 `python verify_data.py` 🚀
