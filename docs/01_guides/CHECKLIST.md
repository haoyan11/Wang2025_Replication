# Wang2025 分析流程检查清单

## 📋 快速检查清单

### 准备阶段
- [ ] 已阅读 [README.md](../../README.md) 项目说明
- [ ] 已安装Python依赖: `pip install -r requirements.txt`
- [ ] 确认 `_config.py` 配置正确（特别是 `USE_FOREST_MASK = False`）

### 数据验证（必须）
- [ ] 运行 `python _verify_data.py`
- [ ] 所有检查项通过（10/10）
- [ ] 特别确认：
  - [ ] 物候数据: 37/37年份完整
  - [ ] TR数据: ERA5-Land格式抽样通过
  - [ ] GPP数据: 日尺度抽样通过
  - [ ] 土壤水分: SMrz深层抽样通过

### 掩膜准备（必须）
- [ ] 运行 `python 00_data_preparation.py`
- [ ] 看到提示: "掩膜策略: 所有陆地（≥30°N）"
- [ ] 生成3个掩膜文件:
  - [ ] `masks/lat_mask.tif`
  - [ ] `masks/forest_mask.tif`
  - [ ] `masks/combined_mask.tif`

### 核心分析（耗时10-15小时）
- [ ] 运行 `python 00_master_pipeline.py`
- [ ] 或分步运行:
  - [ ] `python 02_TRc_calculation.py` (6-8小时)
  - [ ] `python 03_decomposition_original.py` (1小时)
  - [ ] `python 04_decomposition_improved.py` (1小时)
  - [ ] `python 05_statistical_analysis.py` (2小时)
  - [ ] `Rscript 06_SEM_analysis.R` (可选，需R环境)
  - [ ] `python 07_plotting_functions.py`

### 结果检查
- [ ] 检查输出目录 `I:\F\Data4\Wang2025_Analysis\`:
  - [ ] `TRc_annual/` 包含37个年份文件
  - [ ] `Decomposition/` 包含分解结果
  - [ ] `Statistics/` 包含趋势分析
  - [ ] `Figures/` 包含生成的图表

### 植被分层（可选）
- [ ] 运行 `python _utils_vegetation_stratification.py`
- [ ] 检查输出:
  - [ ] `Vegetation_Stratification/TRc_av_by_vegetation.csv`
  - [ ] `Vegetation_Stratification/vegetation_masks/` 包含各植被类型掩膜

---

## 🚨 常见问题快速解决

| 问题 | 快速解决 |
|------|---------|
| "找不到物候数据" | 检查 `I:\F\Data4\Phenology_Output_1\GPP_phenology` 是否存在，文件名是否为小写 `sos_gpp_*.tif` |
| "找不到TR数据" | 检查 TR目录路径，文件名格式是否为 `ERA5L_ET_transp_Daily_mm_YYYYMMDD.tif` |
| "内存不足" | 编辑 `_config.py`: `BLOCK_SIZE=64`, `MAX_WORKERS=1` |
| UnicodeEncodeError | 不影响功能，可忽略 |
| 验证脚本报错 | 查看详细输出，检查缺失的文件和目录 |

---

## 📞 获取帮助

- 详细配置: [数据配置说明.md](数据配置说明.md)
- 所有修改: [项目修改总结汇编.md](../项目修改总结汇编.md) ⭐ (25页完整汇总)
- 版本历史: [CHANGELOG.md](../02_changelogs/CHANGELOG.md)
- 文档导航: [INDEX.md](../INDEX.md)

---

## ⏱️ 预估时间

| 步骤 | 预估时间 |
|------|---------|
| 数据验证 | 2-5分钟 |
| 准备掩膜 | 5-10分钟 |
| TRc计算 | 6-8小时 |
| 分解分析 | 2-3小时 |
| 统计分析 | 2-3小时 |
| 绘图 | 30分钟 |
| **总计** | **10-15小时** |
| 植被分层（可选） | 15-45分钟 |

---

**当前日期**: 2025-12-23
**代码版本**: 已更新以匹配用户数据路径和格式
**掩膜模式**: 全局分析（≥30°N），支持后续植被分层
