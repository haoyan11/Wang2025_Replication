# Wang (2025) 复现代码

> 分析春季物候对森林蒸腾的直接和间接影响
> 基于 Wang et al. (2025) 论文方法

**当前版本**: v1.2.2 (2025-01-10)

---

## 📚 文档导航

| 文档 | 说明 |
|------|------|
| **[使用说明.md](使用说明.md)** | 📖 完整使用指南（从这里开始）⭐ |
| **[docs/项目文档.md](docs/项目文档.md)** | 📋 项目文档（配置/修改记录/分析对比） |
| **[diagnostics/诊断工具说明.md](diagnostics/诊断工具说明.md)** | 🔧 诊断工具使用说明 |

---

## ⚡ 快速开始

```bash
# 1. 验证数据（2-5分钟）
python _verify_data.py

# 2. 准备掩膜（5-10分钟）
python 00_data_preparation.py

# 3. 运行分析（10-15小时）
python 00_master_pipeline.py
```

**详细步骤**: 请查看 [使用说明.md](使用说明.md)

---

## 🎯 核心功能

- ✅ **适配您的数据格式**: ERA5-Land TR、GPP物候（EPSG:4326）、SMrz深层土壤水分
- ✅ **固定窗口方法**: 剥离窗口选择效应，Fixed_Trate可为负值
- ✅ **完整Bug修复**: Bug Fix 3/4已修复，数据保留率100%
- ✅ **性能优化**: Numba JIT加速，统计分析10-50倍提速
- ✅ **诊断工具**: 自动检测数据问题

---

## 📊 最新更新 (v1.2.2)

### 重要修复
- **Bug Fix 4** (2025-01-10): 修复05c错误过滤Fixed_Trate负值（47.80%数据丢失 → 完整保留）
- **Bug Fix 3** (2025-01-05): 修复TR_fixed_window负值过滤（8574像元 → 26409像元）

### 新增功能
- 整合诊断工具模块（`diagnostics/`）
- 统一中文文档结构
- 完善错误检测与提示

详细修改记录: [docs/项目文档.md §4](docs/项目文档.md#4-修改记录)

---

## 🐛 故障排除

| 问题 | 解决方案 |
|------|---------|
| "找不到物候数据" | 运行 `python _verify_data.py` |
| "TR_fixed_window全部为正值" | 运行 `python diagnostics/check_decomposition.py` |
| 内存不足 | 编辑`_config.py`: `BLOCK_SIZE=64`, `MAX_WORKERS=1` |

更多帮助: [使用说明.md](使用说明.md) | [项目文档.md](docs/项目文档.md#7-结果诊断)

---

## 📖 引用

如使用本代码，请引用原文：

```
Wang, X. et al. (2025). Direct and indirect effects of spring phenology
on forest transpiration. [期刊名称], [卷期页码].
```

---

**开始使用**: 请先阅读 [使用说明.md](使用说明.md) 📖
