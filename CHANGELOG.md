# 修改记录 (Changelog)

本文件记录对Wang2025复现代码的所有修改。

---

## [1.2.0] - 2025-12-23

### 🎯 主要更新：土壤水分数据格式修正（月尺度→日尺度）

**修改者**: Claude (应用户数据格式)
**修改原因**: 用户的SMrz和SMs数据是日尺度（SMrz_YYYYMMDD.tif），而非之前配置的月尺度

#### ✅ 代码修改

##### 1. **config.py**
**修改位置**: 第93-96行, 第185-207行

**修改内容**:
- 更新SM文件格式配置（第93-96行）
  ```python
  # 修改前:
  SM_FILE_FORMAT = "SMrz_{year}{month:02d}.tif"  # 月尺度

  # 修改后:
  SM_FILE_FORMAT = "SMrz_{date}.tif"  # 日尺度，{date} = YYYYMMDD
  ```

- 新增 `get_SM_file_path()` 函数（第185-207行）
  ```python
  def get_SM_file_path(date_obj):
      """根据日期获取土壤水分文件路径（日尺度）"""
      yyyymmdd = date_obj.strftime("%Y%m%d")
      file_path = SM_DAILY_DIR / f"SMrz_{yyyymmdd}.tif"
      if file_path.exists():
          return file_path
      return None
  ```

**原因**: 用户数据格式为 `SMrz_19801218.tif`（日尺度），需要相应修改配置和辅助函数

---

##### 2. **verify_data.py**
**修改位置**: 第11-13行, 第117-143行, 第160-166行

**修改内容**:
- 导入 `USE_FOREST_MASK` 配置（第13行）
  ```python
  from config import (
      ROOT, TR_DAILY_DIR, PHENO_DIR, GPP_DAILY_DIR, SM_DAILY_DIR,
      LANDCOVER_FILE, YEAR_START, YEAR_END, USE_FOREST_MASK, get_TR_file_path
  )
  ```

- 更新 `check_SM_data()` 函数为日尺度（第117-143行）
  ```python
  # 修改前: 检查月份
  test_months = [(1982, 1), (2000, 6), (2018, 12)]
  sm_file = SM_DAILY_DIR / f"SMrz_{year}{month:02d}.tif"

  # 修改后: 检查日期
  test_dates = [datetime(1982, 1, 15), datetime(2000, 6, 15), datetime(2018, 12, 15)]
  sm_file = SM_DAILY_DIR / f"SMrz_{test_date.strftime('%Y%m%d')}.tif"
  ```

- 土地覆盖数据检查改为可选（第160-166行）
  ```python
  # 仅在 USE_FOREST_MASK=True 时检查土地覆盖数据
  if USE_FOREST_MASK:
      all_checks.append(check_file(LANDCOVER_FILE, "MODIS IGBP"))
  else:
      print("    ⚠ 跳过检查（当前为全局分析模式，仅在仅森林模式时需要）")
  ```

**原因**:
1. 适配日尺度SM数据格式
2. 土地覆盖数据在当前配置（USE_FOREST_MASK=False）下不需要，避免误报错误

---

#### 📝 文档更新

##### 3. **README.md**
**修改位置**: 第58行, 第159-168行

**修改内容**:
- 更新数据格式表中SM格式（第58行）
  ```markdown
  土壤水分 | `SMrz_YYYYMMDD.tif` ⭐ 日尺度 | SMrz_19801218.tif
  ```
- 更新注意事项（第159-168行）
  - 强调SM日尺度优势：与TR、GPP匹配，无需插值
  - 说明土地覆盖数据在当前配置下不需要

**原因**: 让用户一眼看到数据格式和配置要求

---

##### 4. **数据配置说明.md**
**修改位置**: 第27-43行

**修改内容**:
- 更新SM数据格式说明
  ```markdown
  # 修改前:
  文件格式: SMrz_YYYYMM.tif
  示例: SMrz_198201.tif (1982年1月)

  # 修改后:
  文件格式: SMrz_YYYYMMDD.tif ⭐ **日尺度**
  示例: SMrz_19820101.tif (1982年1月1日)
  优势: 日尺度数据比月尺度更精确，可以捕捉短期土壤水分变化
  ```

**原因**: 更新文档以匹配实际数据格式，强调日尺度的优势

---

##### 5. **土壤水分格式修正说明.md** (新建)
**功能**: 本次修改的总结文档

**内容**:
- 修改的文件清单
- 数据格式对比（修改前 vs 修改后）
- 日尺度数据优势说明
- 土地覆盖数据影响分析
- 测试方法和预期输出
- 常见问题解答

**原因**: 为用户提供清晰的修改总结，方便查看本次更新内容

**注意**: 看完后可删除（内容已整合到CHANGELOG）

---

### 📋 关于土地覆盖数据

**问题**: 用户报告土地覆盖文件 `MCD12Q1_IGBP_2018.tif` 不存在

**影响分析**: ✅ **不影响当前分析**

**原因**:
1. 当前配置 `USE_FOREST_MASK = False`（全局分析模式）
2. 土地覆盖数据**仅在** `USE_FOREST_MASK = True` 时才需要
3. 用户采用"先全局分析，后植被分层"的工作流程
4. 植被分层时会使用 `utils_vegetation_stratification.py`，那时才需要土地覆盖数据

**解决方案**: `verify_data.py` 现在会根据 `USE_FOREST_MASK` 配置智能跳过不需要的检查

---

### 🎯 数据格式总结（更新后）

| 数据类型 | 时间分辨率 | 文件格式 | 示例 |
|---------|-----------|---------|------|
| TR（蒸腾） | 日尺度 | `ERA5L_ET_transp_Daily_mm_YYYYMMDD.tif` | ERA5L_ET_transp_Daily_mm_19820101.tif |
| **SM（土壤水分）** | **日尺度** | **`SMrz_YYYYMMDD.tif`** | **SMrz_19801218.tif** ⭐ 已更新 |
| GPP | 日尺度 | `GPP_YYYYMMDD.tif` | GPP_19820101.tif |
| 物候 | 年尺度 | `sos_gpp_{year}.tif` | sos_gpp_1982.tif |

---

## [1.1.0] - 2025-12-23

### 🎯 主要更新：文档整理和修改记录系统

**修改者**: Claude (应用户要求)
**修改原因**: 解决文档冗余问题（重复率>60%），建立修改记录系统

#### ✅ 新增
- **README.md** - 全新的主说明文档，清晰简洁（7.5KB）
- **CHANGELOG.md** - 本文件，记录所有修改（11KB）
- **START_HERE.md** - 导航页，引导用户查看正确文档（4.5KB）
- **整理完成.md** - 本次整理的总结报告（临时文档，看完可删除）

#### 📝 文档整理
- **归档** `archive/README.md.original` - 原始说明（已过时）
- **归档** `archive/QUICKSTART.md.original` - 原始快速指南（已过时）
- **删除** `修改说明.md` - 内容已整合到"完成总结.md"
- **删除** `WORKFLOW_NO_FOREST_MASK.md` - 内容已整合到"完成总结.md"
- **删除** `文件清单与清理建议.md` - 临时分析文档

#### 🔧 文档结构（整理后）
```
核心文档（5个）:
├── README.md                   ← 主说明文档（新建）
├── START_HERE.md               ← 导航页（新建）
├── CHECKLIST.md                ← 快速检查清单（保留）
├── 完成总结.md                 ← 详细操作指南（保留）
├── 数据配置说明.md             ← 配置详解（保留）
└── CHANGELOG.md                ← 修改记录（新建）

归档文档:
└── archive/
    ├── README.md.original
    └── QUICKSTART.md.original
```

---

## [1.0.0] - 2025-12-23

### 🎯 主要更新：适配用户数据路径和格式

#### ✅ 核心代码修改

##### 1. **config.py**
**修改位置**: 第13-44行, 第55-65行, 第74-105行, 第137-159行

**修改内容**:
- 更新所有数据路径为用户实际路径
  ```python
  ROOT = Path(r"I:\F\Data4")
  TR_DAILY_DIR = ROOT / "Meteorological Data" / "ERA5_Land" / "ET_components" / "ET_transp" / "ET_transp_Daily" / "ET_transp_Daily_2"
  PHENO_DIR = ROOT / "Phenology_Output_1" / "GPP_phenology"
  GPP_DAILY_DIR = ROOT / "GLASS_GPP" / "GLASS_GPP_daily_interpolated"
  SM_DAILY_DIR = GLEAM_ROOT / "SMrz" / "SMrz_Daily"  # 深层土壤水分
  ```
- 添加掩膜策略配置
  ```python
  USE_FOREST_MASK = False  # 先全局分析，后植被分层
  ```
- 更新文件命名格式
  ```python
  TR_FILE_FORMAT = "ERA5L_ET_transp_Daily_mm_{date}.tif"  # ERA5-Land格式
  PHENO_FILE_FORMAT = {
      'SOS': 'sos_gpp_{year}.tif',      # 小写格式
      'POS': 'pos_doy_gpp_{year}.tif',
      'EOS': 'eos_gpp_{year}.tif',
  }
  ```
- 更新 `get_TR_file_path()` 函数
  ```python
  def get_TR_file_path(date_obj):
      yyyymmdd = date_obj.strftime("%Y%m%d")
      file_path = TR_DAILY_DIR / f"ERA5L_ET_transp_Daily_mm_{yyyymmdd}.tif"
      if file_path.exists():
          return file_path
      return None
  ```

**原因**: 适配用户的ERA5-Land TR数据、GPP物候数据、SMrz深层土壤水分

---

##### 2. **00_data_preparation.py**
**修改位置**: 第22-24行, 第229-263行, 第298-358行

**修改内容**:
- 导入 `USE_FOREST_MASK` 配置
  ```python
  from config import (
      ROOT, LAT_MIN, FOREST_CLASSES, NODATA_OUT,
      USE_FOREST_MASK,  # ← 新增
      LANDCOVER_FILE, TR_DAILY_DIR, PHENO_DIR
  )
  ```
- 修复f-string警告（第246-248行）
  ```python
  # 修复前（有警告）:
  print(f"缺失: {['SOS', 'POS', 'EOS'][i] for i, x in enumerate(files_exist) if not x}")

  # 修复后:
  missing_types = [name for i, name in enumerate(['SOS', 'POS', 'EOS']) if not files_exist[i]]
  print(f"    {year}: ✗ (缺失: {missing_types})")
  ```
- 更新物候文件检查为小写格式（第232-234行）
  ```python
  sos_file = PHENO_DIR / f"sos_gpp_{year}.tif"  # 物候代码输出小写
  pos_file = PHENO_DIR / f"pos_doy_gpp_{year}.tif"
  eos_file = PHENO_DIR / f"eos_gpp_{year}.tif"
  ```
- 更新TR文件检查为ERA5-Land格式（第256行）
  ```python
  tr_file = TR_DAILY_DIR / f"ERA5L_ET_transp_Daily_mm_{test_date.strftime('%Y%m%d')}.tif"
  ```
- 添加掩膜模式逻辑（第316-341行）
  ```python
  if USE_FOREST_MASK:
      # 模式A: 森林掩膜
      forest_mask, _ = create_forest_mask(LANDCOVER_FILE, FOREST_CLASSES)
      combined_mask = combine_masks(lat_mask, forest_mask)
  else:
      # 模式B: 仅纬度掩膜
      forest_mask = np.ones_like(lat_mask)  # 占位文件
      combined_mask = lat_mask.copy()
  ```

**原因**:
1. 修复代码警告
2. 适配GPP物候的小写命名格式
3. 适配ERA5-Land TR数据格式
4. 支持两种掩膜策略（全局分析 vs 仅森林）

---

##### 3. **02_TRc_calculation.py**
**修改位置**: 第33-46行, 第73-100行, 第167-189行

**修改内容**:
- 简化 `get_TR_file_path()` 函数（第33-46行）
  ```python
  def get_TR_file_path(date_obj):
      """
      获取ERA5-Land TR日数据文件路径
      格式: ERA5L_ET_transp_Daily_mm_YYYYMMDD.tif
      """
      yyyymmdd = date_obj.strftime("%Y%m%d")
      file_path = TR_DAILY_DIR / f"ERA5L_ET_transp_Daily_mm_{yyyymmdd}.tif"
      if file_path.exists():
          return file_path
      return None
  ```
- 更新物候文件读取为小写格式（第90-91行）
  ```python
  sos_file = PHENO_DIR / f"sos_gpp_{year}.tif"
  pos_file = PHENO_DIR / f"pos_doy_gpp_{year}.tif"  # 注意是pos_doy
  ```
- 块处理版本同样更新（第182-183行）

**原因**: 适配ERA5-Land TR格式和GPP物候小写命名

---

##### 4. **03_decomposition_original.py**
**修改位置**: 第20-24行, 第76-83行, 第119-124行

**修改内容**:
- 更改物候目录为GPP物候（第22行）
  ```python
  PHENO_DIR = ROOT / "Phenology_Output_1" / "GPP_phenology"  # GPP物候
  ```
- 更新所有物候文件读取为小写格式
  ```python
  # 在 calculate_multiyear_mean 函数（第80-81行）
  sos, _ = read_geotiff(PHENO_DIR / f"sos_gpp_{year}.tif")
  pos, profile = read_geotiff(PHENO_DIR / f"pos_doy_gpp_{year}.tif")

  # 在 decompose_TRc_original 函数（第122-123行）
  sos_y, _ = read_geotiff(PHENO_DIR / f"sos_gpp_{year}.tif")
  pos_y, _ = read_geotiff(PHENO_DIR / f"pos_doy_gpp_{year}.tif")
  ```

**原因**: 使用GPP物候数据，适配小写命名格式

---

#### ✅ 新增工具

##### 5. **verify_data.py** (新建)
**功能**: 数据验证脚本

**主要功能**:
1. 检查所有必需目录是否存在
2. 验证物候数据完整性（1982-2018全部年份）
3. 抽样检查TR、GPP、土壤水分数据
4. 生成验证报告

**使用方法**:
```bash
python verify_data.py
```

**原因**: 在运行分析前快速验证数据是否准备就绪，避免运行到一半才发现数据缺失

---

#### 📄 文档创建

##### 6. **数据配置说明.md** (新建, 13KB)
**内容**:
- 所有数据路径和格式说明
- SMrz vs SMs土壤水分对比（推荐使用SMrz深层）
- 从物候代码借鉴的6个关键特性
- 如何切换数据源（GPP/T物候，深层/表层土壤水分）
- 数据验证方法
- 常见问题解答

##### 7. **完成总结.md** (新建, 14KB)
**内容**:
- 修改的文件清单和具体位置
- 下一步操作指南（3步）
- 如何切换数据源
- 植被分层分析方法
- 故障排除

##### 8. **CHECKLIST.md** (新建, 3KB)
**内容**:
- 快速检查清单
- 常见问题快速解决
- 预估时间

---

### 🎯 配置建议

#### 土壤水分选择：SMrz（深层）推荐

| 特性 | SMrz (深层) | SMs (表层) |
|------|------------|-----------|
| 与蒸腾关系 | ⭐⭐⭐⭐⭐ 直接相关 | ⭐⭐ 间接相关 |
| 森林根系 | ⭐⭐⭐⭐⭐ 深根依赖深层水 | ⭐⭐ 浅层为主 |
| 文献标准 | ⭐⭐⭐⭐⭐ 标准做法 | ⭐⭐ 较少使用 |
| Wang论文对应 | ⭐⭐⭐⭐⭐ 完全对应 | ⭐ 不对应 |

**配置**: `config.py` 已设置为 `SM_DAILY_DIR = GLEAM_ROOT / "SMrz" / "SMrz_Daily"`

---

### 🔄 掩膜策略

**当前配置**: `USE_FOREST_MASK = False`

**含义**:
- ✅ 先计算≥30°N所有陆地的结果
- ✅ 完成后运行 `utils_vegetation_stratification.py` 可提取各植被类型统计
- ✅ 优点：一次计算，多次使用；可灵活对比不同植被类型

**如需切换**:
```python
# 编辑 config.py
USE_FOREST_MASK = True  # 改为True，仅计算森林
```

---

### ⚠️ 已知注意事项

1. **物候文件命名**: 必须是小写（`sos_gpp_`, `pos_doy_gpp_`, `eos_gpp_`）
2. **土壤水分时间分辨率**: 月尺度（`SMrz_YYYYMM.tif`），TR是日尺度
3. **Module 01可跳过**: 您已有GPP物候数据，无需运行物候提取

---

## 📋 文件修改总结

### 修改的文件（4个核心代码）
1. ✅ `config.py` - 更新所有数据路径和格式
2. ✅ `00_data_preparation.py` - 修复警告、支持掩膜策略
3. ✅ `02_TRc_calculation.py` - 适配ERA5-Land格式
4. ✅ `03_decomposition_original.py` - 使用GPP物候

### 新增的文件（5个工具+6个文档）
**工具**:
1. ✅ `verify_data.py` - 数据验证脚本

**文档**:
1. ✅ `README.md` - 主说明文档（重建）
2. ✅ `START_HERE.md` - 导航页
3. ✅ `CHANGELOG.md` - 本文件
4. ✅ `数据配置说明.md` - 详细配置指南
5. ✅ `完成总结.md` - 完成报告
6. ✅ `CHECKLIST.md` - 检查清单

### 归档的文件（2个）
1. ✅ `archive/README.md.original` - 原始说明（已过时）
2. ✅ `archive/QUICKSTART.md.original` - 原始快速指南（已过时）

### 未修改的文件
- ✅ `04_decomposition_improved.py` - 改进分解方法（无需修改）
- ✅ `05_statistical_analysis.py` - 统计分析（无需修改）
- ✅ `06_SEM_analysis.R` - SEM分析（无需修改）
- ✅ `07_plotting_functions.py` - 绘图（无需修改）
- ✅ `utils_climatology.py` - 气候平均态（无需修改）
- ✅ `utils_vegetation_stratification.py` - 植被分层（无需修改）
- ✅ `00_master_pipeline.py` - 主流程（无需修改）
- ✅ `requirements.txt` - 依赖包（无需修改）

---

## 🔮 未来修改规范

### 规则1: 优先修改现有文件
- ✅ **优先**: 在现有文件基础上修改
- ❌ **避免**: 轻易创建新文件导致冗余

### 规则2: 每次修改必须记录
修改任何文件后，必须在本文件（CHANGELOG.md）中记录：
1. **修改的文件名**
2. **修改位置**（行号范围）
3. **修改内容**（代码片段）
4. **修改原因**

### 规则3: 记录格式
```markdown
#### 文件名
**修改位置**: 第X-Y行

**修改内容**:
​```python
# 代码片段
​```

**原因**: 简要说明
```

### 规则4: 版本号规则
- **主版本号** (X.0.0): 重大架构变更
- **次版本号** (1.X.0): 功能新增或重要修改
- **修订号** (1.0.X): Bug修复、文档更新

---

## 📞 查看修改

- **查看特定文件的修改**: 在本文件中搜索文件名
- **查看特定日期的修改**: 查看对应日期的章节
- **查看最新修改**: 查看文件顶部的最新版本

---

**最后更新**: 2025-12-23
**当前版本**: 1.1.0
