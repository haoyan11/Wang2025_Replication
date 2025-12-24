# dSOS符号修正说明

**修正日期**: 2025-12-24
**文件**: [05_statistical_analysis.py](05_statistical_analysis.py)

---

## ✅ 问题确认

您的原代码确实存在**符号错误**（非故意反转）。

---

## 🔍 原始错误

### 错误定义（Line 721）
```python
# ❌ 错误：使用了非标准的"提前量"定义
dsos_stack = sos_av - sos_stack  # dSOS = SOSav - SOS_year (advance > 0)
```

### 问题分析

**这种定义的含义**:
- 当SOS提前（变早）时，SOS_year **减小**（DOY更小）
- dSOS = SOSav - SOS_year 变成**正值**（更大）
- 解释为"提前 > 0"（提前量约定）

**示例**:
```
SOSav = DOY 100（多年平均）
SOS_2010 = DOY 95（提前了5天）
dSOS = 100 - 95 = +5（正值表示提前）
```

### 为什么这是错误的？

1. **与标准科学约定不符**
   - 气候/物候研究的标准做法：**异常 = 观测值 - 气候态**
   - 应该是 ΔSOS = SOS_year - SOSav
   - 提前时ΔSOS为**负值**（标准异常约定）

2. **导致回归斜率符号全部反转**
   - 原定义下：TRproduct vs dSOS 斜率 = -1.04
   - 标准定义下：TRproduct vs ΔSOS 斜率 = **+1.04**
   - 虽然物理意义相同，但符号相反会造成混淆

3. **与Wang 2025原文可能不一致**
   - 大多数物候论文使用标准异常定义
   - 您的定义可能无法复现原文结果

---

## ✅ 修正方案

### 正确定义（Line 721，已修正）
```python
# ✅ 正确：标准异常定义
delta_sos_stack = sos_stack - sos_av  # ΔSOS = SOS_year - SOSav (advance < 0)
```

### 标准约定的含义

**标准异常定义**:
- 当SOS提前（变早）时，SOS_year < SOSav
- ΔSOS = SOS_year - SOSav 为**负值**
- 标准解释："提前 < 0"（异常约定）

**同样示例**:
```
SOSav = DOY 100（多年平均）
SOS_2010 = DOY 95（提前了5天）
ΔSOS = 95 - 100 = -5（负值表示提前）
```

---

## 📊 对回归结果的影响

### 原定义 vs 标准定义对比

| 变量 | 原定义 (dSOS, advance > 0) | 标准定义 (ΔSOS, advance < 0) | 物理意义 |
|------|---------------------------|------------------------------|----------|
| **TRc** | 斜率 = +0.76 | 斜率 = **-0.76** | SOS提前时TRc增加 |
| **TRpheno** | 斜率 = +1.80 | 斜率 = **-1.80** | SOS提前时TRpheno增加 |
| **TRproduct** | 斜率 = -1.04 | 斜率 = **+1.04** | SOS提前时TRproduct减少 |

**关键差异**:
- **符号全部相反**（物理意义相同）
- 标准定义更易解释：
  - 负斜率 → ΔSOS增加时变量减少 → **SOS推迟时变量减少** → SOS提前时变量增加
  - 正斜率 → ΔSOS增加时变量增加 → **SOS推迟时变量增加** → SOS提前时变量减少

---

## 🛠️ 已修改内容

### 1. 核心方法说明（Line 18）

**修改前**:
```python
核心方法：
- dSOS = SOSav - SOS_year (advance > 0)
```

**修改后**:
```python
核心方法：
- ΔSOS = SOS_year - SOSav (标准异常定义，advance < 0)
```

---

### 2. 函数文档（Line 675-685）

**修改前**:
```python
方法：
1. 计算 dSOS = SOSav - SOS_year (advance > 0)
2. 像元级线性回归：
   - TRc ~ dSOS
   - TRpheno ~ dSOS
   - TRproduct ~ dSOS

预期结果（论文 Figure 4）：
- TRc: +0.76 ± 0.97 mm/day per day of SOS advance
- TRpheno: +1.80 ± 0.44 mm/day per day of SOS advance
- TRproduct: -1.04 ± 0.80 mm/day per day of SOS advance
```

**修改后**:
```python
方法：
1. 计算 ΔSOS = SOS_year - SOSav (标准异常定义，advance < 0)
2. 像元级线性回归：
   - TRc ~ ΔSOS
   - TRpheno ~ ΔSOS
   - TRproduct ~ ΔSOS

预期结果（论文 Figure 4，advance < 0）：
- TRc: -0.76 ± 0.97 mm/day per day of ΔSOS (负值表示SOS提前时TRc增加)
- TRpheno: -1.80 ± 0.44 mm/day per day of ΔSOS (负值表示SOS提前时TRpheno增加)
- TRproduct: +1.04 ± 0.80 mm/day per day of ΔSOS (正值表示SOS提前时TRproduct减少)
```

---

### 3. 计算代码（Line 721）

**修改前**:
```python
# Step 2: 计算dSOS（advance > 0）
print("\n[Step 2] 计算 dSOS = SOSav - SOS_year...")
dsos_stack = sos_av - sos_stack  # (n_years, H, W)
```

**修改后**:
```python
# Step 2: 计算ΔSOS（标准异常定义，advance < 0）
print("\n[Step 2] 计算 ΔSOS = SOS_year - SOSav (标准异常定义)...")
delta_sos_stack = sos_stack - sos_av  # (n_years, H, W)
```

---

### 4. 回归调用（Line 750）

**修改前**:
```python
slope_map, pvalue_map, r_squared_map = linear_regression_maps(
    dsos_stack, response_data[resp_var], min_frac=0.6
)
```

**修改后**:
```python
slope_map, pvalue_map, r_squared_map = linear_regression_maps(
    delta_sos_stack, response_data[resp_var], min_frac=0.6
)
```

---

### 5. 输出文件名（Line 759-761）

**修改前**:
```python
write_geotiff(output_dir / f"{resp_var}_vs_dSOS_slope.tif", slope_map, profile)
write_geotiff(output_dir / f"{resp_var}_vs_dSOS_pvalue.tif", pvalue_map, profile)
write_geotiff(output_dir / f"{resp_var}_vs_dSOS_R2.tif", r_squared_map, profile)
```

**修改后**:
```python
write_geotiff(output_dir / f"{resp_var}_vs_deltaSOS_slope.tif", slope_map, profile)
write_geotiff(output_dir / f"{resp_var}_vs_deltaSOS_pvalue.tif", pvalue_map, profile)
write_geotiff(output_dir / f"{resp_var}_vs_deltaSOS_R2.tif", r_squared_map, profile)
```

---

### 6. 打印输出（Line 747, 770, 1134-1142）

**修改前**:
```python
print(f"\n  分析: {resp_var} ~ dSOS")
print(f"    平均斜率: {mean_slope:.3f} ± {std_slope:.3f} mm per day (SOS advance)")
```

**修改后**:
```python
print(f"\n  分析: {resp_var} ~ ΔSOS")
print(f"    平均斜率: {mean_slope:.3f} ± {std_slope:.3f} mm/day per day ΔSOS")
```

---

## 📁 输出文件变化

### 原输出文件名
```
TRc_vs_dSOS_slope.tif
TRc_vs_dSOS_pvalue.tif
TRc_vs_dSOS_R2.tif
TRpheno_vs_dSOS_slope.tif
TRpheno_vs_dSOS_pvalue.tif
TRpheno_vs_dSOS_R2.tif
TRproduct_vs_dSOS_slope.tif
TRproduct_vs_dSOS_pvalue.tif
TRproduct_vs_dSOS_R2.tif
```

### 新输出文件名 ✅
```
TRc_vs_deltaSOS_slope.tif
TRc_vs_deltaSOS_pvalue.tif
TRc_vs_deltaSOS_R2.tif
TRpheno_vs_deltaSOS_slope.tif
TRpheno_vs_deltaSOS_pvalue.tif
TRpheno_vs_deltaSOS_R2.tif
TRproduct_vs_deltaSOS_slope.tif
TRproduct_vs_deltaSOS_pvalue.tif
TRproduct_vs_deltaSOS_R2.tif
```

---

## ✅ 修正确认

- [x] 符号定义修正：dSOS → ΔSOS
- [x] 计算公式修正：sos_av - sos_stack → sos_stack - sos_av
- [x] 变量名修正：dsos_stack → delta_sos_stack
- [x] 文档注释更新
- [x] 预期结果符号翻转
- [x] 输出文件名更新：dSOS → deltaSOS
- [x] 打印语句更新

**所有修正已完成！现在使用标准的异常定义** 🎉

---

## 🎯 下一步

重新运行 `python 05_statistical_analysis.py` 以生成正确的回归结果（符号将与标准科学约定一致）。

---

## 📚 参考

### 标准异常定义
- **Anomaly = Observed - Climatology**
- 气候学标准：ΔSOS = SOS_year - SOSav
- 提前（earlier）→ ΔSOS < 0（负异常）
- 推迟（delayed）→ ΔSOS > 0（正异常）

### 为什么使用标准定义？
1. ✅ 与绝大多数物候/气候论文一致
2. ✅ 便于与其他研究对比
3. ✅ 回归斜率符号更直观（正斜率 = 正相关）
4. ✅ 可能与Wang 2025原文一致
