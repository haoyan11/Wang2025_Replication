#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module 02: TRc（累积蒸腾）计算
计算POS-EOS窗口内的累积蒸腾量
TRc_y = Σ[POS to EOS] TR(t)
"""

import numpy as np
import rasterio
from pathlib import Path
import os
from tqdm import tqdm
from datetime import datetime, timedelta
from rasterio.windows import Window
from concurrent.futures import ProcessPoolExecutor, as_completed
import warnings
warnings.filterwarnings('ignore')

# 导入配置
from _config import (
    PHENO_DIR, TR_DAILY_DIR, TRC_ANNUAL_DIR, CLIMATOLOGY_DIR,
    YEAR_START, YEAR_END, BLOCK_SIZE, MAX_WORKERS, NODATA_OUT, TR_FILE_FORMAT,
    PHENO_FILE_FORMAT, TEMPLATE_RASTER, MASK_FILE,
    GPP_DAILY_DIR, GPP_DAILY_FORMAT,  # GPP/NDVI相关配置（已重定向到NDVI）
    OUTPUT_CUMULATIVE_FORMAT, OUTPUT_CLIMATOLOGY_FORMAT,  # 输出文件名格式
    MIDDLE_VAR_NAME  # 用于日志标签
)

# 确保输出目录存在
TRC_ANNUAL_DIR.mkdir(parents=True, exist_ok=True)
CLIMATOLOGY_DIR.mkdir(parents=True, exist_ok=True)

# 向后兼容：保留OUTPUT_DIR别名
OUTPUT_DIR = TRC_ANNUAL_DIR

# 运行模式控制
RUN_MODE = os.getenv("WANG_RUN_MODE", "skip").strip().lower()
OVERWRITE = RUN_MODE == "overwrite"

def should_write(path):
    path = Path(path)
    return OVERWRITE or not path.exists()

def outputs_complete(years):
    trc_done = all((OUTPUT_DIR / f"TRc_{year}.tif").exists() for year in years)
    gppc_done = all((OUTPUT_DIR / OUTPUT_CUMULATIVE_FORMAT.format(year=year)).exists() for year in years)
    clim_done = all(
        (CLIMATOLOGY_DIR / name).exists()
        for name in ("TR_daily_climatology.tif", "POSav.tif", "EOSav.tif",
                     OUTPUT_CLIMATOLOGY_FORMAT['daily'], OUTPUT_CLIMATOLOGY_FORMAT['cumulative_av'])
    )
    return trc_done and gppc_done and clim_done

# ==================== 辅助函数 ====================
def get_TR_file_path(date_obj):
    """
    获取TR日数据文件路径（使用_config.py中的TR_FILE_FORMAT）
    """
    yyyymmdd = date_obj.strftime("%Y%m%d")
    file_path = TR_DAILY_DIR / TR_FILE_FORMAT.format(date=yyyymmdd)
    if file_path.exists():
        return file_path
    return None


def load_template_profile():
    if not TEMPLATE_RASTER.exists():
        raise FileNotFoundError(f"模板栅格不存在: {TEMPLATE_RASTER}")
    with rasterio.open(TEMPLATE_RASTER) as src:
        return src.profile.copy()

def load_mask():
    if not MASK_FILE.exists():
        raise FileNotFoundError(f"掩膜文件不存在: {MASK_FILE}")
    with rasterio.open(MASK_FILE) as src:
        mask_data = src.read(1)
        mask_nodata = src.nodata
    mask = np.isfinite(mask_data)
    if mask_nodata is not None and np.isfinite(mask_nodata):
        mask &= mask_data != mask_nodata
    mask &= mask_data > 0
    return mask


def check_spatial_consistency(ref_profile, file_path, data_profile, var_name="data"):
    if ref_profile["width"] != data_profile["width"] or ref_profile["height"] != data_profile["height"]:
        raise ValueError(
            f"{var_name}尺寸不一致: {file_path}\n"
            f"  Expected: {ref_profile['height']}x{ref_profile['width']}\n"
            f"  Got: {data_profile['height']}x{data_profile['width']}"
        )
    if ref_profile.get("crs") != data_profile.get("crs"):
        raise ValueError(
            f"{var_name} CRS不一致: {file_path}\n"
            f"  Expected: {ref_profile.get('crs')}\n"
            f"  Got: {data_profile.get('crs')}"
        )
    ref_transform = ref_profile.get("transform")
    data_transform = data_profile.get("transform")
    if ref_transform is not None and data_transform is not None:
        transform_match = all(
            abs(ref_transform[i] - data_transform[i]) < 1e-6
            for i in range(6)
        )
        if not transform_match:
            raise ValueError(
                f"{var_name} Transform不一致: {file_path}\n"
                f"  Expected: {ref_transform}\n"
                f"  Got: {data_transform}"
            )

def is_leap_year(year):
    """判断是否为闰年"""
    return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)

# 物候DOY使用无闰日（1-365）
PHENO_MAX_DOY = 365
MIN_VALID_FRAC = 0.60  # 有效天数阈值（POS-EOS窗口内）

def noleap_doy_from_date(date_obj):
    """
    将真实日期映射为无闰日DOY（1-365）。
    闰年跳过2月29日，3月1日及之后的DOY减1。
    返回None表示2月29日（需跳过）。
    """
    doy = date_obj.timetuple().tm_yday
    if is_leap_year(date_obj.year):
        if doy == 60:
            return None
        if doy > 60:
            return doy - 1
    return doy

# ==================== 块处理版本（优化内存+性能） ====================
def calculate_TRc_block_optimized(year, mask, template_profile):
    """
    块处理优化版本：先天后块循环 + NODATA修复 + valid_pheno检查

    优化要点：
    1. 循环顺序改为"先天后块"，每天只打开TR文件一次（性能提升18倍+）
    2. 使用栅格真实nodata值，而非硬编码（修复数据正确性）
    3. 增加valid_pheno检查，避免"无效物候输出为0"（修复数据正确性）

    Parameters:
    -----------
    year : int
        年份
    mask : ndarray
        有效掩膜（布尔数组）
    template_profile : dict
        模板栅格profile（统一网格）

    Returns:
    --------
    TRc : ndarray
        累积蒸腾栅格
    profile : dict
        栅格配置
    """
    print(f"\n=== 计算TRc（块处理-加速版）: {year} ===")

    # 读取物候数据（物候代码输出格式）
    pos_file = PHENO_DIR / PHENO_FILE_FORMAT['POS'].format(year=year)
    eos_file = PHENO_DIR / PHENO_FILE_FORMAT['EOS'].format(year=year)

    if not pos_file.exists() or not eos_file.exists():
        print(f"  ✗ 错误：缺少物候数据")
        print(f"    POS: {pos_file} ({'存在' if pos_file.exists() else '不存在'})")
        print(f"    EOS: {eos_file} ({'存在' if eos_file.exists() else '不存在'})")
        return None, None

    with rasterio.open(pos_file) as src:
        pos = src.read(1)
        height, width = src.height, src.width
        nodata_pos = src.nodata  # ✅ 读取真实nodata
    with rasterio.open(eos_file) as src:
        eos = src.read(1)
        nodata_eos = src.nodata  # ✅ 读取真实nodata

    # 使用统一模板profile
    profile = template_profile.copy()
    test_date = datetime(year, 1, 15)
    test_tr_file = get_TR_file_path(test_date)
    if test_tr_file is None or not test_tr_file.exists():
        print(f"  ✗ 错误：找不到TR数据文件用于检查profile一致性")
        return None, None
    with rasterio.open(test_tr_file) as src:
        check_spatial_consistency(profile, test_tr_file, src.profile, f"TR({year})")

    days_in_year = 366 if is_leap_year(year) else 365
    dates_year = [datetime(year, 1, 1) + timedelta(days=i) for i in range(days_in_year)]

    # ✅ 修复：增加valid_pheno检查，避免"无效物候输出为0"
    def _is_valid_pheno(arr, nodata):
        """检查物候值是否有效"""
        if nodata is None:
            return np.isfinite(arr)
        # nodata 可能是 nan
        if np.isnan(nodata):
            return np.isfinite(arr)
        return (arr != nodata) & np.isfinite(arr)

    valid_pheno = (
        mask
        & _is_valid_pheno(pos, nodata_pos)
        & _is_valid_pheno(eos, nodata_eos)
        & (pos > 0) & (eos > 0)
        & (eos >= pos)
        & (pos <= PHENO_MAX_DOY) & (eos <= PHENO_MAX_DOY)
    )

    # 初始化TRc：只对valid_pheno初始化为0，其他为NODATA
    TRc = np.full((height, width), NODATA_OUT, dtype=np.float32)
    TRc[valid_pheno] = 0.0
    TRc_count = np.zeros((height, width), dtype=np.int16)

    # 生成块窗口
    block_windows = [
        Window(col_off, row_off,
               min(BLOCK_SIZE, width - col_off),
               min(BLOCK_SIZE, height - row_off))
        for row_off in range(0, height, BLOCK_SIZE)
        for col_off in range(0, width, BLOCK_SIZE)
    ]
    print(f"  总块数: {len(block_windows)}")

    missing_files = 0
    nodata_tr = None

    # ✅ 优化：改为"先天后块"循环（性能提升18倍+）
    for date_obj in tqdm(dates_year, desc="逐日累加（先天后块）", leave=False):
        tr_file = get_TR_file_path(date_obj)
        if tr_file is None or not tr_file.exists():
            missing_files += 1
            continue

        doy = noleap_doy_from_date(date_obj)
        if doy is None:
            continue

        try:
            with rasterio.open(tr_file) as src:
                # ✅ 每文件读取nodata，避免跨文件不一致
                nodata_tr = src.nodata

                # 逐块读取并累加
                for win in block_windows:
                    r0, r1 = win.row_off, win.row_off + win.height
                    c0, c1 = win.col_off, win.col_off + win.width

                    pos_b = pos[r0:r1, c0:c1]
                    eos_b = eos[r0:r1, c0:c1]
                    vph_b = valid_pheno[r0:r1, c0:c1]

                    # 如果块内没有有效像元，跳过
                    if not np.any(vph_b):
                        continue

                    # 读取块内TR数据
                    tr_b = src.read(1, window=win)

                    # 累加条件：POS ≤ doy ≤ EOS
                    in_window = vph_b & (pos_b <= doy) & (doy <= eos_b)

                    # ✅ 使用真实nodata判断TR有效性
                    if nodata_tr is None:
                        valid_tr = np.isfinite(tr_b)
                    else:
                        if np.isnan(nodata_tr):
                            valid_tr = np.isfinite(tr_b)
                        else:
                            valid_tr = (tr_b != nodata_tr) & np.isfinite(tr_b)
                    valid_tr &= (tr_b >= 0)

                    acc = in_window & valid_tr
                    if np.any(acc):
                        trc_block = TRc[r0:r1, c0:c1]
                        trc_block[acc] += tr_b[acc]
                        count_block = TRc_count[r0:r1, c0:c1]
                        count_block[acc] += 1

        except Exception as e:
            missing_files += 1
            continue

    print(f"  缺失日文件数: {missing_files}/{days_in_year}")

    # 有效天数比例过滤（窗口内有效天数需 >= MIN_VALID_FRAC * 窗口长度）
    pos_i = np.rint(pos).astype(np.int16)
    eos_i = np.rint(eos).astype(np.int16)
    window_len = (eos_i - pos_i + 1).astype(np.int16)
    min_required = np.ceil(MIN_VALID_FRAC * window_len).astype(np.int16)
    valid_count = valid_pheno & (TRc_count >= min_required)
    TRc[valid_pheno & ~valid_count] = NODATA_OUT

    # 质量检查
    TRc_valid = TRc[valid_count & (TRc != NODATA_OUT)]
    if TRc_valid.size > 0:
        print(f"  TRc范围: {TRc_valid.min():.2f} - {TRc_valid.max():.2f} mm")
        print(f"  TRc平均: {TRc_valid.mean():.2f} mm")
        print(f"  有效像元: {TRc_valid.size}")
    else:
        print("  ⚠ 警告：没有有效TRc值（请检查物候/掩膜/TR数据）")

    return TRc, profile

# ==================== 气候态计算（多年平均）====================
def calculate_daily_TR_climatology(years, mask, template_profile):
    """
    计算多年平均日TR时间序列（日气候态）

    对于每个day-of-year (1-365)，计算所有年份的平均TR值
    用于Wang (2025)原版分解方法的TRpheno计算

    Parameters:
    -----------
    years : list
        年份列表
    mask : ndarray
        有效掩膜（布尔数组）
    template_profile : dict
        模板栅格profile（统一网格）

    Returns:
    --------
    TR_climatology : ndarray, shape (365, H, W)
        365天的多年平均TR（单位：mm/day）
    profile : dict
        栅格配置
    """
    print("\n=== 计算多年平均日TR气候态 ===")
    print(f"  年份范围: {min(years)}-{max(years)} ({len(years)}年)")

    profile = template_profile.copy()
    test_date = datetime(years[0], 1, 15)
    test_file = get_TR_file_path(test_date)
    if test_file is None or not test_file.exists():
        print(f"  ✗ 错误：找不到TR数据文件")
        return None, None
    with rasterio.open(test_file) as src:
        check_spatial_consistency(profile, test_file, src.profile, "TR样本")
        height, width = src.height, src.width
        nodata_tr = src.nodata

    # 初始化累加器：(365天, H, W)
    TR_sum = np.zeros((365, height, width), dtype=np.float64)
    TR_count = np.zeros((365, height, width), dtype=np.int16)

    print(f"  栅格尺寸: {height} × {width}")
    print(f"  数据大小: {365 * height * width * 8 / 1e9:.2f} GB (float64)")

    # 逐年累加
    for year in tqdm(years, desc="逐年累加TR"):
        days_in_year = 366 if is_leap_year(year) else 365
        dates_year = [datetime(year, 1, 1) + timedelta(days=i) for i in range(days_in_year)]

        for date_obj in dates_year:
            doy = date_obj.timetuple().tm_yday

            # 闰年处理：跳过2月29日，3月1日之后DOY-2映射到0-based索引
            if is_leap_year(year):
                if doy == 60:
                    continue  # 跳过2月29日（DOY=60）
                elif doy > 60:
                    doy_idx = doy - 2  # 61->59(3月1日), 62->60, ..., 366->364
                else:
                    doy_idx = doy - 1  # 1->0, ..., 59->58(2月28日)
            else:
                doy_idx = doy - 1  # 1-based to 0-based

            tr_file = get_TR_file_path(date_obj)
            if tr_file is None or not tr_file.exists():
                continue

            try:
                with rasterio.open(tr_file) as src:
                    tr_daily = src.read(1)
                    nodata_tr = src.nodata

                # 判断有效TR值
                if nodata_tr is None:
                    valid_tr = mask & np.isfinite(tr_daily)
                else:
                    if np.isnan(nodata_tr):
                        valid_tr = mask & np.isfinite(tr_daily)
                    else:
                        valid_tr = mask & (tr_daily != nodata_tr) & np.isfinite(tr_daily)
                valid_tr &= (tr_daily >= 0)

                # 累加
                tr_sum_block = TR_sum[doy_idx]
                tr_sum_block[valid_tr] += tr_daily[valid_tr]
                tr_count_block = TR_count[doy_idx]
                tr_count_block[valid_tr] += 1

            except Exception as e:
                continue

    # 计算平均值
    print("  计算平均值...")
    TR_climatology = np.full((365, height, width), NODATA_OUT, dtype=np.float32)

    with np.errstate(divide='ignore', invalid='ignore'):
        for doy_idx in range(365):
            valid = (TR_count[doy_idx] > 0) & mask
            tr_clim_block = TR_climatology[doy_idx]
            tr_clim_block[valid] = (TR_sum[doy_idx][valid] /
                                    TR_count[doy_idx][valid]).astype(np.float32)

    # 统计信息
    print(f"\n  统计信息:")
    for doy_idx in [0, 90, 180, 270, 364]:  # 采样几个日期
        valid_data = TR_climatology[doy_idx][mask & (TR_climatology[doy_idx] != NODATA_OUT)]
        if valid_data.size > 0:
            print(f"    DOY {doy_idx+1:3d}: 平均={valid_data.mean():.3f} mm/day, "
                  f"范围=[{valid_data.min():.3f}, {valid_data.max():.3f}], "
                  f"有效像元={valid_data.size}")

    return TR_climatology, profile

def calculate_phenology_climatology(years, mask, template_profile):
    """
    计算多年平均物候（POSav, EOSav）

    Parameters:
    -----------
    years : list
        年份列表
    mask : ndarray
        有效掩膜（布尔数组），用于保持与TR气候态的空间一致性
    template_profile : dict
        模板栅格profile（统一网格）

    Returns:
    --------
    POSav, EOSav : ndarray, shape (H, W)
        多年平均POS和EOS（仅在mask区域内有效）
    profile : dict
        栅格配置
    """
    print("\n=== 计算多年平均物候（POSav, EOSav）===")
    print(f"  年份范围: {min(years)}-{max(years)} ({len(years)}年)")

    pos_stack = []
    eos_stack = []
    profile = template_profile.copy()
    profile_checked = False

    for year in tqdm(years, desc="读取物候数据"):
        pos_file = PHENO_DIR / PHENO_FILE_FORMAT['POS'].format(year=year)
        eos_file = PHENO_DIR / PHENO_FILE_FORMAT['EOS'].format(year=year)

        if not pos_file.exists() or not eos_file.exists():
            print(f"  ⚠ 跳过 {year}：缺少物候数据")
            continue

        with rasterio.open(pos_file) as src:
            pos = src.read(1)
            if not profile_checked:
                check_spatial_consistency(profile, pos_file, src.profile, f"POS({year})")
            nodata_pos = src.nodata
        with rasterio.open(eos_file) as src:
            eos = src.read(1)
            if not profile_checked:
                check_spatial_consistency(profile, eos_file, src.profile, f"EOS({year})")
            nodata_eos = src.nodata
            profile_checked = True

        # 构建有效物候掩膜（与TRc计算逻辑一致）
        def _is_valid_pheno_value(arr, nodata):
            """检查物候值是否有效"""
            if nodata is None:
                return np.isfinite(arr)
            if np.isnan(nodata):
                return np.isfinite(arr)
            return (arr != nodata) & np.isfinite(arr)

        valid = (
            mask &
            _is_valid_pheno_value(pos, nodata_pos) &
            _is_valid_pheno_value(eos, nodata_eos) &
            (pos > 0) & (pos <= 365) &  # 合法DOY范围
            (eos > 0) & (eos <= 365) &  # 合法DOY范围
            (eos >= pos)  # EOS必须>=POS
        )

        # 转换为NaN（无效值设为NaN，不参与多年平均）
        pos_masked = np.where(valid, pos, np.nan)
        eos_masked = np.where(valid, eos, np.nan)

        pos_stack.append(pos_masked)
        eos_stack.append(eos_masked)

    # 计算多年平均
    pos_stack = np.stack(pos_stack, axis=0)
    eos_stack = np.stack(eos_stack, axis=0)

    with np.errstate(invalid='ignore'):
        POSav = np.nanmean(pos_stack, axis=0)
        EOSav = np.nanmean(eos_stack, axis=0)

        # 应用掩膜：仅保留mask区域的数据，保持与TR气候态一致
        POSav = np.where(mask & np.isfinite(POSav), POSav, NODATA_OUT)
        EOSav = np.where(mask & np.isfinite(EOSav), EOSav, NODATA_OUT)

    # 统计
    pos_valid = POSav[(POSav != NODATA_OUT) & mask]
    eos_valid = EOSav[(EOSav != NODATA_OUT) & mask]

    print(f"  POSav统计: 平均={pos_valid.mean():.1f} DOY, "
          f"范围=[{pos_valid.min():.0f}, {pos_valid.max():.0f}], "
          f"有效像元={pos_valid.size}")
    print(f"  EOSav统计: 平均={eos_valid.mean():.1f} DOY, "
          f"范围=[{eos_valid.min():.0f}, {eos_valid.max():.0f}], "
          f"有效像元={eos_valid.size}")

    return POSav, EOSav, profile

def save_climatology_data():
    """
    计算并保存气候态数据（用于原版分解方法）

    输出文件:
    ---------
    Wang2025_Analysis/Climatology/
        ├── TR_daily_climatology.tif  (365波段，每波段为一个DOY的多年平均TR)
        ├── POSav.tif                 (多年平均POS)
        └── EOSav.tif                 (多年平均EOS)
    """
    print("\n" + "="*70)
    print("计算多年平均气候态数据（用于TRc原版分解）")
    print("="*70)

    # 使用配置文件中的输出目录（已在脚本开头创建）
    climatology_dir = CLIMATOLOGY_DIR

    years = list(range(YEAR_START, YEAR_END + 1))
    template_profile = load_template_profile()

    output_tr = climatology_dir / "TR_daily_climatology.tif"
    output_pos = climatology_dir / "POSav.tif"
    output_eos = climatology_dir / "EOSav.tif"
    if not OVERWRITE and output_tr.exists() and output_pos.exists() and output_eos.exists():
        print("  ⚠ 气候态文件已存在，按skip模式跳过计算")
        return

    # 读取掩膜
    print("\n[1/3] 读取掩膜...")
    mask_file = MASK_FILE

    if not mask_file.exists():
        print(f"  ✗ 错误：掩膜文件不存在: {mask_file}")
        print("  请先运行 Module 00: 00_data_preparation.py")
        return

    mask = load_mask()

    print(f"  ✓ 有效像元数: {np.sum(mask)}")

    # 计算TR日气候态
    print("\n[2/3] 计算TR日气候态（365天）...")
    TR_climatology, profile = calculate_daily_TR_climatology(years, mask, template_profile)

    if TR_climatology is None:
        print("  ✗ 计算失败")
        return

    # 保存TR气候态（365波段）
    output_tr = climatology_dir / "TR_daily_climatology.tif"
    profile_tr = profile.copy()
    profile_tr.update(
        dtype=rasterio.float32,
        count=365,  # 365个波段
        compress='lzw',
        nodata=NODATA_OUT
    )

    print(f"  保存TR气候态: {output_tr.name}")
    if should_write(output_tr):
        with rasterio.open(output_tr, 'w', **profile_tr) as dst:
            for doy_idx in range(365):
                dst.write(TR_climatology[doy_idx].astype(np.float32), doy_idx + 1)
        print(f"  ✓ 已保存: {output_tr}")
    else:
        print(f"  ⚠ 跳过已存在: {output_tr}")

    # 计算物候气候态
    print("\n[3/3] 计算物候气候态（POSav, EOSav）...")
    POSav, EOSav, _ = calculate_phenology_climatology(years, mask, template_profile)

    # 保存POSav/EOSav（使用profile_tr以确保CRS和空间参考与TR_daily_climatology一致）
    profile_pheno = profile_tr.copy()
    profile_pheno.update(
        dtype=rasterio.float32,
        count=1,
        compress='lzw',
        nodata=NODATA_OUT
    )
    output_pos = climatology_dir / "POSav.tif"
    if should_write(output_pos):
        with rasterio.open(output_pos, 'w', **profile_pheno) as dst:
            dst.write(POSav.astype(np.float32), 1)
        print(f"  ✓ 已保存: {output_pos}")
    else:
        print(f"  ⚠ 跳过已存在: {output_pos}")
    # 保存EOSav
    if should_write(output_eos):
        with rasterio.open(output_eos, 'w', **profile_pheno) as dst:
            dst.write(EOSav.astype(np.float32), 1)
        print(f"  ✓ 已保存: {output_eos}")
    else:
        print(f"  ⚠ 跳过已存在: {output_eos}")

    print("\n" + "="*70)
    print("✓ 气候态数据计算完成！")
    print(f"输出目录: {climatology_dir}")
    print(f"  - TR_daily_climatology.tif (365波段)")
    print(f"  - POSav.tif")
    print(f"  - EOSav.tif")
    print("="*70)

# ==================== 主处理流程 ====================
def main(use_block_processing=True):
    """
    主处理流程 - 计算TRc和{MIDDLE_VAR_NAME}c

    Parameters:
    -----------
    use_block_processing : bool
        是否使用块处理（普通版已移除，参数仅保留兼容）
    """
    print("\n" + "="*70)
    print(f"Module 02: TRc和{MIDDLE_VAR_NAME}c计算")
    print("="*70)

    years = list(range(YEAR_START, YEAR_END + 1))

    if RUN_MODE == "skip" and outputs_complete(years):
        print("  ✓ 输出齐全，跳过 Module 02")
        return

    # 读取掩膜
    print("\n[1/3] 读取掩膜...")
    mask_file = MASK_FILE

    if not mask_file.exists():
        print(f"  ✗ 错误：掩膜文件不存在: {mask_file}")
        print("  请先运行 Module 00: 00_data_preparation.py")
        return

    mask = load_mask()

    print(f"  ✓ 有效像元数: {np.sum(mask)}")

    # 加载模板配置（用于统一输出网格）
    template_profile = load_template_profile()

    # 检查TR数据
    print("\n[2/3] 检查TR数据...")
    test_date = datetime(YEAR_START, 1, 15)  # ✅ 使用YEAR_START年份测试
    test_file = get_TR_file_path(test_date)

    if test_file is None or not test_file.exists():
        print(f"  ✗ 错误：找不到TR数据")
        print(f"    测试日期: {test_date.strftime('%Y-%m-%d')}")
        print(f"    预期路径示例: {TR_DAILY_DIR / TR_FILE_FORMAT.format(date=f'{YEAR_START}0115')}")
        print(f"\n请检查以下路径是否正确:")
        print(f"  TR_DAILY_DIR = {TR_DAILY_DIR}")
        print(f"\n或修改 get_TR_file_path() 函数以匹配您的文件命名格式")
        return
    else:
        print(f"  ✓ 找到TR数据: {test_file.name}")

    # 逐年计算TRc和{MIDDLE_VAR_NAME}c
    print(f"\n[3/3] 逐年计算TRc和{MIDDLE_VAR_NAME}c ({YEAR_START}-{YEAR_END})...")

    # 统计
    processed_trc = 0
    processed_gppc = 0
    skipped_trc = 0
    skipped_gppc = 0
    failed = 0

    for year in years:
        trc_file = OUTPUT_DIR / f"TRc_{year}.tif"
        gppc_file = OUTPUT_DIR / OUTPUT_CUMULATIVE_FORMAT.format(year=year)

        year_failed = False

        # ========== 计算TRc ==========
        if trc_file.exists() and not OVERWRITE:
            print(f"  ⚠ TRc已存在: {year} ({trc_file.name})")
            skipped_trc += 1
        else:
            # 统一使用块处理（普通版已移除）
            if not use_block_processing:
                print("  ⚠ 已移除普通版，强制使用块处理加速版")
            TRc, profile = calculate_TRc_block_optimized(year, mask, template_profile)

            if TRc is None:
                print(f"  ✗ 跳过年份 {year}（缺少TRc数据）")
                year_failed = True
            else:
                # 保存TRc结果
                profile.update(dtype=rasterio.float32, count=1, compress='lzw', nodata=NODATA_OUT)
                with rasterio.open(trc_file, 'w', **profile) as dst:
                    dst.write(TRc.astype(np.float32), 1)
                print(f"  ✓ 已保存TRc: {trc_file.name}")
                processed_trc += 1

        # ========== 计算{MIDDLE_VAR_NAME}c ==========
        if not year_failed:
            if gppc_file.exists() and not OVERWRITE:
                print(f"  ⚠ {MIDDLE_VAR_NAME}c已存在: {year} ({gppc_file.name})")
                skipped_gppc += 1
            else:
                try:
                    GPPc, gppc_profile = calculate_GPPc_block_optimized(year, mask, template_profile)

                    if GPPc is None:
                        print(f"  ✗ 跳过{MIDDLE_VAR_NAME}c计算 {year}（缺少{MIDDLE_VAR_NAME}数据）")
                    else:
                        # 保存{MIDDLE_VAR_NAME}c结果
                        gppc_profile.update(dtype=rasterio.float32, count=1, compress='lzw', nodata=NODATA_OUT)
                        with rasterio.open(gppc_file, 'w', **gppc_profile) as dst:
                            dst.write(GPPc.astype(np.float32), 1)
                        print(f"  ✓ 已保存{MIDDLE_VAR_NAME}c: {gppc_file.name}")
                        processed_gppc += 1
                except Exception as e:
                    print(f"  ✗ {MIDDLE_VAR_NAME}c计算失败 {year}: {str(e)}")
        else:
            failed += 1

    print("\n" + "="*70)
    print(f"✓ TRc和{MIDDLE_VAR_NAME}c计算完成！")
    print(f"输出目录: {OUTPUT_DIR}")
    print(f"\n统计:")
    print(f"  TRc - 新处理: {processed_trc}, 已跳过: {skipped_trc}")
    print(f"  {MIDDLE_VAR_NAME}c - 新处理: {processed_gppc}, 已跳过: {skipped_gppc}")
    if failed > 0:
        print(f"  - 失败: {failed} 个年份（数据缺失）")
    print(f"  - 总计: {len(years)} 个年份")
    print("="*70)

# ==================== 并行处理版本 ====================
def process_single_year(args):
    """
    处理单个年份（用于并行）- 同时计算TRc和{MIDDLE_VAR_NAME}c

    Parameters:
    -----------
    args : tuple
        (year, mask, use_block_processing)

    Returns:
    --------
    result : dict
        处理结果字典
    """
    year, mask, use_block_processing = args
    trc_file = OUTPUT_DIR / f"TRc_{year}.tif"
    gppc_file = OUTPUT_DIR / OUTPUT_CUMULATIVE_FORMAT.format(year=year)

    messages = []
    trc_status = 'skipped'
    gppc_status = 'skipped'

    try:
        template_profile = load_template_profile()

        # ========== 处理TRc ==========
        if trc_file.exists() and not OVERWRITE:
            messages.append(f"TRc已存在: {trc_file.name}")
            trc_status = 'skipped'
        else:
            TRc, profile = calculate_TRc_block_optimized(year, mask, template_profile)

            if TRc is None:
                messages.append('TRc: 缺少物候数据')
                trc_status = 'failed'
            else:
                # 保存TRc结果
                profile.update(dtype=rasterio.float32, count=1, compress='lzw', nodata=NODATA_OUT)
                with rasterio.open(trc_file, 'w', **profile) as dst:
                    dst.write(TRc.astype(np.float32), 1)
                messages.append(f"TRc已保存: {trc_file.name}")
                trc_status = 'success'

        # ========== 处理GPPc ==========
        if gppc_file.exists() and not OVERWRITE:
            messages.append(f"{MIDDLE_VAR_NAME}c已存在: {gppc_file.name}")
            gppc_status = 'skipped'
        else:
            try:
                GPPc, gppc_profile = calculate_GPPc_block_optimized(year, mask, template_profile)

                if GPPc is None:
                    messages.append(f'{MIDDLE_VAR_NAME}c: 缺少{MIDDLE_VAR_NAME}数据')
                    gppc_status = 'failed'
                else:
                    # 保存GPPc结果
                    gppc_profile.update(dtype=rasterio.float32, count=1, compress='lzw', nodata=NODATA_OUT)
                    with rasterio.open(gppc_file, 'w', **gppc_profile) as dst:
                        dst.write(GPPc.astype(np.float32), 1)
                    messages.append(f"{MIDDLE_VAR_NAME}c已保存: {gppc_file.name}")
                    gppc_status = 'success'
            except Exception as e:
                messages.append(f"{MIDDLE_VAR_NAME}c出错: {str(e)}")
                gppc_status = 'error'

        # 综合状态
        if trc_status == 'success' or gppc_status == 'success':
            overall_status = 'success'
        elif trc_status == 'failed' or gppc_status == 'failed':
            overall_status = 'failed'
        else:
            overall_status = 'skipped'

        return {
            'year': year,
            'status': overall_status,
            'message': '; '.join(messages)
        }

    except Exception as e:
        return {
            'year': year,
            'status': 'error',
            'message': f"处理出错: {str(e)}"
        }

def main_parallel(use_block_processing=True, max_workers=None):
    """
    主处理流程（并行版本）- 同时计算TRc和{MIDDLE_VAR_NAME}c

    Parameters:
    -----------
    use_block_processing : bool
        是否使用块处理（普通版已移除，参数仅保留兼容）
    max_workers : int
        并行进程数（默认使用MAX_WORKERS配置）
    """
    if max_workers is None:
        max_workers = MAX_WORKERS

    print("\n" + "="*70)
    print(f"Module 02: TRc和{MIDDLE_VAR_NAME}c计算（并行版本，{max_workers}个进程）")
    print("="*70)
    if not use_block_processing:
        print("⚠ 已移除普通版，强制使用块处理加速版")

    years = list(range(YEAR_START, YEAR_END + 1))

    if RUN_MODE == "skip" and outputs_complete(years):
        print("  ✓ 输出齐全，跳过 Module 02（并行）")
        return

    # 读取掩膜
    print("\n[1/3] 读取掩膜...")
    mask_file = MASK_FILE

    if not mask_file.exists():
        print(f"  ✗ 错误：掩膜文件不存在: {mask_file}")
        print("  请先运行 Module 00: 00_data_preparation.py")
        return

    mask = load_mask()

    print(f"  ✓ 有效像元数: {np.sum(mask)}")

    # 检查TR数据
    print("\n[2/3] 检查TR数据...")
    test_date = datetime(YEAR_START, 1, 15)
    test_file = get_TR_file_path(test_date)

    if test_file is None or not test_file.exists():
        print(f"  ✗ 错误：找不到TR数据")
        print(f"    测试日期: {test_date.strftime('%Y-%m-%d')}")
        print(f"    预期路径示例: {TR_DAILY_DIR / TR_FILE_FORMAT.format(date=f'{YEAR_START}0115')}")
        print(f"\n请检查以下路径是否正确:")
        print(f"  TR_DAILY_DIR = {TR_DAILY_DIR}")
        return
    else:
        print(f"  ✓ 找到TR数据: {test_file.name}")

    # 并行计算TRc和{MIDDLE_VAR_NAME}c
    print(f"\n[3/3] 并行计算TRc和{MIDDLE_VAR_NAME}c ({YEAR_START}-{YEAR_END}, {max_workers}进程)...")

    # 准备参数
    args_list = [(year, mask, use_block_processing) for year in years]

    # 统计
    processed = 0
    skipped = 0
    failed = 0

    # 并行处理
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_single_year, args): args[0]
                   for args in args_list}

        for future in tqdm(as_completed(futures), total=len(futures), desc="处理进度"):
            result = future.result()

            if result['status'] == 'success':
                print(f"  ✓ {result['year']}: {result['message']}")
                processed += 1
            elif result['status'] == 'skipped':
                print(f"  ⚠ {result['year']}: {result['message']}")
                skipped += 1
            elif result['status'] == 'failed':
                print(f"  ✗ {result['year']}: {result['message']}")
                failed += 1
            else:  # error
                print(f"  ✗ {result['year']}: {result['message']}")
                failed += 1

    print("\n" + "="*70)
    print(f"✓ TRc和{MIDDLE_VAR_NAME}c计算完成！")
    print(f"输出目录: {OUTPUT_DIR}")
    print(f"\n统计:")
    print(f"  - 新处理: {processed} 个年份")
    print(f"  - 已跳过: {skipped} 个年份（文件已存在）")
    if failed > 0:
        print(f"  - 失败: {failed} 个年份")
    print(f"  - 总计: {len(years)} 个年份")
    print("="*70)


# ==================== {MIDDLE_VAR_NAME}相关函数（新增） ====================

def get_GPP_file_path(date_obj):
    """
    获取{MIDDLE_VAR_NAME}日数据文件路径（使用_config.py中的GPP_DAILY_FORMAT）
    """
    yyyymmdd = date_obj.strftime("%Y%m%d")
    file_path = GPP_DAILY_DIR / GPP_DAILY_FORMAT.format(date=yyyymmdd)
    if file_path.exists():
        return file_path
    return None


def calculate_GPPc_block_optimized(year, mask, template_profile):
    """
    块处理优化版本：计算GPPc/NDVIc（POS-EOS累积）

    完全对齐calculate_TRc_block_optimized的逻辑：
    - 使用无闰日DOY映射（noleap_doy_from_date）
    - 使用块处理（先天后块）
    - 有效天数阈值过滤（MIN_VALID_FRAC）

    Parameters:
    -----------
    year : int
        年份
    mask : ndarray
        研究区掩膜
    template_profile : dict
        rasterio配置模板

    Returns:
    --------
    GPPc : ndarray
        累积GPP/NDVI
    output_profile : dict
        输出配置
    """
    print(f"\n=== 计算{MIDDLE_VAR_NAME}c（块处理-加速版）: {year} ===")

    # 读取物候数据（EOS版本：POS和EOS）
    pos_file = PHENO_DIR / PHENO_FILE_FORMAT['POS'].format(year=year)
    eos_file = PHENO_DIR / PHENO_FILE_FORMAT['EOS'].format(year=year)

    if not pos_file.exists() or not eos_file.exists():
        print(f"  ✗ 错误：缺少物候数据")
        print(f"    POS: {pos_file} ({'存在' if pos_file.exists() else '不存在'})")
        print(f"    EOS: {eos_file} ({'存在' if eos_file.exists() else '不存在'})")
        return None, None

    with rasterio.open(pos_file) as src:
        pos = src.read(1)
        height, width = src.height, src.width
        nodata_pos = src.nodata

    with rasterio.open(eos_file) as src:
        eos = src.read(1)
        nodata_eos = src.nodata

    # 使用统一模板profile
    profile = template_profile.copy()
    test_date = datetime(year, 1, 15)
    test_gpp_file = get_GPP_file_path(test_date)
    if test_gpp_file is None or not test_gpp_file.exists():
        print(f"  ✗ 错误：找不到{MIDDLE_VAR_NAME}数据文件用于检查profile一致性")
        return None, None
    with rasterio.open(test_gpp_file) as src:
        check_spatial_consistency(profile, test_gpp_file, src.profile, f"{MIDDLE_VAR_NAME}({year})")

    days_in_year = 366 if is_leap_year(year) else 365
    dates_year = [datetime(year, 1, 1) + timedelta(days=i) for i in range(days_in_year)]

    # 有效物候检查（与TRc一致）
    def _is_valid_pheno(arr, nodata):
        """检查物候值是否有效"""
        if nodata is None:
            return np.isfinite(arr)
        if np.isnan(nodata):
            return np.isfinite(arr)
        return (arr != nodata) & np.isfinite(arr)

    valid_pheno = (
        mask
        & _is_valid_pheno(pos, nodata_pos)
        & _is_valid_pheno(eos, nodata_eos)
        & (pos > 0) & (eos > 0)
        & (eos >= pos)
        & (pos <= PHENO_MAX_DOY) & (eos <= PHENO_MAX_DOY)
    )

    # 初始化GPPc：只对valid_pheno初始化为0，其他为NODATA
    GPPc = np.full((height, width), NODATA_OUT, dtype=np.float32)
    GPPc[valid_pheno] = 0.0
    GPPc_count = np.zeros((height, width), dtype=np.int16)

    # 生成块窗口
    block_windows = [
        Window(col_off, row_off,
               min(BLOCK_SIZE, width - col_off),
               min(BLOCK_SIZE, height - row_off))
        for row_off in range(0, height, BLOCK_SIZE)
        for col_off in range(0, width, BLOCK_SIZE)
    ]
    print(f"  总块数: {len(block_windows)}")

    missing_files = 0
    nodata_gpp = None

    # 优化：改为"先天后块"循环（与TRc一致）
    for date_obj in tqdm(dates_year, desc=f"逐日累加{MIDDLE_VAR_NAME}（先天后块）", leave=False):
        gpp_file = get_GPP_file_path(date_obj)
        if gpp_file is None or not gpp_file.exists():
            missing_files += 1
            continue

        # 使用无闰日DOY映射（与TRc一致）
        doy = noleap_doy_from_date(date_obj)
        if doy is None:
            continue  # 跳过2月29日

        try:
            with rasterio.open(gpp_file) as src:
                nodata_gpp = src.nodata

                # 逐块读取并累加
                for win in block_windows:
                    r0, r1 = win.row_off, win.row_off + win.height
                    c0, c1 = win.col_off, win.col_off + win.width

                    pos_b = pos[r0:r1, c0:c1]
                    eos_b = eos[r0:r1, c0:c1]
                    vph_b = valid_pheno[r0:r1, c0:c1]

                    # 如果块内没有有效像元，跳过
                    if not np.any(vph_b):
                        continue

                    # 读取块内GPP数据
                    gpp_b = src.read(1, window=win)

                    # 累加条件：POS ≤ doy ≤ EOS
                    in_window = vph_b & (pos_b <= doy) & (doy <= eos_b)

                    # 使用真实nodata判断GPP有效性
                    if nodata_gpp is None:
                        valid_gpp = np.isfinite(gpp_b)
                    else:
                        if np.isnan(nodata_gpp):
                            valid_gpp = np.isfinite(gpp_b)
                        else:
                            valid_gpp = (gpp_b != nodata_gpp) & np.isfinite(gpp_b)
                    valid_gpp &= (gpp_b >= 0)  # NDVI负值过滤，与其他模块保持一致

                    acc = in_window & valid_gpp
                    if np.any(acc):
                        gppc_block = GPPc[r0:r1, c0:c1]
                        gppc_block[acc] += gpp_b[acc]
                        count_block = GPPc_count[r0:r1, c0:c1]
                        count_block[acc] += 1

        except Exception as e:
            missing_files += 1
            continue

    print(f"  缺失日文件数: {missing_files}/{days_in_year}")

    # 有效天数比例过滤（与TRc一致：窗口内有效天数需 >= MIN_VALID_FRAC * 窗口长度）
    pos_i = np.rint(pos).astype(np.int16)
    eos_i = np.rint(eos).astype(np.int16)
    window_len = (eos_i - pos_i + 1).astype(np.int16)
    min_required = np.ceil(MIN_VALID_FRAC * window_len).astype(np.int16)
    valid_count = valid_pheno & (GPPc_count >= min_required)
    GPPc[valid_pheno & ~valid_count] = NODATA_OUT

    # 质量检查
    GPPc_valid = GPPc[valid_count & (GPPc != NODATA_OUT)]
    if GPPc_valid.size > 0:
        print(f"  {MIDDLE_VAR_NAME}c范围: {GPPc_valid.min():.4f} - {GPPc_valid.max():.4f}")
        print(f"  {MIDDLE_VAR_NAME}c平均: {GPPc_valid.mean():.4f}")
        print(f"  有效像元: {GPPc_valid.size}")
    else:
        print(f"  ⚠ 警告：没有有效{MIDDLE_VAR_NAME}c值（请检查物候/掩膜/{MIDDLE_VAR_NAME}数据）")

    # 更新profile
    output_profile = profile.copy()
    output_profile.update({
        'dtype': 'float32',
        'nodata': NODATA_OUT
    })

    return GPPc, output_profile


def save_climatology_data_GPP():
    """
    计算并保存GPP/NDVI气候态数据（EOS版本）

    对齐calculate_daily_TR_climatology的逻辑：
    - 使用无闰日DOY映射（闰年跳过2月29日，3月1日后DOY-1）

    输出：
    -------
    - {MIDDLE_VAR_NAME}_daily_climatology.tif : 365波段，每个DOY的多年平均值
    - {MIDDLE_VAR_NAME}c_av.tif : 基于[POSav, EOSav]的多年平均累积值
    """
    print("\n" + "="*70)
    print(f"计算{MIDDLE_VAR_NAME}气候态数据（EOS版本）")
    print("="*70)

    years = list(range(YEAR_START, YEAR_END + 1))

    # 读取模板和掩膜
    template_profile = load_template_profile()
    mask = load_mask()
    height, width = template_profile['height'], template_profile['width']

    # 1. 计算365天的气候态（与TR气候态逻辑对齐）
    print(f"\n[1/3] 计算{MIDDLE_VAR_NAME}日气候态...")

    # 初始化累加器：(365天, H, W)
    GPP_sum = np.zeros((365, height, width), dtype=np.float64)
    GPP_count = np.zeros((365, height, width), dtype=np.int16)

    print(f"  栅格尺寸: {height} × {width}")
    print(f"  数据大小: {365 * height * width * 8 / 1e9:.2f} GB (float64)")

    # 逐年累加（与TR气候态逻辑对齐）
    for year in tqdm(years, desc=f"逐年累加{MIDDLE_VAR_NAME}"):
        days_in_year = 366 if is_leap_year(year) else 365
        dates_year = [datetime(year, 1, 1) + timedelta(days=i) for i in range(days_in_year)]

        for date_obj in dates_year:
            doy = date_obj.timetuple().tm_yday

            # 闰年处理：跳过2月29日，3月1日之后DOY-2映射到0-based索引（与TR一致）
            if is_leap_year(year):
                if doy == 60:
                    continue  # 跳过2月29日（DOY=60）
                elif doy > 60:
                    doy_idx = doy - 2  # 61->59(3月1日), 62->60, ..., 366->364
                else:
                    doy_idx = doy - 1  # 1->0, ..., 59->58(2月28日)
            else:
                doy_idx = doy - 1  # 1-based to 0-based

            gpp_file = get_GPP_file_path(date_obj)
            if gpp_file is None or not gpp_file.exists():
                continue

            try:
                with rasterio.open(gpp_file) as src:
                    gpp_daily = src.read(1)
                    nodata_gpp = src.nodata

                # 判断有效GPP值
                if nodata_gpp is None:
                    valid_gpp = mask & np.isfinite(gpp_daily)
                else:
                    if np.isnan(nodata_gpp):
                        valid_gpp = mask & np.isfinite(gpp_daily)
                    else:
                        valid_gpp = mask & (gpp_daily != nodata_gpp) & np.isfinite(gpp_daily)
                valid_gpp &= (gpp_daily >= 0)  # NDVI负值过滤，与其他模块保持一致

                # 累加
                gpp_sum_block = GPP_sum[doy_idx]
                gpp_sum_block[valid_gpp] += gpp_daily[valid_gpp]
                gpp_count_block = GPP_count[doy_idx]
                gpp_count_block[valid_gpp] += 1

            except Exception as e:
                continue

    # 计算平均值
    print("  计算平均值...")
    gpp_daily_av = np.full((365, height, width), NODATA_OUT, dtype=np.float32)

    with np.errstate(divide='ignore', invalid='ignore'):
        for doy_idx in range(365):
            valid = (GPP_count[doy_idx] > 0) & mask
            gpp_clim_block = gpp_daily_av[doy_idx]
            gpp_clim_block[valid] = (GPP_sum[doy_idx][valid] /
                                     GPP_count[doy_idx][valid]).astype(np.float32)

    # 统计信息（与TR气候态格式一致）
    print(f"\n  统计信息:")
    for doy_idx in [0, 90, 180, 270, 364]:  # 采样几个日期
        valid_data = gpp_daily_av[doy_idx][mask & (gpp_daily_av[doy_idx] != NODATA_OUT)]
        if valid_data.size > 0:
            print(f"    DOY {doy_idx+1:3d}: 平均={valid_data.mean():.4f}, "
                  f"范围=[{valid_data.min():.4f}, {valid_data.max():.4f}], "
                  f"有效像元={valid_data.size}")

    print(f"  ✓ {MIDDLE_VAR_NAME}日气候态计算完成")
    print(f"    有效DOY: {np.sum(np.any(gpp_daily_av != NODATA_OUT, axis=(1,2)))} / 365")

    # 保存日气候态
    gpp_clim_file = CLIMATOLOGY_DIR / OUTPUT_CLIMATOLOGY_FORMAT['daily']
    output_profile = template_profile.copy()
    output_profile.update({
        'count': 365,
        'dtype': 'float32',
        'nodata': NODATA_OUT
    })

    with rasterio.open(gpp_clim_file, 'w', **output_profile) as dst:
        for doy in range(365):
            band_data = gpp_daily_av[doy].copy()
            band_data[~np.isfinite(band_data)] = NODATA_OUT
            dst.write(band_data, doy + 1)

    print(f"  ✓ 已保存: {gpp_clim_file}")

    # 2. 读取物候气候态
    print("\n[2/3] 读取物候气候态...")
    pos_av_file = CLIMATOLOGY_DIR / "POSav.tif"
    eos_av_file = CLIMATOLOGY_DIR / "EOSav.tif"

    if not pos_av_file.exists() or not eos_av_file.exists():
        print("  警告：物候气候态文件不存在，需要先运行save_climatology_data()")
        return

    with rasterio.open(pos_av_file) as src:
        pos_av = src.read(1)
        pos_nodata = src.nodata

    with rasterio.open(eos_av_file) as src:
        eos_av = src.read(1)
        eos_nodata = src.nodata

    # 3. 计算累积气候态（基于POSav和EOSav）
    print(f"\n[3/3] 计算{MIDDLE_VAR_NAME}c_av...")
    GPPc_av = np.full((height, width), NODATA_OUT, dtype=np.float32)

    valid_pheno = np.isfinite(pos_av) & np.isfinite(eos_av)
    if pos_nodata is not None:
        valid_pheno &= (pos_av != pos_nodata)
    if eos_nodata is not None:
        valid_pheno &= (eos_av != eos_nodata)
    valid_pheno &= (pos_av > 0) & (pos_av <= 365)
    valid_pheno &= (eos_av > 0) & (eos_av <= 365)
    valid_pheno &= (pos_av < eos_av)
    valid_pheno &= mask

    # 逐像元累加（添加MIN_VALID_FRAC过滤，与TRc一致）
    for i in tqdm(range(height), desc=f"逐像元累加{MIDDLE_VAR_NAME}"):
        for j in range(width):
            if not valid_pheno[i, j]:
                continue

            pos = int(np.round(pos_av[i, j]))
            eos = int(np.round(eos_av[i, j]))

            if pos < 1 or pos > 365 or eos < 1 or eos > 365 or pos >= eos:
                continue

            # 累加[POS, EOS]范围的GPP
            gpp_sum = 0.0
            valid_days = 0
            window_len = eos - pos + 1
            for doy in range(pos, eos + 1):
                gpp_val = gpp_daily_av[doy - 1, i, j]
                if gpp_val != NODATA_OUT and np.isfinite(gpp_val):
                    gpp_sum += gpp_val
                    valid_days += 1

            # MIN_VALID_FRAC过滤：有效天数需达到窗口长度的MIN_VALID_FRAC比例
            min_required = int(np.ceil(MIN_VALID_FRAC * window_len))
            if valid_days >= min_required:
                GPPc_av[i, j] = gpp_sum

    valid_gppav = (GPPc_av != NODATA_OUT) & np.isfinite(GPPc_av)
    print(f"  ✓ {MIDDLE_VAR_NAME}c_av计算完成")
    print(f"    有效像元: {np.sum(valid_gppav)}")
    if np.any(valid_gppav):
        print(f"    {MIDDLE_VAR_NAME}c_av范围: {np.min(GPPc_av[valid_gppav]):.2f} - {np.max(GPPc_av[valid_gppav]):.2f}")

    # 保存累积气候态
    gppav_file = CLIMATOLOGY_DIR / OUTPUT_CLIMATOLOGY_FORMAT['cumulative_av']
    output_profile_single = template_profile.copy()
    output_profile_single.update({
        'count': 1,
        'dtype': 'float32',
        'nodata': NODATA_OUT
    })

    with rasterio.open(gppav_file, 'w', **output_profile_single) as dst:
        dst.write(GPPc_av, 1)

    print(f"  ✓ 已保存: {gppav_file}")
    print("\n" + "="*70)
    print(f"✓ {MIDDLE_VAR_NAME}气候态数据计算完成！")
    print(f"输出目录: {CLIMATOLOGY_DIR}")
    print("="*70)


# ==================== 主程序 ====================
if __name__ == "__main__":
    years = list(range(YEAR_START, YEAR_END + 1))
    if RUN_MODE == "skip" and outputs_complete(years):
        print("✓ Module 02 输出齐全，已跳过全部计算")
    else:
        # 步骤1：计算年度TRc（POS-EOS累积蒸腾）
        # 方式1：串行处理（推荐，稳定可靠）
        main(use_block_processing=True)

        # 方式2：并行处理（可选，SSD环境下可能更快）
        # ⚠️ 注意：
        # - 并行适用于SSD，HDD可能适得其反
        # - MAX_WORKERS=2-4 比较稳妥（I/O密集任务）
        # - 自动支持断点续算
        # main_parallel(use_block_processing=True, max_workers=2)

        # ======================================================================
        # 步骤2：计算气候态数据（用于03_decomposition_original.py原版分解）
        # ⚠️ 重要：仅在完成步骤1后，且准备进行分解时再运行此步骤
        # ======================================================================
        save_climatology_data()
        #
        # 此步骤将计算并保存：
        #   - TR_daily_climatology.tif (365波段，每个DOY的多年平均TR)
        #   - POSav.tif (多年平均POS)
        #   - EOSav.tif (多年平均EOS)
        #
        # 输出目录：Wang2025_Analysis/Climatology/
        #
        # 注意：此步骤计算密集，预计需要30-60分钟（取决于数据规模）

        # ======================================================================
        # 步骤3：计算GPP气候态数据（新增）
        # ======================================================================
        save_climatology_data_GPP()
        #
        # 此步骤将计算并保存：
        #   - NDVI_daily_climatology.tif (365波段，每个DOY的多年平均NDVI)
        #   - NDVIc_av.tif (基于POSav和EOSav的多年平均累积NDVI)
        #
        # 输出目录：Wang2025_Analysis_EOS_GPP_NDVI_GLEAM/Climatology/
        #
        # 注意：此步骤需要在save_climatology_data()之后运行（依赖POSav和EOSav）

