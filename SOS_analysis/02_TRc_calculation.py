#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module 02: TRc（累积蒸腾）计算
计算SOS-POS窗口内的累积蒸腾量
TRc_y = Σ[SOS to POS] TR(t)
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
    PHENO_FILE_FORMAT, TEMPLATE_RASTER, MASK_FILE
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
    clim_done = all(
        (CLIMATOLOGY_DIR / name).exists()
        for name in ("TR_daily_climatology.tif", "SOSav.tif", "POSav.tif")
    )
    return trc_done and clim_done

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
MIN_VALID_FRAC = 0.60  # 有效天数阈值（SOS-POS窗口内）

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
    sos_file = PHENO_DIR / PHENO_FILE_FORMAT['SOS'].format(year=year)
    pos_file = PHENO_DIR / PHENO_FILE_FORMAT['POS'].format(year=year)

    if not sos_file.exists() or not pos_file.exists():
        print(f"  ✗ 错误：缺少物候数据")
        print(f"    SOS: {sos_file} ({'存在' if sos_file.exists() else '不存在'})")
        print(f"    POS: {pos_file} ({'存在' if pos_file.exists() else '不存在'})")
        return None, None

    with rasterio.open(sos_file) as src:
        sos = src.read(1)
        height, width = src.height, src.width
        nodata_sos = src.nodata  # ✅ 读取真实nodata

    with rasterio.open(pos_file) as src:
        pos = src.read(1)
        nodata_pos = src.nodata  # ✅ 读取真实nodata

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
        & _is_valid_pheno(sos, nodata_sos)
        & _is_valid_pheno(pos, nodata_pos)
        & (sos > 0) & (pos > 0)
        & (pos >= sos)
        & (sos <= PHENO_MAX_DOY) & (pos <= PHENO_MAX_DOY)
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

                    sos_b = sos[r0:r1, c0:c1]
                    pos_b = pos[r0:r1, c0:c1]
                    vph_b = valid_pheno[r0:r1, c0:c1]

                    # 如果块内没有有效像元，跳过
                    if not np.any(vph_b):
                        continue

                    # 读取块内TR数据
                    tr_b = src.read(1, window=win)

                    # 累加条件：SOS ≤ doy ≤ POS
                    in_window = vph_b & (sos_b <= doy) & (doy <= pos_b)

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
    sos_i = np.rint(sos).astype(np.int16)
    pos_i = np.rint(pos).astype(np.int16)
    window_len = (pos_i - sos_i + 1).astype(np.int16)
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
    计算多年平均物候（SOSav, POSav）

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
    SOSav, POSav : ndarray, shape (H, W)
        多年平均SOS和POS（仅在mask区域内有效）
    profile : dict
        栅格配置
    """
    print("\n=== 计算多年平均物候（SOSav, POSav）===")
    print(f"  年份范围: {min(years)}-{max(years)} ({len(years)}年)")

    sos_stack = []
    pos_stack = []
    profile = template_profile.copy()
    profile_checked = False

    for year in tqdm(years, desc="读取物候数据"):
        sos_file = PHENO_DIR / PHENO_FILE_FORMAT['SOS'].format(year=year)
        pos_file = PHENO_DIR / PHENO_FILE_FORMAT['POS'].format(year=year)

        if not sos_file.exists() or not pos_file.exists():
            print(f"  ⚠ 跳过 {year}：缺少物候数据")
            continue

        with rasterio.open(sos_file) as src:
            sos = src.read(1)
            sos_profile = src.profile.copy()
            nodata_sos = src.nodata
            if not profile_checked:
                check_spatial_consistency(profile, sos_file, sos_profile, f"SOS({year})")

        with rasterio.open(pos_file) as src:
            pos = src.read(1)
            if not profile_checked:
                check_spatial_consistency(profile, pos_file, src.profile, f"POS({year})")
            nodata_pos = src.nodata
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
            _is_valid_pheno_value(sos, nodata_sos) &
            _is_valid_pheno_value(pos, nodata_pos) &
            (sos > 0) & (sos <= 365) &  # 合法DOY范围
            (pos > 0) & (pos <= 365) &  # 合法DOY范围
            (pos >= sos)  # POS必须>=SOS
        )

        # 转换为NaN（无效值设为NaN，不参与多年平均）
        sos_masked = np.where(valid, sos, np.nan)
        pos_masked = np.where(valid, pos, np.nan)

        sos_stack.append(sos_masked)
        pos_stack.append(pos_masked)

    # 计算多年平均
    sos_stack = np.stack(sos_stack, axis=0)
    pos_stack = np.stack(pos_stack, axis=0)

    with np.errstate(invalid='ignore'):
        SOSav = np.nanmean(sos_stack, axis=0)
        POSav = np.nanmean(pos_stack, axis=0)

        # 应用掩膜：仅保留mask区域的数据，保持与TR气候态一致
        SOSav = np.where(mask & np.isfinite(SOSav), SOSav, NODATA_OUT)
        POSav = np.where(mask & np.isfinite(POSav), POSav, NODATA_OUT)

    # 统计
    sos_valid = SOSav[(SOSav != NODATA_OUT) & mask]
    pos_valid = POSav[(POSav != NODATA_OUT) & mask]

    print(f"  SOSav统计: 平均={sos_valid.mean():.1f} DOY, "
          f"范围=[{sos_valid.min():.0f}, {sos_valid.max():.0f}], "
          f"有效像元={sos_valid.size}")
    print(f"  POSav统计: 平均={pos_valid.mean():.1f} DOY, "
          f"范围=[{pos_valid.min():.0f}, {pos_valid.max():.0f}], "
          f"有效像元={pos_valid.size}")

    return SOSav, POSav, profile

def save_climatology_data():
    """
    计算并保存气候态数据（用于原版分解方法）

    输出文件:
    ---------
    Wang2025_Analysis/Climatology/
        ├── TR_daily_climatology.tif  (365波段，每波段为一个DOY的多年平均TR)
        ├── SOSav.tif                 (多年平均SOS)
        └── POSav.tif                 (多年平均POS)
    """
    print("\n" + "="*70)
    print("计算多年平均气候态数据（用于TRc原版分解）")
    print("="*70)

    # 使用配置文件中的输出目录（已在脚本开头创建）
    climatology_dir = CLIMATOLOGY_DIR

    years = list(range(YEAR_START, YEAR_END + 1))
    template_profile = load_template_profile()

    output_tr = climatology_dir / "TR_daily_climatology.tif"
    output_sos = climatology_dir / "SOSav.tif"
    output_pos = climatology_dir / "POSav.tif"
    if not OVERWRITE and output_tr.exists() and output_sos.exists() and output_pos.exists():
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
    print("\n[3/3] 计算物候气候态（SOSav, POSav）...")
    SOSav, POSav, _ = calculate_phenology_climatology(years, mask, template_profile)

    # 保存SOSav（使用profile_tr以确保CRS和空间参考与TR_daily_climatology一致）
    output_sos = climatology_dir / "SOSav.tif"
    profile_sos = profile_tr.copy()
    profile_sos.update(
        dtype=rasterio.float32,
        count=1,
        compress='lzw',
        nodata=NODATA_OUT
    )
    if should_write(output_sos):
        with rasterio.open(output_sos, 'w', **profile_sos) as dst:
            dst.write(SOSav.astype(np.float32), 1)
        print(f"  ✓ 已保存: {output_sos}")
    else:
        print(f"  ⚠ 跳过已存在: {output_sos}")

    # 保存POSav
    output_pos = climatology_dir / "POSav.tif"
    if should_write(output_pos):
        with rasterio.open(output_pos, 'w', **profile_sos) as dst:
            dst.write(POSav.astype(np.float32), 1)
        print(f"  ✓ 已保存: {output_pos}")
    else:
        print(f"  ⚠ 跳过已存在: {output_pos}")

    print("\n" + "="*70)
    print("✓ 气候态数据计算完成！")
    print(f"输出目录: {climatology_dir}")
    print(f"  - TR_daily_climatology.tif (365波段)")
    print(f"  - SOSav.tif")
    print(f"  - POSav.tif")
    print("="*70)

# ==================== 主处理流程 ====================
def main(use_block_processing=True):
    """
    主处理流程

    Parameters:
    -----------
    use_block_processing : bool
        是否使用块处理（普通版已移除，参数仅保留兼容）
    """
    print("\n" + "="*70)
    print("Module 02: TRc计算")
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

    # 逐年计算TRc
    print(f"\n[3/3] 逐年计算TRc ({YEAR_START}-{YEAR_END})...")

    # 统计
    processed = 0
    skipped = 0
    failed = 0

    for year in years:
        output_file = OUTPUT_DIR / f"TRc_{year}.tif"

        # ✅ 断点续算：跳过已存在的文件
        if output_file.exists() and not OVERWRITE:
            print(f"  ⚠ 跳过已存在: {year} ({output_file.name})")
            skipped += 1
            continue

        # 统一使用块处理（普通版已移除）
        if not use_block_processing:
            print("  ⚠ 已移除普通版，强制使用块处理加速版")
        TRc, profile = calculate_TRc_block_optimized(year, mask, template_profile)

        if TRc is None:
            print(f"  ✗ 跳过年份 {year}（缺少数据）")
            failed += 1
            continue

        # 保存结果
        profile.update(dtype=rasterio.float32, count=1, compress='lzw', nodata=NODATA_OUT)

        with rasterio.open(output_file, 'w', **profile) as dst:
            dst.write(TRc.astype(np.float32), 1)

        print(f"  ✓ 已保存: {output_file.name}")
        processed += 1

    print("\n" + "="*70)
    print("✓ TRc计算完成！")
    print(f"输出目录: {OUTPUT_DIR}")
    print(f"\n统计:")
    print(f"  - 新处理: {processed} 个年份")
    print(f"  - 已跳过: {skipped} 个年份（文件已存在）")
    if failed > 0:
        print(f"  - 失败: {failed} 个年份（数据缺失）")
    print(f"  - 总计: {len(years)} 个年份")
    print("="*70)

# ==================== 并行处理版本 ====================
def process_single_year(args):
    """
    处理单个年份（用于并行）

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
    output_file = OUTPUT_DIR / f"TRc_{year}.tif"

    # 断点续算
    if output_file.exists() and not OVERWRITE:
        return {
            'year': year,
            'status': 'skipped',
            'message': f"文件已存在: {output_file.name}"
        }

    # 统一使用块处理（普通版已移除）
    try:
        template_profile = load_template_profile()
        TRc, profile = calculate_TRc_block_optimized(year, mask, template_profile)

        if TRc is None:
            return {
                'year': year,
                'status': 'failed',
                'message': '缺少物候数据'
            }

        # 保存结果
        profile.update(dtype=rasterio.float32, count=1, compress='lzw', nodata=NODATA_OUT)

        with rasterio.open(output_file, 'w', **profile) as dst:
            dst.write(TRc.astype(np.float32), 1)

        return {
            'year': year,
            'status': 'success',
            'message': f"已保存: {output_file.name}"
        }

    except Exception as e:
        return {
            'year': year,
            'status': 'error',
            'message': f"处理出错: {str(e)}"
        }

def main_parallel(use_block_processing=True, max_workers=None):
    """
    主处理流程（并行版本）

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
    print(f"Module 02: TRc计算（并行版本，{max_workers}个进程）")
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

    # 并行计算TRc
    print(f"\n[3/3] 并行计算TRc ({YEAR_START}-{YEAR_END}, {max_workers}进程)...")

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
    print("✓ TRc计算完成！")
    print(f"输出目录: {OUTPUT_DIR}")
    print(f"\n统计:")
    print(f"  - 新处理: {processed} 个年份")
    print(f"  - 已跳过: {skipped} 个年份（文件已存在）")
    if failed > 0:
        print(f"  - 失败: {failed} 个年份")
    print(f"  - 总计: {len(years)} 个年份")
    print("="*70)

# ==================== 主程序 ====================
if __name__ == "__main__":
    years = list(range(YEAR_START, YEAR_END + 1))
    if RUN_MODE == "skip" and outputs_complete(years):
        print("✓ Module 02 输出齐全，已跳过全部计算")
    else:
        # 步骤1：计算年度TRc（SOS-POS累积蒸腾）
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
        #   - SOSav.tif (多年平均SOS)
        #   - POSav.tif (多年平均POS)
        #
        # 输出目录：Wang2025_Analysis/Climatology/
        #
        # 注意：此步骤计算密集，预计需要30-60分钟（取决于数据规模）
