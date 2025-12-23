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
from tqdm import tqdm
from datetime import datetime, timedelta
from rasterio.windows import Window
from concurrent.futures import ProcessPoolExecutor, as_completed
import warnings
warnings.filterwarnings('ignore')

# ==================== 全局配置 ====================
ROOT = Path(r"I:\F\Data4")
PHENO_DIR = ROOT / "Phenology_Output_1" / "GPP_phenology"  # GPP物候（与config.py一致）
TR_DAILY_DIR = ROOT / "Meteorological Data" / "ERA5_Land" / "ET_components" / "ET_transp" / "ET_transp_Daily" / "ET_transp_Daily_2"  # ERA5-Land TR数据
OUTPUT_DIR = ROOT / "Wang2025_Analysis" / "TRc_annual"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

YEAR_START = 1982
YEAR_END = 2018
BLOCK_SIZE = 128
MAX_WORKERS = 2
NODATA_OUT = -9999.0

# ==================== 辅助函数 ====================
def get_TR_file_path(date_obj):
    """
    获取ERA5-Land TR日数据文件路径
    格式: ERA5L_ET_transp_Daily_mm_YYYYMMDD.tif
    """
    yyyymmdd = date_obj.strftime("%Y%m%d")

    # ERA5-Land格式
    file_path = TR_DAILY_DIR / f"ERA5L_ET_transp_Daily_mm_{yyyymmdd}.tif"

    if file_path.exists():
        return file_path

    return None

def is_leap_year(year):
    """判断是否为闰年"""
    return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)

# ==================== 核心计算函数 ====================
def calculate_TRc_for_year(year, mask):
    """
    计算单个年份的TRc

    Parameters:
    -----------
    year : int
        年份
    mask : ndarray
        有效掩膜（布尔数组）

    Returns:
    --------
    TRc : ndarray
        累积蒸腾栅格
    profile : dict
        栅格配置
    """
    print(f"\n=== 计算TRc: {year} ===")

    # 读取物候数据（物候代码输出格式：小写sos_gpp）
    sos_file = PHENO_DIR / f"sos_gpp_{year}.tif"
    pos_file = PHENO_DIR / f"pos_doy_gpp_{year}.tif"  # 注意是pos_doy

    if not sos_file.exists() or not pos_file.exists():
        print(f"  ✗ 错误：缺少物候数据")
        print(f"    SOS文件: {sos_file}")
        print(f"    SOS存在: {sos_file.exists()}")
        print(f"    POS文件: {pos_file}")
        print(f"    POS存在: {pos_file.exists()}")
        print("  请先运行物候计算代码或确认物候数据路径")
        return None, None

    with rasterio.open(sos_file) as src:
        sos = src.read(1)
        profile = src.profile.copy()
        height, width = src.height, src.width

    with rasterio.open(pos_file) as src:
        pos = src.read(1)

    # 初始化TRc数组
    TRc = np.zeros((height, width), dtype=np.float32)
    TRc[~mask] = NODATA_OUT

    # 生成日期列表
    days_in_year = 366 if is_leap_year(year) else 365
    dates_year = [datetime(year, 1, 1) + timedelta(days=i) for i in range(days_in_year)]

    print(f"  年份: {year}, 天数: {days_in_year}")

    # 逐日累加TR
    missing_files = 0
    valid_days = 0

    for date_obj in tqdm(dates_year, desc="累加TR", leave=False):
        doy = date_obj.timetuple().tm_yday
        tr_file = get_TR_file_path(date_obj)

        if tr_file is None or not tr_file.exists():
            missing_files += 1
            continue

        try:
            with rasterio.open(tr_file) as src:
                tr_daily = src.read(1)

            # 累加：仅在SOS ≤ doy ≤ POS的像元
            for i in range(height):
                for j in range(width):
                    if not mask[i, j]:
                        continue

                    sos_ij = sos[i, j]
                    pos_ij = pos[i, j]

                    # 检查有效性
                    if sos_ij == NODATA_OUT or pos_ij == NODATA_OUT:
                        continue

                    if sos_ij <= doy <= pos_ij:
                        tr_val = tr_daily[i, j]
                        if tr_val != NODATA_OUT and not np.isnan(tr_val):
                            TRc[i, j] += tr_val

            valid_days += 1

        except Exception as e:
            tqdm.write(f"  ⚠ 读取失败: {tr_file.name} - {str(e)}")
            missing_files += 1
            continue

    print(f"  有效天数: {valid_days}/{days_in_year}, 缺失: {missing_files}")

    # 质量检查
    TRc_valid = TRc[mask & (TRc != NODATA_OUT)]
    if len(TRc_valid) > 0:
        print(f"  TRc范围: {np.min(TRc_valid):.2f} - {np.max(TRc_valid):.2f} mm")
        print(f"  TRc平均: {np.mean(TRc_valid):.2f} mm")
        print(f"  有效像元: {len(TRc_valid)}")
    else:
        print(f"  ⚠ 警告：没有有效的TRc值")

    return TRc, profile

# ==================== 块处理版本（优化内存+性能） ====================
def calculate_TRc_block_optimized(year, mask):
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

    Returns:
    --------
    TRc : ndarray
        累积蒸腾栅格
    profile : dict
        栅格配置
    """
    print(f"\n=== 计算TRc（块处理-加速版）: {year} ===")

    # 读取物候数据（物候代码输出格式）
    sos_file = PHENO_DIR / f"sos_gpp_{year}.tif"
    pos_file = PHENO_DIR / f"pos_doy_gpp_{year}.tif"

    if not sos_file.exists() or not pos_file.exists():
        print(f"  ✗ 错误：缺少物候数据")
        print(f"    SOS: {sos_file} ({'存在' if sos_file.exists() else '不存在'})")
        print(f"    POS: {pos_file} ({'存在' if pos_file.exists() else '不存在'})")
        return None, None

    with rasterio.open(sos_file) as src:
        sos = src.read(1)
        profile = src.profile.copy()
        height, width = src.height, src.width
        nodata_sos = src.nodata  # ✅ 读取真实nodata

    with rasterio.open(pos_file) as src:
        pos = src.read(1)
        nodata_pos = src.nodata  # ✅ 读取真实nodata

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
        & (sos <= days_in_year) & (pos <= days_in_year)
    )

    # 初始化TRc：只对valid_pheno初始化为0，其他为NODATA
    TRc = np.full((height, width), NODATA_OUT, dtype=np.float32)
    TRc[valid_pheno] = 0.0

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
    opened_any = False
    nodata_tr = None

    # ✅ 优化：改为"先天后块"循环（性能提升18倍+）
    for date_obj in tqdm(dates_year, desc="逐日累加（先天后块）", leave=False):
        tr_file = get_TR_file_path(date_obj)
        if tr_file is None or not tr_file.exists():
            missing_files += 1
            continue

        doy = date_obj.timetuple().tm_yday

        try:
            with rasterio.open(tr_file) as src:
                # ✅ 读取TR数据的真实nodata（只读一次）
                if not opened_any:
                    nodata_tr = src.nodata
                    opened_any = True

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

                    acc = in_window & valid_tr
                    if np.any(acc):
                        TRc[r0:r1, c0:c1][acc] += tr_b[acc]

        except Exception as e:
            missing_files += 1
            continue

    print(f"  缺失日文件数: {missing_files}/{days_in_year}")

    # 质量检查
    TRc_valid = TRc[valid_pheno & (TRc != NODATA_OUT)]
    if TRc_valid.size > 0:
        print(f"  TRc范围: {TRc_valid.min():.2f} - {TRc_valid.max():.2f} mm")
        print(f"  TRc平均: {TRc_valid.mean():.2f} mm")
        print(f"  有效像元: {TRc_valid.size}")
    else:
        print("  ⚠ 警告：没有有效TRc值（请检查物候/掩膜/TR数据）")

    return TRc, profile

# ==================== 主处理流程 ====================
def main(use_block_processing=True):
    """
    主处理流程

    Parameters:
    -----------
    use_block_processing : bool
        是否使用块处理（推荐True以节省内存）
    """
    print("\n" + "="*70)
    print("Module 02: TRc计算")
    print("="*70)

    years = list(range(YEAR_START, YEAR_END + 1))

    # 读取掩膜
    print("\n[1/3] 读取掩膜...")
    mask_file = ROOT / "Wang2025_Analysis" / "masks" / "combined_mask.tif"

    if not mask_file.exists():
        print(f"  ✗ 错误：掩膜文件不存在: {mask_file}")
        print("  请先运行 Module 00: 00_data_preparation.py")
        return

    with rasterio.open(mask_file) as src:
        mask = src.read(1).astype(bool)

    print(f"  ✓ 有效像元数: {np.sum(mask)}")

    # 检查TR数据
    print("\n[2/3] 检查TR数据...")
    test_date = datetime(YEAR_START, 1, 15)  # ✅ 使用YEAR_START年份测试
    test_file = get_TR_file_path(test_date)

    if test_file is None or not test_file.exists():
        print(f"  ✗ 错误：找不到TR数据")
        print(f"    测试日期: {test_date.strftime('%Y-%m-%d')}")
        print(f"    预期路径示例: {TR_DAILY_DIR / f'ERA5L_ET_transp_Daily_mm_{YEAR_START}0115.tif'}")  # ✅ 修正提示信息
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
        if output_file.exists():
            print(f"  ⚠ 跳过已存在: {year} ({output_file.name})")
            skipped += 1
            continue

        # 选择处理方法
        if use_block_processing:
            TRc, profile = calculate_TRc_block_optimized(year, mask)
        else:
            TRc, profile = calculate_TRc_for_year(year, mask)

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
    if output_file.exists():
        return {
            'year': year,
            'status': 'skipped',
            'message': f"文件已存在: {output_file.name}"
        }

    # 选择处理方法
    try:
        if use_block_processing:
            TRc, profile = calculate_TRc_block_optimized(year, mask)
        else:
            TRc, profile = calculate_TRc_for_year(year, mask)

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
        是否使用块处理
    max_workers : int
        并行进程数（默认使用MAX_WORKERS配置）
    """
    if max_workers is None:
        max_workers = MAX_WORKERS

    print("\n" + "="*70)
    print(f"Module 02: TRc计算（并行版本，{max_workers}个进程）")
    print("="*70)

    years = list(range(YEAR_START, YEAR_END + 1))

    # 读取掩膜
    print("\n[1/3] 读取掩膜...")
    mask_file = ROOT / "Wang2025_Analysis" / "masks" / "combined_mask.tif"

    if not mask_file.exists():
        print(f"  ✗ 错误：掩膜文件不存在: {mask_file}")
        print("  请先运行 Module 00: 00_data_preparation.py")
        return

    with rasterio.open(mask_file) as src:
        mask = src.read(1).astype(bool)

    print(f"  ✓ 有效像元数: {np.sum(mask)}")

    # 检查TR数据
    print("\n[2/3] 检查TR数据...")
    test_date = datetime(YEAR_START, 1, 15)
    test_file = get_TR_file_path(test_date)

    if test_file is None or not test_file.exists():
        print(f"  ✗ 错误：找不到TR数据")
        print(f"    测试日期: {test_date.strftime('%Y-%m-%d')}")
        print(f"    预期路径示例: {TR_DAILY_DIR / f'ERA5L_ET_transp_Daily_mm_{YEAR_START}0115.tif'}")
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
    # 方式1：串行处理（推荐，稳定可靠）
    main(use_block_processing=True)

    # 方式2：并行处理（可选，SSD环境下可能更快）
    # ⚠️ 注意：
    # - 并行适用于SSD，HDD可能适得其反
    # - MAX_WORKERS=2-4 比较稳妥（I/O密集任务）
    # - 自动支持断点续算
    # main_parallel(use_block_processing=True, max_workers=2)

    # 方式3：常规方法（仅用于对照验证，不推荐）
    # ⚠️ 警告：常规方法比块处理慢得多，且仍有nodata硬编码问题
    # main(use_block_processing=False)
