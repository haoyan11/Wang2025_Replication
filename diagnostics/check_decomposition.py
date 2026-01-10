#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
分解数据诊断工具集 (Decomposition Diagnostics)

整合了以下功能:
1. TR_fixed_window值分布检查 (check_tr_fixed_values.py的功能)
2. 物候和分解数据有效像元统计
3. Enhancement_factor计算诊断 (diagnose_enhancement_factor.py的功能)

用途:
- 验证03c固定窗口分解的正确性
- 检查TR_fixed_window是否存在负值
- 诊断enhancement_factor计算逻辑

作者: Wang2025 Replication Project
版本: 1.0.0
日期: 2025-01-05
"""

import numpy as np
from pathlib import Path
import sys
import argparse

# 导入配置
sys.path.append(str(Path(__file__).parent.parent))
from _config import (
    DECOMPOSITION_FIXED_DIR,
    TRC_ANNUAL_DIR,
    CLIMATOLOGY_DIR,
    PHENO_DIR,
    NODATA_OUT
)

# ==================== 辅助函数 ====================

def read_geotiff(filepath):
    """读取GeoTIFF（单波段）"""
    try:
        from osgeo import gdal
        gdal.UseExceptions()
        ds = gdal.Open(str(filepath))
        if ds is None:
            return None, None
        band = ds.GetRasterBand(1)
        data = band.ReadAsArray()
        nodata = band.GetNoDataValue()
        ds = None
        return data, nodata
    except Exception as e:
        print(f"❌ 读取失败: {e}")
        return None, None

def read_multiband_geotiff(filepath, band_count=365):
    """读取多波段GeoTIFF（气候态）"""
    try:
        from osgeo import gdal
        gdal.UseExceptions()
        ds = gdal.Open(str(filepath))
        if ds is None:
            return None, None
        data = np.zeros((band_count, ds.RasterYSize, ds.RasterXSize), dtype=np.float32)
        for i in range(band_count):
            band = ds.GetRasterBand(i + 1)
            data[i] = band.ReadAsArray()
        nodata = ds.GetRasterBand(1).GetNoDataValue()
        ds = None
        return data, nodata
    except Exception as e:
        print(f"❌ 读取失败: {e}")
        return None, None

def _sum_cum_range(cum_arr, start_doy, end_doy):
    """
    累积和计算（与03c中的_sum_cum_range一致）
    cum_arr: shape (365, H, W)
    start_doy, end_doy: shape (H, W) - DOY值 (1-365)
    """
    flat = cum_arr.reshape(cum_arr.shape[0], -1)
    start = start_doy.ravel()
    end = end_doy.ravel()
    idx = np.arange(start.size)
    end_vals = flat[end - 1, idx]
    start_vals = np.where(start > 1, flat[start - 2, idx], 0.0)
    return (end_vals - start_vals).reshape(start_doy.shape)

# ==================== 功能1: 检查TR_fixed_window值分布 ====================

def check_tr_values(years=[1982, 2000, 2018]):
    """
    检查TR_fixed_window的正负值分布

    目的: 验证TR_fixed_window确实存在负值（干旱年份）

    Parameters
    ----------
    years : list
        要检查的年份列表
    """
    print(f"\n{'='*70}")
    print(f"TR_fixed_window 值分布检查")
    print(f"{'='*70}\n")
    print(f"输出目录: {DECOMPOSITION_FIXED_DIR}\n")

    if not DECOMPOSITION_FIXED_DIR.exists():
        print(f"❌ 输出目录不存在: {DECOMPOSITION_FIXED_DIR}")
        return

    for year in years:
        filepath = DECOMPOSITION_FIXED_DIR / f"TR_fixed_window_{year}.tif"

        print(f"{'='*70}")
        print(f"检查文件: {filepath.name}")
        print(f"{'='*70}\n")

        if not filepath.exists():
            print(f"❌ 文件不存在: {filepath}\n")
            continue

        data, nodata = read_geotiff(filepath)
        if data is None:
            continue

        # 掩膜有效数据
        if nodata is not None:
            valid_mask = (data != nodata) & np.isfinite(data)
        else:
            valid_mask = np.isfinite(data)

        valid_data = data[valid_mask]

        if len(valid_data) == 0:
            print("⚠️  没有有效数据！\n")
            continue

        # 统计信息
        print(f"📊 数据统计:")
        print(f"  总像元数: {data.size:,}")
        print(f"  有效像元: {len(valid_data):,} ({len(valid_data)/data.size*100:.2f}%)")
        print(f"  NODATA值: {nodata}")

        print(f"\n📈 值域分布:")
        print(f"  最小值: {valid_data.min():.6f}")
        print(f"  最大值: {valid_data.max():.6f}")
        print(f"  平均值: {valid_data.mean():.6f}")
        print(f"  中位数: {np.median(valid_data):.6f}")
        print(f"  标准差: {valid_data.std():.6f}")

        # 关键检查：负值
        n_negative = np.sum(valid_data < 0)
        n_zero = np.sum(valid_data == 0)
        n_positive = np.sum(valid_data > 0)

        print(f"\n🔍 正负值分布:")
        print(f"  负值像元: {n_negative:,} ({n_negative/len(valid_data)*100:.2f}%)")
        print(f"  零值像元: {n_zero:,} ({n_zero/len(valid_data)*100:.2f}%)")
        print(f"  正值像元: {n_positive:,} ({n_positive/len(valid_data)*100:.2f}%)")

        if n_negative > 0:
            print(f"\n✅ 确认存在负值！")
            print(f"  负值范围: [{valid_data[valid_data < 0].min():.6f}, {valid_data[valid_data < 0].max():.6f}]")
            print(f"  负值均值: {valid_data[valid_data < 0].mean():.6f}")
        else:
            print(f"\n⚠️  未发现负值！这可能表明计算有问题。")

        # 百分位数
        print(f"\n📉 百分位数:")
        percentiles = [1, 5, 10, 25, 50, 75, 90, 95, 99]
        for p in percentiles:
            val = np.percentile(valid_data, p)
            print(f"  {p:2d}%: {val:10.6f}")

        print()

    print(f"{'='*70}")
    print(f"✅ 检查完成！")
    print(f"{'='*70}\n")

# ==================== 功能2: 检查物候有效像元数 ====================

def check_phenology_pixels(year=1982):
    """
    检查物候和分解数据的有效像元数

    目的: 验证SOS/POS/TR的有效像元数是否合理

    Parameters
    ----------
    year : int
        要检查的年份
    """
    print(f"\n{'='*70}")
    print(f"物候和分解数据有效像元检查 ({year}年)")
    print(f"{'='*70}\n")

    # 读取气候态
    sos_av, _ = read_geotiff(CLIMATOLOGY_DIR / "SOSav.tif")
    pos_av, _ = read_geotiff(CLIMATOLOGY_DIR / "POSav.tif")

    if sos_av is None or pos_av is None:
        print("❌ 无法读取气候态数据")
        return

    total_pixels = sos_av.size

    # 统计气候态
    sos_av_valid = np.isfinite(sos_av) & (sos_av > 0) & (sos_av <= 365)
    pos_av_valid = np.isfinite(pos_av) & (pos_av > 0) & (pos_av <= 365)

    print(f"📊 气候态物候:")
    print(f"  SOSav有效: {np.sum(sos_av_valid):,} ({np.sum(sos_av_valid)/total_pixels*100:.2f}%)")
    print(f"  POSav有效: {np.sum(pos_av_valid):,} ({np.sum(pos_av_valid)/total_pixels*100:.2f}%)")
    print(f"  两者交集: {np.sum(sos_av_valid & pos_av_valid):,} ({np.sum(sos_av_valid & pos_av_valid)/total_pixels*100:.2f}%)")

    # 读取当年数据
    sos_y, _ = read_geotiff(PHENO_DIR / f"sos_gpp_{year}.tif")
    pos_y, _ = read_geotiff(PHENO_DIR / f"pos_doy_gpp_{year}.tif")
    tr_fixed, _ = read_geotiff(DECOMPOSITION_FIXED_DIR / f"TR_fixed_window_{year}.tif")

    if sos_y is None or pos_y is None or tr_fixed is None:
        print("\n❌ 无法读取年度数据")
        return

    sos_y_valid = np.isfinite(sos_y) & (sos_y > 0) & (sos_y <= 365)
    pos_y_valid = np.isfinite(pos_y) & (pos_y > 0) & (pos_y <= 365)
    tr_valid = np.isfinite(tr_fixed) & (tr_fixed != NODATA_OUT) & (np.abs(tr_fixed) < 9000)

    print(f"\n📊 {year}年数据:")
    print(f"  SOS有效: {np.sum(sos_y_valid):,} ({np.sum(sos_y_valid)/total_pixels*100:.2f}%)")
    print(f"  POS有效: {np.sum(pos_y_valid):,} ({np.sum(pos_y_valid)/total_pixels*100:.2f}%)")
    print(f"  两者交集: {np.sum(sos_y_valid & pos_y_valid):,} ({np.sum(sos_y_valid & pos_y_valid)/total_pixels*100:.2f}%)")
    print(f"\n  TR_fixed_window有效: {np.sum(tr_valid):,} ({np.sum(tr_valid)/total_pixels*100:.2f}%)")

    # 完整交集
    complete = sos_y_valid & pos_y_valid & tr_valid
    print(f"\n  完整交集（SOS∩POS∩TR）: {np.sum(complete):,} ({np.sum(complete)/total_pixels*100:.2f}%)")

    print(f"\n{'='*70}\n")

# ==================== 功能3: 诊断enhancement_factor（高级）====================

def diagnose_enhancement_factor(row, col, year):
    """
    详细诊断单个像元的enhancement_factor计算

    目的: 手动重现03c的计算过程，验证逻辑正确性

    Parameters
    ----------
    row, col : int
        像元坐标
    year : int
        年份
    """
    print(f"\n{'='*80}")
    print(f"Enhancement Factor 诊断 - 像元[{row},{col}] 年份{year}")
    print(f"{'='*80}\n")

    # 1. 读取当年数据
    trc_y, _ = read_geotiff(TRC_ANNUAL_DIR / f"TRc_{year}.tif")
    sos_y, _ = read_geotiff(PHENO_DIR / f"sos_gpp_{year}.tif")
    pos_y, _ = read_geotiff(PHENO_DIR / f"pos_doy_gpp_{year}.tif")

    if trc_y is None or sos_y is None or pos_y is None:
        print("❌ 无法读取当年数据")
        return

    trc_y_val = trc_y[row, col]
    sos_y_val = sos_y[row, col]
    pos_y_val = pos_y[row, col]

    print(f"📊 当年数据:")
    print(f"  TRc_y = {trc_y_val:.3f} mm")
    print(f"  SOS_y = {sos_y_val:.1f} (DOY)")
    print(f"  POS_y = {pos_y_val:.1f} (DOY)")
    print(f"  窗口长度 = {pos_y_val - sos_y_val + 1:.0f} days")

    # 2. 读取气候态
    tr_daily_av, _ = read_multiband_geotiff(CLIMATOLOGY_DIR / "TR_daily_climatology.tif", 365)
    sos_av, _ = read_geotiff(CLIMATOLOGY_DIR / "SOSav.tif")
    pos_av, _ = read_geotiff(CLIMATOLOGY_DIR / "POSav.tif")

    if tr_daily_av is None or sos_av is None or pos_av is None:
        print("❌ 无法读取气候态")
        return

    sos_av_val = sos_av[row, col]
    pos_av_val = pos_av[row, col]

    print(f"\n📈 气候态:")
    print(f"  SOSav = {sos_av_val:.1f} (DOY)")
    print(f"  POSav = {pos_av_val:.1f} (DOY)")
    print(f"  固定窗口长度 = {pos_av_val - sos_av_val + 1:.0f} days")

    # 3. 计算累积气候态
    tr_daily_av[~np.isfinite(tr_daily_av)] = 0.0
    tr_daily_av[tr_daily_av < 0] = 0.0
    tr_cum = np.cumsum(tr_daily_av[:, row, col])

    print(f"\n  TR气候态累积 (采样):")
    for d in [0, 30, 60, 90, 120]:
        if d < len(tr_cum):
            print(f"    DOY {d+1:3d}: cum={tr_cum[d]:.2f} mm")

    # 4. 计算TRc_av（固定窗口的气候态累积）
    sos_av_i = int(np.clip(np.rint(sos_av_val), 1, 365))
    pos_av_i = int(np.clip(np.rint(pos_av_val), 1, 365))

    if sos_av_i > 1:
        trc_av_val = tr_cum[pos_av_i - 1] - tr_cum[sos_av_i - 2]
    else:
        trc_av_val = tr_cum[pos_av_i - 1]

    print(f"\n  TRc_av = {trc_av_val:.3f} mm")
    print(f"    计算方式: tr_cum[{pos_av_i-1}] - tr_cum[{sos_av_i-2 if sos_av_i > 1 else 'N/A'}]")
    print(f"    = {tr_cum[pos_av_i - 1]:.3f} - {tr_cum[sos_av_i - 2] if sos_av_i > 1 else 0:.3f}")

    # 5. 计算TRc_y_clim（当年窗口的气候态累积）
    sos_y_i = int(np.clip(np.rint(sos_y_val), 1, 365))
    pos_y_i = int(np.clip(np.rint(pos_y_val), 1, 365))

    if sos_y_i > 1:
        trc_y_clim_val = tr_cum[pos_y_i - 1] - tr_cum[sos_y_i - 2]
    else:
        trc_y_clim_val = tr_cum[pos_y_i - 1]

    print(f"\n  TRc_y_clim = {trc_y_clim_val:.3f} mm")
    print(f"    计算方式: tr_cum[{pos_y_i-1}] - tr_cum[{sos_y_i-2 if sos_y_i > 1 else 'N/A'}]")
    print(f"    = {tr_cum[pos_y_i - 1]:.3f} - {tr_cum[sos_y_i - 2] if sos_y_i > 1 else 0:.3f}")

    # 6. 计算enhancement_factor
    if trc_y_clim_val > 1e-6:
        enhancement_factor_val = trc_y_val / trc_y_clim_val
    else:
        enhancement_factor_val = 1.0

    print(f"\n🔍 增强因子:")
    print(f"  enhancement_factor = TRc_y / TRc_y_clim")
    print(f"                     = {trc_y_val:.3f} / {trc_y_clim_val:.3f}")
    print(f"                     = {enhancement_factor_val:.6f}")

    if enhancement_factor_val < 1.0:
        print(f"  ✅ 增强因子 < 1: 当年低于气候态（干旱）")
    elif enhancement_factor_val > 1.0:
        print(f"  ✅ 增强因子 > 1: 当年高于气候态（湿润）")
    else:
        print(f"  ⚠️  增强因子 = 1: 当年等于气候态")

    # 7. 计算TRc_y_fixed
    trc_y_fixed_val = enhancement_factor_val * trc_av_val

    print(f"\n  TRc_y_fixed = enhancement_factor × TRc_av")
    print(f"              = {enhancement_factor_val:.6f} × {trc_av_val:.3f}")
    print(f"              = {trc_y_fixed_val:.3f} mm")

    # 8. 计算TR_fixed_window
    tr_fixed_window_val = trc_y_fixed_val - trc_av_val

    print(f"\n🎯 TR_fixed_window:")
    print(f"  TR_fixed_window = TRc_y_fixed - TRc_av")
    print(f"                  = {trc_y_fixed_val:.3f} - {trc_av_val:.3f}")
    print(f"                  = {tr_fixed_window_val:.3f} mm")

    if tr_fixed_window_val < 0:
        print(f"  ✅ TR_fixed_window < 0: 固定窗口内强度低于气候态！")
    elif tr_fixed_window_val > 0:
        print(f"  ⚠️  TR_fixed_window > 0: 固定窗口内强度高于气候态")
    else:
        print(f"  ⚠️  TR_fixed_window = 0: 固定窗口内强度等于气候态")

    # 9. 验证：读取03c输出的TR_fixed_window
    tr_fixed_file = DECOMPOSITION_FIXED_DIR / f"TR_fixed_window_{year}.tif"
    tr_fixed_output, _ = read_geotiff(tr_fixed_file)

    if tr_fixed_output is not None:
        tr_fixed_output_val = tr_fixed_output[row, col]
        print(f"\n📁 03c输出文件中的值:")
        print(f"  TR_fixed_window = {tr_fixed_output_val:.3f} mm")

        diff = abs(tr_fixed_window_val - tr_fixed_output_val)
        if diff < 0.01:
            print(f"  ✅ 手动计算与文件一致！(误差: {diff:.6f})")
        else:
            print(f"  ❌ 手动计算与文件不一致！(误差: {diff:.6f})")
            print(f"  可能存在计算错误！")

    print(f"\n{'='*80}\n")

# ==================== 主函数 ====================

def main():
    """命令行入口"""
    parser = argparse.ArgumentParser(
        description="分解数据诊断工具 - 验证03c固定窗口分解的正确性",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  # 快速检查TR值分布（默认）
  python diagnostics/check_decomposition.py

  # 检查物候像元数
  python diagnostics/check_decomposition.py --mode pixels

  # 详细诊断enhancement_factor
  python diagnostics/check_decomposition.py --mode factor --pixel 500 500 --years 1990

  # 全面检查
  python diagnostics/check_decomposition.py --mode all --years 1982 2000 2018
        """
    )

    parser.add_argument('--mode',
                        choices=['values', 'pixels', 'factor', 'all'],
                        default='values',
                        help='诊断模式: values=值分布, pixels=像元统计, factor=因子诊断, all=全部')
    parser.add_argument('--years',
                        nargs='+',
                        type=int,
                        default=[1982, 2000, 2018],
                        help='检查的年份列表（默认: 1982 2000 2018）')
    parser.add_argument('--pixel',
                        nargs=2,
                        type=int,
                        metavar=('ROW', 'COL'),
                        help='诊断像元坐标 (仅factor模式需要)')

    args = parser.parse_args()

    print("\n" + "="*80)
    print("分解数据诊断工具 v1.0.0")
    print("="*80)

    if args.mode in ['values', 'all']:
        check_tr_values(args.years)

    if args.mode in ['pixels', 'all']:
        check_phenology_pixels(args.years[0])

    if args.mode == 'factor':
        if args.pixel is None:
            print("\n❌ 错误: factor模式需要指定--pixel参数")
            print("示例: python diagnostics/check_decomposition.py --mode factor --pixel 500 500\n")
            return
        diagnose_enhancement_factor(args.pixel[0], args.pixel[1], args.years[0])

    if args.mode == 'all' and args.pixel is not None:
        # 如果all模式且提供了pixel，也运行factor诊断
        print("\n" + "="*80)
        print("额外诊断: Enhancement Factor")
        print("="*80)
        diagnose_enhancement_factor(args.pixel[0], args.pixel[1], args.years[0])

if __name__ == "__main__":
    main()
