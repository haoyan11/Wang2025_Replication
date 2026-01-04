#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
一次性脚本: 将物候数据从Clarke 1866重投影到EPSG:4326
运行此脚本后，所有物候文件将被原地替换为EPSG:4326版本
"""

import numpy as np
import rasterio
from rasterio.warp import calculate_default_transform, reproject
from rasterio.enums import Resampling
from pathlib import Path
from tqdm import tqdm
import shutil
import sys

# 导入配置
from _config import PHENO_DIR, PHENO_FILE_FORMAT, YEAR_START, YEAR_END, ROOT

# 目标CRS
TARGET_CRS = 'EPSG:4326'

# 掩膜文件路径（用于裁剪重投影后的数据）
MASK_FILE = ROOT / "Wang2025_Analysis" / "masks" / "combined_mask.tif"

# 输出目录（保留原始数据，重投影到新目录）
OUTPUT_SUFFIX = "_EPSG4326"

# =========================
# 固定输入/输出路径（仅新增这一段）
# 输入：原始 Clarke 1866 目录
INPUT_PHENO_DIR = ROOT / "Phenology_Output_1" / "GPP_phenology"
# 输出：EPSG:4326 目录
OUTPUT_PHENO_DIR = ROOT / "Phenology_Output_1" / "GPP_phenology_EPSG4326"
# =========================

def reproject_file_to_new_dir(file_path, output_dir, mask_array=None, mask_transform=None):
    """
    将单个文件重投影到EPSG:4326并应用掩膜（保存到新目录）

    Parameters:
    -----------
    file_path : Path
        要重投影的文件路径
    output_dir : Path
        输出目录路径
    mask_array : ndarray, optional
        掩膜数组（布尔类型），用于过滤out-of-mask像元
    mask_transform : Affine, optional
        掩膜的transform（用于验证对齐）

    Returns:
    --------
    success : bool
        是否成功重投影
    """
    if not file_path.exists():
        print(f"  [SKIP] 文件不存在: {file_path.name}")
        return False

    # 输出文件路径
    output_file = output_dir / file_path.name

    try:
        with rasterio.open(file_path) as src:
            src_crs = src.crs

            # 检查是否已经是EPSG:4326
            if src_crs is not None and str(src_crs) == TARGET_CRS:
                print(f"  [SKIP] 已是{TARGET_CRS}: {file_path.name}")
                # 如果已经是目标CRS，直接复制到输出目录
                if not output_file.exists():
                    shutil.copy2(file_path, output_file)
                return True

            # 计算变换参数
            transform, width, height = calculate_default_transform(
                src_crs, TARGET_CRS, src.width, src.height, *src.bounds
            )

            # 读取源数据
            src_data = src.read(1)
            nodata = src.nodata

            # 创建目标数组
            dst_data = np.empty((height, width), dtype=src_data.dtype)

            # 重投影
            reproject(
                source=src_data,
                destination=dst_data,
                src_transform=src.transform,
                src_crs=src_crs,
                src_nodata=nodata,
                dst_transform=transform,
                dst_crs=TARGET_CRS,
                dst_nodata=nodata,
                resampling=Resampling.nearest  # 物候数据使用最近邻（整数DOY）
            )

            # 应用掩膜（如果提供）
            if mask_array is not None:
                # 验证尺寸和transform匹配
                if mask_array.shape != dst_data.shape:
                    print(f"  [WARN] 掩膜尺寸不匹配: mask {mask_array.shape} vs data {dst_data.shape}")
                elif mask_transform != transform:
                    print(f"  [WARN] 掩膜transform不匹配")
                else:
                    # 应用掩膜：将mask为False的像元设为nodata
                    if nodata is None:
                        nodata = -9999  # 如果原文件没有nodata，设置一个

                    out_of_mask_count = np.sum(~mask_array & (dst_data != nodata) & np.isfinite(dst_data) & (dst_data > 0))
                    dst_data[~mask_array] = nodata

                    if out_of_mask_count > 0:
                        print(f"  [MASK] 已移除 {out_of_mask_count} 个out-of-mask像元")

            # 更新profile
            profile = src.profile.copy()
            profile.update({
                'crs': TARGET_CRS,
                'transform': transform,
                'width': width,
                'height': height,
                'nodata': nodata
            })

        # 保存到输出目录（不修改原文件）
        with rasterio.open(output_file, 'w', **profile) as dst:
            dst.write(dst_data, 1)

        print(f"  [OK] 重投影成功: {file_path.name} ({src_crs} -> {TARGET_CRS})")
        return True

    except Exception as e:
        print(f"  [FAIL] 重投影失败: {file_path.name}")
        print(f"        错误: {e}")
        return False

def main():
    """批量重投影所有物候文件"""
    print("\n" + "="*70)
    print("物候数据CRS统一工具: Clarke 1866 -> EPSG:4326 + 掩膜过滤")
    print("="*70)

    # 检查物候目录
    if not INPUT_PHENO_DIR.exists():
        print(f"\n[ERROR] 物候目录不存在: {INPUT_PHENO_DIR}")
        return

    # 加载combined_mask
    mask_array = None
    mask_transform = None

    if MASK_FILE.exists():
        print(f"\n加载掩膜文件: {MASK_FILE.name}")
        try:
            with rasterio.open(MASK_FILE) as mask_src:
                mask_data = mask_src.read(1)
                mask_array = mask_data.astype(bool)  # 转换为布尔数组
                mask_transform = mask_src.transform
                mask_crs = mask_src.crs

                # 验证掩膜CRS是否为目标CRS
                if str(mask_crs) != TARGET_CRS:
                    print(f"  [WARN] 掩膜CRS ({mask_crs}) 与目标CRS ({TARGET_CRS}) 不匹配")
                    print(f"         将跳过掩膜过滤")
                    mask_array = None
                else:
                    print(f"  [OK] 掩膜已加载: {np.sum(mask_array)} 有效像元")
        except Exception as e:
            print(f"  [WARN] 掩膜加载失败: {e}")
            print(f"         将继续但不应用掩膜过滤")
    else:
        print(f"\n[INFO] 掩膜文件不存在: {MASK_FILE}")
        print(f"       将重投影但不应用掩膜过滤")

    # 创建输出目录
    output_dir = OUTPUT_PHENO_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n源目录: {INPUT_PHENO_DIR}")
    print(f"输出目录: {output_dir}")
    print(f"年份范围: {YEAR_START}-{YEAR_END}")
    print(f"目标CRS: {TARGET_CRS}")

    # 收集所有需要处理的文件
    files_to_process = []

    # 添加年度物候文件
    pheno_vars = ['SOS', 'POS', 'EOS']  # 只处理核心物候指标
    for year in range(YEAR_START, YEAR_END + 1):
        for var in pheno_vars:
            if var in PHENO_FILE_FORMAT:
                filename = PHENO_FILE_FORMAT[var].format(year=year)
                file_path = INPUT_PHENO_DIR / filename
                files_to_process.append(file_path)

    print(f"\n待处理文件数: {len(files_to_process)}")

    # 询问用户确认
    print("\n[INFO] 此操作将保留原始文件，重投影结果保存到:")
    print(f"       {output_dir}")

    # 检查是否有--yes参数
    if '--yes' not in sys.argv:
        response = input("\n是否继续？(yes/no): ").strip().lower()
        if response not in ['yes', 'y']:
            print("\n[CANCELLED] 操作已取消")
            return
    else:
        print("\n[AUTO] 自动确认 (--yes参数)")

    # 批量处理
    print(f"\n{'='*70}")
    print("开始批量重投影...")
    print(f"{'='*70}\n")

    success_count = 0
    skip_count = 0
    fail_count = 0

    for file_path in tqdm(files_to_process, desc="重投影+掩膜进度"):
        if not file_path.exists():
            skip_count += 1
            continue

        if reproject_file_to_new_dir(file_path, output_dir, mask_array, mask_transform):
            success_count += 1
        else:
            fail_count += 1

    # 输出统计
    print(f"\n{'='*70}")
    print("处理完成")
    print(f"{'='*70}")
    print(f"  [OK] 成功重投影: {success_count} 文件")
    print(f"  [SKIP] 跳过（不存在）: {skip_count} 文件")
    print(f"  [FAIL] 失败: {fail_count} 文件")
    print(f"  总计: {len(files_to_process)} 文件")

    if fail_count > 0:
        print(f"\n[WARNING] 有 {fail_count} 个文件处理失败，请检查错误信息")
    elif success_count > 0:
        print(f"\n[SUCCESS] 所有物候数据已统一为 {TARGET_CRS}")
    else:
        print(f"\n[INFO] 所有文件已是 {TARGET_CRS}，已复制到输出目录")

if __name__ == "__main__":
    main()
