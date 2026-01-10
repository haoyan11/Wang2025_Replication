#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module 00: 数据准备
创建分析所需的掩膜文件、重投影物候并检查数据完整性
1. 创建纬度掩膜（≥30°N）
2. [可选] 创建森林类型掩膜（IGBP分类1-5）
3. 合并掩膜
4. 物候数据重投影到EPSG:4326（可选）
5. 检查数据完整性
"""

import numpy as np
import rasterio
from rasterio.warp import reproject
from rasterio.enums import Resampling
from pathlib import Path
from datetime import datetime, timedelta
from tqdm import tqdm
import shutil
import sys
import warnings
warnings.filterwarnings('ignore')

# 导入配置（固定使用_config.py）
from _config import (
    ROOT, LAT_MIN, FOREST_CLASSES, NODATA_OUT,
    USE_FOREST_MASK, LANDCOVER_FILE, TR_DAILY_DIR, PHENO_DIR,
    YEAR_START, YEAR_END, PHENO_FILE_FORMAT,
    GPP_DAILY_DIR, SM_DAILY_DIR, TR_FILE_FORMAT, get_TR_file_path,
    TEMPLATE_RASTER, MASK_FILE, get_GPP_file_path
)

# 运行开关
RUN_REPROJECT_PHENO = True
RUN_VERIFY_DATA = True
RUN_VERIFY_STRICT = False  # True时验证失败直接退出
REPROJECT_AUTO_YES = False  # True时跳过交互确认

# 重投影配置
TARGET_CRS = 'EPSG:4326'
PHENO_DIR_NAME = PHENO_DIR.name
if PHENO_DIR_NAME.endswith("_EPSG4326"):
    INPUT_PHENO_DIR = PHENO_DIR.parent / PHENO_DIR_NAME.replace("_EPSG4326", "")
    OUTPUT_PHENO_DIR = PHENO_DIR
else:
    INPUT_PHENO_DIR = PHENO_DIR
    OUTPUT_PHENO_DIR = PHENO_DIR.parent / f"{PHENO_DIR_NAME}_EPSG4326"

# 输出目录
OUTPUT_DIR = ROOT / "Wang2025_Analysis" / "masks"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 交集掩膜配置（按首个有效文件做交集）
INTERSECT_DATA_MASK = True

# ==================== 辅助函数 ====================
def ensure_template_raster():
    """
    确保模板栅格存在（全流程只认 TEMPLATE_RASTER）
    若缺失，则从TR日尺度文件复制网格信息生成空模板。
    """
    if TEMPLATE_RASTER.exists():
        return TEMPLATE_RASTER
    candidate_years = [YEAR_START, (YEAR_START + YEAR_END) // 2, YEAR_END]
    candidate_dates = [(1, 15), (7, 1), (10, 1)]
    src_file = None
    for year in candidate_years:
        for month, day in candidate_dates:
            test_date = datetime(year, month, day)
            tr_file = get_TR_file_path(test_date)
            if tr_file and tr_file.exists():
                src_file = tr_file
                break
        if src_file:
            break
    if src_file is None:
        raise FileNotFoundError("找不到可用的TR日尺度文件来生成模板栅格")
    TEMPLATE_RASTER.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(src_file) as src:
        profile = src.profile.copy()
        data = np.full((profile["height"], profile["width"]), NODATA_OUT, dtype=np.float32)
        profile.update(dtype=rasterio.float32, count=1, compress='lzw', nodata=NODATA_OUT)
    with rasterio.open(TEMPLATE_RASTER, "w", **profile) as dst:
        dst.write(data, 1)
    print(f"  ✓ 已创建模板栅格: {TEMPLATE_RASTER}")
    return TEMPLATE_RASTER


def get_template_file():
    """获取模板文件（统一网格）"""
    return ensure_template_raster()


def load_template_profile(template_file):
    with rasterio.open(template_file) as src:
        profile = src.profile.copy()
    return profile


def _compare_profile(ref_profile, other_profile, label, file_path):
    if ref_profile["height"] != other_profile["height"] or ref_profile["width"] != other_profile["width"]:
        raise ValueError(
            f"{label}尺寸不一致: {file_path}\n"
            f"  Expected: {ref_profile['height']}x{ref_profile['width']}\n"
            f"  Got: {other_profile['height']}x{other_profile['width']}"
        )
    if ref_profile.get("crs") != other_profile.get("crs"):
        raise ValueError(
            f"{label} CRS不一致: {file_path}\n"
            f"  Expected: {ref_profile.get('crs')}\n"
            f"  Got: {other_profile.get('crs')}"
        )
    ref_transform = ref_profile.get("transform")
    other_transform = other_profile.get("transform")
    if ref_transform is not None and other_transform is not None:
        transform_match = all(
            abs(ref_transform[i] - other_transform[i]) < 1e-6
            for i in range(6)
        )
        if not transform_match:
            raise ValueError(
                f"{label} Transform不一致: {file_path}\n"
                f"  Expected: {ref_transform}\n"
                f"  Got: {other_transform}"
            )


def _check_profile_match(ref_profile, file_path, label):
    if not file_path.exists():
        print(f"  ⚠ 未找到{label}样本: {file_path.name}")
        return False
    with rasterio.open(file_path) as src:
        _compare_profile(ref_profile, src.profile, label, file_path)
    print(f"  ✓ {label}一致: {file_path.name}")
    return True


def fast_grid_consistency_check(template_profile):
    """
    快速一致性检查（fail-fast）：
    mask/SOS/POS/GPP/SM/土地覆盖 与模板网格一致性
    """
    print("\n[0] 快速一致性检查（fail-fast）:")

    candidate_years = [YEAR_START, (YEAR_START + YEAR_END) // 2, YEAR_END]

    # 物候样本（SOS/POS）
    for var_key, label in [("SOS", "SOS物候"), ("POS", "POS物候")]:
        sample = None
        for year in candidate_years:
            pheno_file = PHENO_DIR / PHENO_FILE_FORMAT[var_key].format(year=year)
            if pheno_file.exists():
                sample = pheno_file
                break
        if sample is None:
            print(f"  ⚠ 未找到{label}样本，跳过一致性检查")
        else:
            _check_profile_match(template_profile, sample, label)

    # GPP日尺度样本
    gpp_dates = [datetime(1982, 1, 15), datetime(2000, 7, 1), datetime(2018, 10, 1)]
    gpp_sample = None
    for date_obj in gpp_dates:
        cand = GPP_DAILY_DIR / f"GPP_{date_obj.strftime('%Y%m%d')}.tif"
        if cand.exists():
            gpp_sample = cand
            break
    if gpp_sample is None:
        print("  ⚠ 未找到GPP样本，跳过一致性检查")
    else:
        _check_profile_match(template_profile, gpp_sample, "GPP日尺度")

    # SM日尺度样本
    sm_dates = [datetime(1982, 1, 15), datetime(2000, 6, 15), datetime(2018, 12, 15)]
    sm_sample = None
    for date_obj in sm_dates:
        cand = SM_DAILY_DIR / f"SMrz_{date_obj.strftime('%Y%m%d')}.tif"
        if cand.exists():
            sm_sample = cand
            break
    if sm_sample is None:
        print("  ⚠ 未找到SM样本，跳过一致性检查")
    else:
        _check_profile_match(template_profile, sm_sample, "SM日尺度")

    # 土地覆盖样本
    if USE_FOREST_MASK and LANDCOVER_FILE.exists():
        _check_profile_match(template_profile, LANDCOVER_FILE, "土地覆盖")
    elif USE_FOREST_MASK:
        print(f"  ⚠ 未找到土地覆盖文件: {LANDCOVER_FILE.name}")


def reproject_to_template(src_data, src_profile, template_profile, nodata_out):
    """将栅格重投影/重采样到模板网格"""
    dst_height = template_profile["height"]
    dst_width = template_profile["width"]
    dst_transform = template_profile["transform"]
    dst_crs = template_profile["crs"]

    dst_data = np.empty((dst_height, dst_width), dtype=src_data.dtype)
    reproject(
        source=src_data,
        destination=dst_data,
        src_transform=src_profile["transform"],
        src_crs=src_profile["crs"],
        src_nodata=src_profile.get("nodata"),
        dst_transform=dst_transform,
        dst_crs=dst_crs,
        dst_nodata=nodata_out,
        resampling=Resampling.nearest
    )
    return dst_data

# ==================== 1. 创建纬度掩膜 ====================
def create_lat_mask(template_file, lat_min=30.0):
    """
    创建纬度掩膜

    Parameters:
    -----------
    template_file : Path
        模板文件路径
    lat_min : float
        最小纬度（度）

    Returns:
    --------
    mask : ndarray
        纬度掩膜（布尔数组）
    profile : dict
        栅格配置
    """
    print(f"\n[1] 创建纬度掩膜（≥{lat_min}°N）...")

    with rasterio.open(template_file) as src:
        profile = src.profile.copy()
        height, width = src.height, src.width
        transform = src.transform
        crs = src.crs

    # CRS检查：确保是地理坐标系（EPSG:4326或类似）
    if crs is not None:
        if not crs.is_geographic:
            print(f"  ⚠ 警告：模板CRS不是地理坐标系: {crs}")
            print(f"  纬度掩膜可能不准确，建议使用EPSG:4326数据")
        elif str(crs) != 'EPSG:4326':
            print(f"  ⚠ 注意：模板CRS为 {crs}，非标准EPSG:4326")
            print(f"  纬度计算可能略有偏差")

    # 生成纬度数组
    # 对于每一行，计算中心点纬度
    lat_array = np.zeros((height, width), dtype=np.float32)

    for row in range(height):
        # 计算该行的地理纬度
        lon, lat = rasterio.transform.xy(transform, row, 0, offset='center')
        lat_array[row, :] = lat

    # 创建掩膜
    mask = lat_array >= lat_min

    print(f"  栅格尺寸: {height} × {width}")
    print(f"  纬度范围: {np.min(lat_array):.2f}° - {np.max(lat_array):.2f}°")
    print(f"  有效像元数: {np.sum(mask)} ({np.sum(mask) / (height * width) * 100:.1f}%)")

    return mask.astype(np.uint8), profile

# ==================== 2. 创建森林掩膜 ====================
def create_forest_mask(landcover_file, template_profile, forest_classes=[1, 2, 3, 4, 5]):
    """
    创建森林类型掩膜

    Parameters:
    -----------
    landcover_file : Path
        MODIS土地覆盖文件路径
    forest_classes : list
        森林类型代码列表
        1: Evergreen Needleleaf Forest (ENF)
        2: Evergreen Broadleaf Forest (EBF)
        3: Deciduous Needleleaf Forest (DNF)
        4: Deciduous Broadleaf Forest (DBF)
        5: Mixed Forest (MF)

    Returns:
    --------
    mask : ndarray
        森林掩膜（布尔数组）
    profile : dict
        栅格配置
    """
    print(f"\n[2] 创建森林类型掩膜...")

    if not landcover_file.exists():
        print(f"  ⚠ 警告：土地覆盖文件不存在: {landcover_file}")
        print("  将创建全1掩膜（跳过森林类型筛选）")

        # 使用模板创建全1掩膜
        template = get_template_file()
        with rasterio.open(template) as src:
            profile = src.profile.copy()
            height, width = src.height, src.width

        mask = np.ones((height, width), dtype=np.uint8)
        return mask, profile

    with rasterio.open(landcover_file) as src:
        lc_data = src.read(1)
        profile = src.profile.copy()

    if profile["crs"] != template_profile["crs"] or profile["transform"] != template_profile["transform"] or \
       profile["width"] != template_profile["width"] or profile["height"] != template_profile["height"]:
        lc_data = reproject_to_template(lc_data, profile, template_profile, nodata_out=NODATA_OUT)
        profile = template_profile.copy()

    # 创建森林掩膜
    mask = np.isin(lc_data, forest_classes)

    print(f"  森林类型: {forest_classes}")
    print(f"  有效像元数: {np.sum(mask)} ({np.sum(mask) / mask.size * 100:.1f}%)")

    # 统计各森林类型
    for fc in forest_classes:
        count = np.sum(lc_data == fc)
        fc_names = {
            1: 'ENF (常绿针叶林)',
            2: 'EBF (常绿阔叶林)',
            3: 'DNF (落叶针叶林)',
            4: 'DBF (落叶阔叶林)',
            5: 'MF (混交林)'
        }
        print(f"    {fc_names.get(fc, f'类型{fc}')}: {count} 像元")

    return mask.astype(np.uint8), profile

# ==================== 3. 合并掩膜 ====================
def combine_masks(lat_mask, forest_mask):
    """
    合并纬度掩膜和森林掩膜

    Parameters:
    -----------
    lat_mask : ndarray
        纬度掩膜
    forest_mask : ndarray
        森林掩膜

    Returns:
    --------
    combined_mask : ndarray
        组合掩膜
    """
    print(f"\n[3] 合并掩膜...")

    # 确保尺寸一致
    if lat_mask.shape != forest_mask.shape:
        raise ValueError(
            f"掩膜尺寸不一致: 纬度掩膜 {lat_mask.shape}, 森林掩膜 {forest_mask.shape}"
        )

    # 合并（逻辑与）
    combined_mask = (lat_mask.astype(bool) & forest_mask.astype(bool)).astype(np.uint8)

    print(f"  纬度掩膜有效像元: {np.sum(lat_mask)}")
    print(f"  森林掩膜有效像元: {np.sum(forest_mask)}")
    print(f"  组合掩膜有效像元: {np.sum(combined_mask)}")
    print(f"  保留比例: {np.sum(combined_mask) / np.sum(lat_mask) * 100:.1f}%")

    return combined_mask


def _find_first_tr_file():
    # 从YEAR_START年起，寻找最早可用的TR日文件
    for day_offset in range(0, 31):
        date_obj = datetime(YEAR_START, 1, 1) + timedelta(days=day_offset)
        tr_file = get_TR_file_path(date_obj)
        if tr_file and tr_file.exists():
            return tr_file
    return None


def _find_first_gpp_file():
    # 从YEAR_START年起，寻找最早可用的GPP日文件
    for day_offset in range(0, 31):
        date_obj = datetime(YEAR_START, 1, 1) + timedelta(days=day_offset)
        gpp_file = get_GPP_file_path(date_obj, daily=True)
        if gpp_file and gpp_file.exists():
            return gpp_file
    return None


def _valid_mask_from_file(file_path, template_profile, label,
                          value_range=None, non_negative=False):
    if file_path is None or not file_path.exists():
        raise FileNotFoundError(f"{label}文件不存在: {file_path}")
    with rasterio.open(file_path) as src:
        data = src.read(1)
        profile = src.profile.copy()
        nodata = src.nodata
    _compare_profile(template_profile, profile, label, file_path)
    valid = np.isfinite(data)
    if nodata is not None and not np.isnan(nodata):
        valid &= (data != nodata)
    else:
        valid &= (data > -9000)
    if non_negative:
        valid &= (data >= 0)
    if value_range is not None:
        vmin, vmax = value_range
        valid &= (data >= vmin) & (data <= vmax)
    print(f"  {label}有效像元: {np.sum(valid)}")
    return valid


def _apply_data_intersection_mask(combined_mask, template_profile):
    if not INTERSECT_DATA_MASK:
        return combined_mask

    print("\n  [DATA] 交集掩膜：首个有效文件")
    sos_file = PHENO_DIR / PHENO_FILE_FORMAT["SOS"].format(year=YEAR_START)
    pos_file = PHENO_DIR / PHENO_FILE_FORMAT["POS"].format(year=YEAR_START)
    tr_file = _find_first_tr_file()
    gpp_file = _find_first_gpp_file()

    if tr_file is None:
        raise FileNotFoundError("未找到可用TR日尺度文件用于掩膜交集")
    if gpp_file is None:
        raise FileNotFoundError("未找到可用GPP日尺度文件用于掩膜交集")
    if not sos_file.exists():
        raise FileNotFoundError(f"未找到SOS文件: {sos_file}")
    if not pos_file.exists():
        raise FileNotFoundError(f"未找到POS文件: {pos_file}")

    sos_valid = _valid_mask_from_file(sos_file, template_profile, "SOS",
                                      value_range=(1, 365))
    pos_valid = _valid_mask_from_file(pos_file, template_profile, "POS",
                                      value_range=(1, 365))
    tr_valid = _valid_mask_from_file(tr_file, template_profile, "TR",
                                     non_negative=True)
    gpp_valid = _valid_mask_from_file(gpp_file, template_profile, "GPP",
                                      non_negative=True)

    before = np.sum(combined_mask)
    combined_mask = combined_mask & sos_valid & pos_valid & tr_valid & gpp_valid
    after = np.sum(combined_mask)
    print(f"  数据交集: {before} -> {after} (移除 {before - after})")
    return combined_mask

# ==================== 4. 数据完整性检查 ====================
def check_data_availability(years=[2000, 2010, 2018]):
    """
    检查关键数据的可用性

    Parameters:
    -----------
    years : list
        检查的年份列表
    """
    print(f"\n[4] 数据完整性检查...")

    results = {
        'pheno': {'available': 0, 'missing': 0},
        'tr': {'available': 0, 'missing': 0}
    }

    # 检查物候数据
    print("\n  物候数据 (SOS/POS/EOS):")
    for year in years:
        # 根据配置的格式检查物候文件
        sos_file = PHENO_DIR / PHENO_FILE_FORMAT['SOS'].format(year=year)
        pos_file = PHENO_DIR / PHENO_FILE_FORMAT['POS'].format(year=year)
        eos_file = PHENO_DIR / PHENO_FILE_FORMAT['EOS'].format(year=year)

        files_exist = [
            sos_file.exists(),
            pos_file.exists(),
            eos_file.exists()
        ]

        if all(files_exist):
            print(f"    {year}: ✓")
            results['pheno']['available'] += 1
        else:
            # 修复警告：先计算缺失列表再插入f-string
            missing_types = [name for i, name in enumerate(['SOS', 'POS', 'EOS']) if not files_exist[i]]
            print(f"    {year}: ✗ (缺失: {missing_types})")
            results['pheno']['missing'] += 1

    # 检查TR数据（抽样检查每年1月15日）
    print("\n  TR数据 (每年1月15日抽样 - ERA5-Land):")
    for year in years:
        test_date = datetime(year, 1, 15)
        tr_file = get_TR_file_path(test_date)

        if tr_file and tr_file.exists():
            print(f"    {year}: ✓ ({tr_file.name})")
            results['tr']['available'] += 1
        else:
            expected_name = TR_FILE_FORMAT.format(date=f"{year}0115")
            print(f"    {year}: ✗ (缺失: {expected_name})")
            results['tr']['missing'] += 1

    # 总结
    print("\n  总结:")
    print(f"    物候数据: {results['pheno']['available']}/{len(years)} 年份可用")
    print(f"    TR数据:   {results['tr']['available']}/{len(years)} 年份可用")

    if results['pheno']['missing'] > 0:
        print(f"\n  ⚠ 提示：缺少物候数据，需运行 Module 01: 01_phenology_extraction.py")

    if results['tr']['missing'] > 0:
        print(f"\n  ⚠ 提示：缺少TR数据，请检查路径: {TR_DAILY_DIR}")

    return results

# ==================== 5. 物候重投影（整合 01_reproject_phenology.py） ====================
def reproject_file_to_new_dir(file_path, output_dir, template_profile,
                              mask_array=None):
    """
    将单个文件重投影到EPSG:4326并应用掩膜（保存到新目录）
    """
    if not file_path.exists():
        print(f"  [SKIP] 文件不存在: {file_path.name}")
        return False

    output_file = output_dir / file_path.name

    try:
        with rasterio.open(file_path) as src:
            src_crs = src.crs

            # 已是目标CRS仍需重采样到模板网格并应用掩膜
            if src_crs is not None and str(src_crs) == TARGET_CRS:
                print(f"  [INFO] 已是{TARGET_CRS}: {file_path.name}（将对齐网格并套掩膜）")

            src_data = src.read(1)
            nodata = src.nodata
            dst_data = np.empty((template_profile["height"], template_profile["width"]),
                                dtype=src_data.dtype)

            reproject(
                source=src_data,
                destination=dst_data,
                src_transform=src.transform,
                src_crs=src_crs,
                src_nodata=nodata,
                dst_transform=template_profile["transform"],
                dst_crs=template_profile["crs"],
                dst_nodata=nodata,
                resampling=Resampling.nearest
            )

            # 掩膜过滤（可选）
            if mask_array is not None:
                if mask_array.shape != dst_data.shape:
                    print(f"  [WARN] 掩膜尺寸不匹配: mask {mask_array.shape} vs data {dst_data.shape}")
                else:
                    if nodata is None:
                        nodata = -9999
                    out_of_mask_count = np.sum(
                        ~mask_array & (dst_data != nodata) & np.isfinite(dst_data) & (dst_data > 0)
                    )
                    dst_data[~mask_array] = nodata
                    if out_of_mask_count > 0:
                        print(f"  [MASK] 已移除 {out_of_mask_count} 个out-of-mask像元")

            profile = src.profile.copy()
            profile.update({
                'crs': template_profile["crs"],
                'transform': template_profile["transform"],
                'width': template_profile["width"],
                'height': template_profile["height"],
                'nodata': nodata
            })

        with rasterio.open(output_file, 'w', **profile) as dst:
            dst.write(dst_data, 1)

        print(f"  [OK] 重投影成功: {file_path.name} ({src_crs} -> {TARGET_CRS})")
        return True

    except Exception as e:
        print(f"  [FAIL] 重投影失败: {file_path.name}")
        print(f"        错误: {e}")
        return False


def reproject_phenology_data(template_profile, mask_array):
    """批量重投影物候数据到EPSG:4326"""
    print("\n" + "="*70)
    print("物候数据CRS统一工具: Clarke 1866 -> EPSG:4326 + 掩膜过滤")
    print("="*70)

    if not INPUT_PHENO_DIR.exists():
        print(f"\n[ERROR] 物候目录不存在: {INPUT_PHENO_DIR}")
        return

    output_dir = OUTPUT_PHENO_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n源目录: {INPUT_PHENO_DIR}")
    print(f"输出目录: {output_dir}")
    print(f"年份范围: {YEAR_START}-{YEAR_END}")
    print(f"目标CRS: {TARGET_CRS}")

    files_to_process = []
    pheno_vars = ['SOS', 'POS', 'EOS']
    for year in range(YEAR_START, YEAR_END + 1):
        for var in pheno_vars:
            if var in PHENO_FILE_FORMAT:
                filename = PHENO_FILE_FORMAT[var].format(year=year)
                file_path = INPUT_PHENO_DIR / filename
                files_to_process.append(file_path)

    print(f"\n待处理文件数: {len(files_to_process)}")
    print("\n[INFO] 此操作将保留原始文件，重投影结果保存到:")
    print(f"       {output_dir}")

    if '--yes' not in sys.argv and not REPROJECT_AUTO_YES:
        response = input("\n是否继续？(yes/no): ").strip().lower()
        if response not in ['yes', 'y']:
            print("\n[CANCELLED] 操作已取消")
            return
    else:
        auto_reason = "--yes参数" if '--yes' in sys.argv else "REPROJECT_AUTO_YES=True"
        print(f"\n[AUTO] 自动确认 ({auto_reason})")

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

        if reproject_file_to_new_dir(file_path, output_dir, template_profile, mask_array):
            success_count += 1
        else:
            fail_count += 1

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

# ==================== 6. 数据验证（整合 _verify_data.py） ====================
def check_directory(dir_path, name):
    """检查目录是否存在"""
    if dir_path.exists():
        print(f"  ✓ {name}: {dir_path}")
        return True
    print(f"  ✗ {name}: {dir_path}")
    print("    目录不存在!")
    return False


def check_file(file_path, name):
    """检查文件是否存在"""
    if file_path.exists():
        print(f"  ✓ {name}: {file_path.name}")
        return True
    print(f"  ✗ {name}: {file_path}")
    print("    文件不存在!")
    return False


def check_phenology_data():
    """检查物候数据完整性"""
    print("\n[2] 检查物候数据 (SOS/POS/EOS):")

    years = list(range(YEAR_START, YEAR_END + 1))
    available = 0
    missing_years = []

    for year in years:
        sos_file = PHENO_DIR / PHENO_FILE_FORMAT['SOS'].format(year=year)
        pos_file = PHENO_DIR / PHENO_FILE_FORMAT['POS'].format(year=year)
        eos_file = PHENO_DIR / PHENO_FILE_FORMAT['EOS'].format(year=year)

        if all([sos_file.exists(), pos_file.exists(), eos_file.exists()]):
            available += 1
        else:
            missing_years.append(year)

    print(f"  可用年份: {available}/{len(years)}")

    if missing_years:
        print(f"  缺失年份: {missing_years[:5]}" + (" ..." if len(missing_years) > 5 else ""))
        return False

    print("  ✓ 所有年份完整")
    return True


def check_TR_data():
    """检查TR数据（抽样检查）"""
    print("\n[3] 检查TR数据 (ERA5-Land格式，抽样检查):")

    test_years = [1982, 1990, 2000, 2010, 2018]
    available = 0

    for year in test_years:
        test_date = datetime(year, 1, 15)
        tr_file = get_TR_file_path(test_date)

        if tr_file and tr_file.exists():
            print(f"  ✓ {year}: {tr_file.name}")
            available += 1
        else:
            expected_name = TR_FILE_FORMAT.format(date=f"{year}0115")
            print(f"  ✗ {year}: 未找到 {expected_name}")

    if available == len(test_years):
        print(f"  ✓ 抽样检查通过 ({available}/{len(test_years)})")
        return True

    print(f"  ⚠ 部分年份缺失 ({available}/{len(test_years)})")
    return False


def check_GPP_data():
    """检查GPP数据（抽样检查）"""
    print("\n[4] 检查GPP数据 (日尺度，抽样检查):")

    test_dates = [
        datetime(1982, 1, 1),
        datetime(2000, 6, 15),
        datetime(2018, 12, 31)
    ]

    available = 0
    for test_date in test_dates:
        gpp_file = GPP_DAILY_DIR / f"GPP_{test_date.strftime('%Y%m%d')}.tif"

        if gpp_file.exists():
            print(f"  ✓ {test_date.date()}: {gpp_file.name}")
            available += 1
        else:
            print(f"  ✗ {test_date.date()}: {gpp_file.name}")

    if available == len(test_dates):
        print(f"  ✓ 抽样检查通过 ({available}/{len(test_dates)})")
        return True

    print(f"  ⚠ 部分日期缺失 ({available}/{len(test_dates)})")
    return False


def check_SM_data():
    """检查土壤水分数据（抽样检查）"""
    print("\n[5] 检查土壤水分数据 (SMrz深层，日尺度):")

    test_dates = [
        datetime(1982, 1, 15),
        datetime(2000, 6, 15),
        datetime(2018, 12, 15)
    ]

    available = 0
    for test_date in test_dates:
        sm_file = SM_DAILY_DIR / f"SMrz_{test_date.strftime('%Y%m%d')}.tif"

        if sm_file.exists():
            print(f"  ✓ {test_date.date()}: {sm_file.name}")
            available += 1
        else:
            print(f"  ✗ {test_date.date()}: {sm_file.name}")

    if available == len(test_dates):
        print(f"  ✓ 抽样检查通过 ({available}/{len(test_dates)})")
        return True

    print(f"  ⚠ 部分日期缺失 ({available}/{len(test_dates)})")
    return False


def run_data_verification():
    """运行完整的数据验证流程"""
    print("\n" + "=" * 70)
    print("Wang2025 数据验证")
    print("=" * 70)

    all_checks = []

    print("\n[1] 检查关键目录:")
    all_checks.append(check_directory(ROOT, "根目录"))
    all_checks.append(check_directory(TR_DAILY_DIR, "TR数据目录"))
    all_checks.append(check_directory(PHENO_DIR, "物候数据目录"))
    all_checks.append(check_directory(GPP_DAILY_DIR, "GPP数据目录"))
    all_checks.append(check_directory(SM_DAILY_DIR, "土壤水分目录"))

    print(f"\n  土地覆盖数据 (USE_FOREST_MASK={USE_FOREST_MASK}):")
    if USE_FOREST_MASK:
        all_checks.append(check_file(LANDCOVER_FILE, "MODIS IGBP"))
    else:
        print("    ⚠ 跳过检查（当前为全局分析模式，仅在仅森林模式时需要）")

    all_checks.append(check_phenology_data())
    all_checks.append(check_TR_data())
    all_checks.append(check_GPP_data())
    all_checks.append(check_SM_data())

    print("\n" + "=" * 70)
    print("验证结果总结")
    print("=" * 70)

    passed = sum(all_checks)
    total = len(all_checks)
    print(f"\n通过检查: {passed}/{total}")

    if passed == total:
        print("\n✓ 所有数据验证通过！")
        return True

    print(f"\n✗ 有 {total - passed} 项检查未通过")
    print("\n请检查:")
    print("  1. _config.py 中的路径配置是否正确")
    print("  2. 数据文件是否已下载到正确位置")
    print("  3. 文件命名格式是否与配置匹配")
    if RUN_VERIFY_STRICT:
        sys.exit(1)
    return False

# ==================== 保存掩膜 ====================
def save_mask(file_path, mask, profile):
    """保存掩膜为GeoTIFF"""
    profile.update(
        dtype=rasterio.uint8,
        count=1,
        compress='lzw',
        nodata=0
    )

    with rasterio.open(file_path, 'w', **profile) as dst:
        dst.write(mask.astype(np.uint8), 1)

    print(f"  ✓ 已保存: {file_path.name}")

# ==================== 主程序 ====================
def main():
    print("\n" + "="*70)
    print("Module 00: 数据准备")
    print("="*70)

    # 显示掩膜策略
    print(f"\n掩膜策略: {'仅森林类型' if USE_FOREST_MASK else '所有陆地（≥30°N）'}")
    if not USE_FOREST_MASK:
        print("  → 将仅使用纬度掩膜，不筛选森林类型")
        print("  → 后续可按植被类型分层分析")

    # 获取模板文件（统一网格：TR日尺度）
    print("\n获取模板文件...")
    template_file = get_template_file()

    if template_file is None:
        print("  ✗ 错误：找不到TR日尺度模板文件")
        print(f"  请检查TR数据路径: {TR_DAILY_DIR}")
        return

    print(f"  ✓ 使用模板: {template_file.name}")
    template_profile = load_template_profile(template_file)

    # CRS检查：必须是地理坐标系
    if template_profile.get("crs") is None or not template_profile["crs"].is_geographic:
        raise ValueError(f"模板CRS不是地理坐标系: {template_profile.get('crs')}")

    # 快速一致性检查（fail-fast）
    fast_grid_consistency_check(template_profile)

    # 1. 创建纬度掩膜
    lat_mask, profile = create_lat_mask(template_file, lat_min=LAT_MIN)
    save_mask(OUTPUT_DIR / "lat_mask.tif", lat_mask, profile)

    # 2. 创建森林掩膜（可选）
    if USE_FOREST_MASK:
        print("\n  使用森林掩膜模式...")
        forest_mask, forest_profile = create_forest_mask(LANDCOVER_FILE, template_profile, FOREST_CLASSES)
        profile = forest_profile

        save_mask(OUTPUT_DIR / "forest_mask.tif", forest_mask, profile)

        # 3. 合并掩膜
        combined_mask = combine_masks(lat_mask, forest_mask)
    else:
        print("\n  跳过森林掩膜（USE_FOREST_MASK=False）")
        print("  组合掩膜将等于纬度掩膜")

        # 创建一个占位的森林掩膜（全1）用于记录
        forest_mask = np.ones_like(lat_mask)
        save_mask(OUTPUT_DIR / "forest_mask.tif", forest_mask, profile)

        # 组合掩膜直接等于纬度掩膜
        combined_mask = lat_mask.copy()

    combined_mask = _apply_data_intersection_mask(combined_mask.astype(bool), template_profile).astype(np.uint8)
    save_mask(MASK_FILE, combined_mask, profile)

    # 4. 物候重投影（可选）
    if RUN_REPROJECT_PHENO:
        reproject_phenology_data(template_profile, combined_mask.astype(bool))

    # 5. 数据完整性检查（快速）
    check_data_availability(years=[1982, 2000, 2010, 2018])

    # 6. 数据验证（全面）
    if RUN_VERIFY_DATA:
        run_data_verification()

    print("\n" + "="*70)
    print("✓ 数据准备完成！")
    print(f"输出目录: {OUTPUT_DIR}")
    print("\n掩膜文件:")
    print(f"  - lat_mask.tif       (纬度 ≥{LAT_MIN}°N)")
    if USE_FOREST_MASK:
        print(f"  - forest_mask.tif    (森林类型: {FOREST_CLASSES})")
        print(f"  - combined_mask.tif  (纬度 AND 森林 AND 数据交集)")
    else:
        print(f"  - forest_mask.tif    (占位文件，全1)")
        print(f"  - combined_mask.tif  (纬度 AND 数据交集)")
    print("="*70)

# ==================== 运行 ====================
if __name__ == "__main__":
    main()
