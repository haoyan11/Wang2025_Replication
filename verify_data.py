#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
数据验证脚本：检查所有必需的数据文件是否存在
运行此脚本以验证数据准备是否完整
"""

import sys
from pathlib import Path
from datetime import datetime
from config import (
    ROOT, TR_DAILY_DIR, PHENO_DIR, GPP_DAILY_DIR, SM_DAILY_DIR,
    LANDCOVER_FILE, YEAR_START, YEAR_END, USE_FOREST_MASK, get_TR_file_path
)

def check_directory(dir_path, name):
    """检查目录是否存在"""
    if dir_path.exists():
        print(f"  ✓ {name}: {dir_path}")
        return True
    else:
        print(f"  ✗ {name}: {dir_path}")
        print(f"    目录不存在!")
        return False

def check_file(file_path, name):
    """检查文件是否存在"""
    if file_path.exists():
        print(f"  ✓ {name}: {file_path.name}")
        return True
    else:
        print(f"  ✗ {name}: {file_path}")
        print(f"    文件不存在!")
        return False

def check_phenology_data():
    """检查物候数据完整性"""
    print("\n[2] 检查物候数据 (SOS/POS/EOS):")

    years = list(range(YEAR_START, YEAR_END + 1))
    available = 0
    missing_years = []

    for year in years:
        sos_file = PHENO_DIR / f"sos_gpp_{year}.tif"
        pos_file = PHENO_DIR / f"pos_doy_gpp_{year}.tif"
        eos_file = PHENO_DIR / f"eos_gpp_{year}.tif"

        if all([sos_file.exists(), pos_file.exists(), eos_file.exists()]):
            available += 1
        else:
            missing_years.append(year)

    print(f"  可用年份: {available}/{len(years)}")

    if missing_years:
        print(f"  缺失年份: {missing_years[:5]}" + (" ..." if len(missing_years) > 5 else ""))
        return False
    else:
        print(f"  ✓ 所有年份完整")
        return True

def check_TR_data():
    """检查TR数据（抽样检查）"""
    print("\n[3] 检查TR数据 (ERA5-Land格式，抽样检查):")

    # 检查几个代表性年份的1月15日
    test_years = [1982, 1990, 2000, 2010, 2018]
    available = 0

    for year in test_years:
        test_date = datetime(year, 1, 15)
        tr_file = get_TR_file_path(test_date)

        if tr_file and tr_file.exists():
            print(f"  ✓ {year}: {tr_file.name}")
            available += 1
        else:
            expected_name = f"ERA5L_ET_transp_Daily_mm_{year}0115.tif"
            print(f"  ✗ {year}: 未找到 {expected_name}")

    if available == len(test_years):
        print(f"  ✓ 抽样检查通过 ({available}/{len(test_years)})")
        return True
    else:
        print(f"  ⚠ 部分年份缺失 ({available}/{len(test_years)})")
        return False

def check_GPP_data():
    """检查GPP数据（抽样检查）"""
    print("\n[4] 检查GPP数据 (日尺度，抽样检查):")

    # 检查几个代表性日期
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
    else:
        print(f"  ⚠ 部分日期缺失 ({available}/{len(test_dates)})")
        return False

def check_SM_data():
    """检查土壤水分数据（抽样检查）"""
    print("\n[5] 检查土壤水分数据 (SMrz深层，日尺度):")

    # 检查几个代表性日期
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
    else:
        print(f"  ⚠ 部分日期缺失 ({available}/{len(test_dates)})")
        return False

def main():
    print("=" * 70)
    print("Wang2025 数据验证脚本")
    print("=" * 70)

    all_checks = []

    # [1] 检查关键目录
    print("\n[1] 检查关键目录:")
    all_checks.append(check_directory(ROOT, "根目录"))
    all_checks.append(check_directory(TR_DAILY_DIR, "TR数据目录"))
    all_checks.append(check_directory(PHENO_DIR, "物候数据目录"))
    all_checks.append(check_directory(GPP_DAILY_DIR, "GPP数据目录"))
    all_checks.append(check_directory(SM_DAILY_DIR, "土壤水分目录"))

    # 检查土地覆盖文件（仅在USE_FOREST_MASK=True时需要）
    print(f"\n  土地覆盖数据 (USE_FOREST_MASK={USE_FOREST_MASK}):")
    if USE_FOREST_MASK:
        all_checks.append(check_file(LANDCOVER_FILE, "MODIS IGBP"))
    else:
        print(f"    ⚠ 跳过检查（当前为全局分析模式，仅在仅森林模式时需要）")
        # 不添加到all_checks，避免影响验证结果

    # [2-5] 检查各类数据
    all_checks.append(check_phenology_data())
    all_checks.append(check_TR_data())
    all_checks.append(check_GPP_data())
    all_checks.append(check_SM_data())

    # 总结
    print("\n" + "=" * 70)
    print("验证结果总结")
    print("=" * 70)

    passed = sum(all_checks)
    total = len(all_checks)

    print(f"\n通过检查: {passed}/{total}")

    if passed == total:
        print("\n✓ 所有数据验证通过！")
        print("\n您可以开始运行分析流程:")
        print("  1. python 00_data_preparation.py    # 准备掩膜")
        print("  2. python 00_master_pipeline.py     # 运行完整分析")
        return 0
    else:
        print(f"\n✗ 有 {total - passed} 项检查未通过")
        print("\n请检查:")
        print("  1. config.py 中的路径配置是否正确")
        print("  2. 数据文件是否已下载到正确位置")
        print("  3. 文件命名格式是否与配置匹配")
        print("\n详细配置说明请查看: 数据配置说明.md")
        return 1

if __name__ == "__main__":
    sys.exit(main())
