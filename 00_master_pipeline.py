#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Master Pipeline: Wang (2025) 完整分析流程
按顺序执行所有分析模块
"""

import sys
import subprocess
from pathlib import Path
from datetime import datetime
import logging

# ==================== 配置日志 ====================
LOG_DIR = Path(r"I:\F\Data4\Wang2025_Analysis\Logs")
LOG_DIR.mkdir(parents=True, exist_ok=True)

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
log_file = LOG_DIR / f"pipeline_{timestamp}.log"

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.FileHandler(log_file, encoding='utf-8'),
        logging.StreamHandler(sys.stdout)
    ]
)

logger = logging.getLogger(__name__)

# ==================== 模块列表 ====================
MODULES = [
    {
        'name': 'Module 00: 数据准备',
        'script': '00_data_preparation.py',
        'required': True,
        'description': '创建掩膜、气候平均态等'
    },
    {
        'name': 'Module 01: 物候提取',
        'script': '01_phenology_extraction.py',
        'required': True,
        'description': '从SIF提取SOS/POS/EOS'
    },
    {
        'name': 'Module 02: TRc计算',
        'script': '02_TRc_calculation.py',
        'required': True,
        'description': '计算SOS-POS窗口内累积蒸腾'
    },
    {
        'name': 'Module 03a: Wang 2025分解（GPP物候）',
        'script': '03a_decomposition_wang2025.py',
        'required': True,
        'description': 'TRpheno + TRproduct 分解（Wang 2025原始方法）'
    },
    {
        'name': 'Module 03b: Timing/Shape分解（可选）',
        'script': '03b_decomposition_timing_shape.py',
        'required': False,
        'description': 'TRtiming + TRshape 分解（新方法）'
    },
    {
        'name': 'Module 04a: 统计分析（Wang 2025方法）',
        'script': '04a_statistical_wang2025.py',
        'required': True,
        'description': 'ΔSOS回归分析（对应03a分解）'
    },
    {
        'name': 'Module 04b: 统计分析（Timing/Shape）',
        'script': '04b_statistical_timing_shape.py',
        'required': False,
        'description': 'Timing/Shape回归分析（对应03b分解）'
    },
    {
        'name': 'Module 05: SEM分析',
        'script': '05_SEM_analysis.R',
        'required': False,
        'description': '结构方程模型（需要R环境）',
        'interpreter': 'Rscript'
    },
    {
        'name': 'Module 06: 绘图',
        'script': '06_plotting.py',
        'required': True,
        'description': '生成所有图表'
    }
]

# ==================== 执行函数 ====================
def run_module(module_info):
    """执行单个模块"""
    logger.info("="*70)
    logger.info(f"开始执行: {module_info['name']}")
    logger.info(f"描述: {module_info['description']}")
    logger.info("="*70)

    script_path = Path(__file__).parent / module_info['script']

    if not script_path.exists():
        logger.error(f"脚本文件不存在: {script_path}")
        return False

    # 确定解释器
    interpreter = module_info.get('interpreter', 'python')

    try:
        # 执行脚本
        if interpreter == 'Rscript':
            cmd = ['Rscript', str(script_path)]
        else:
            cmd = [sys.executable, str(script_path)]

        logger.info(f"执行命令: {' '.join(cmd)}")

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='replace'
        )

        # 输出结果
        if result.stdout:
            logger.info("标准输出:")
            for line in result.stdout.splitlines():
                logger.info(f"  {line}")

        if result.stderr:
            logger.warning("标准错误:")
            for line in result.stderr.splitlines():
                logger.warning(f"  {line}")

        # 检查返回码
        if result.returncode != 0:
            logger.error(f"模块执行失败，返回码: {result.returncode}")
            return False

        logger.info(f"✓ {module_info['name']} 执行成功")
        return True

    except Exception as e:
        logger.error(f"执行出错: {str(e)}")
        return False

def check_dependencies():
    """检查依赖环境"""
    logger.info("\n检查依赖环境...")

    # 检查Python包
    python_packages = [
        'numpy', 'rasterio', 'scipy', 'matplotlib',
        'cartopy', 'pandas', 'tqdm'
    ]

    missing_packages = []
    for pkg in python_packages:
        try:
            __import__(pkg)
            logger.info(f"  ✓ {pkg}")
        except ImportError:
            logger.warning(f"  ✗ {pkg} (缺失)")
            missing_packages.append(pkg)

    if missing_packages:
        logger.warning(f"\n缺少Python包: {', '.join(missing_packages)}")
        logger.warning("请运行: pip install " + ' '.join(missing_packages))

    # 检查R环境
    try:
        result = subprocess.run(['Rscript', '--version'],
                              capture_output=True, text=True)
        if result.returncode == 0:
            logger.info("  ✓ R环境")
        else:
            logger.warning("  ✗ R环境 (未安装或不在PATH中)")
    except FileNotFoundError:
        logger.warning("  ✗ R环境 (未找到Rscript命令)")

    logger.info("")

# ==================== 主流程 ====================
def main(skip_modules=None, only_modules=None):
    """
    主执行流程

    Parameters:
    -----------
    skip_modules : list
        跳过的模块编号列表（如 [4, 7]）
    only_modules : list
        仅执行的模块编号列表（如 [3, 5]）
    """
    logger.info("\n" + "="*70)
    logger.info("Wang (2025) 完整分析流程")
    logger.info("="*70)
    logger.info(f"开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"日志文件: {log_file}")
    logger.info("")

    # 检查依赖
    check_dependencies()

    # 执行模块
    success_count = 0
    failed_modules = []

    for idx, module_info in enumerate(MODULES):
        module_num = idx  # 模块编号从0开始

        # 跳过逻辑
        if skip_modules and module_num in skip_modules:
            logger.info(f"跳过 {module_info['name']} (用户指定)")
            continue

        if only_modules and module_num not in only_modules:
            logger.info(f"跳过 {module_info['name']} (非指定模块)")
            continue

        # 执行模块
        success = run_module(module_info)

        if success:
            success_count += 1
        else:
            failed_modules.append(module_info['name'])

            # 如果是必需模块且失败，询问是否继续
            if module_info.get('required', True):
                logger.error(f"\n必需模块 {module_info['name']} 执行失败")
                user_input = input("是否继续执行后续模块？[y/N]: ")
                if user_input.lower() != 'y':
                    logger.error("用户中止流程")
                    break

        logger.info("")

    # 总结
    logger.info("="*70)
    logger.info("流程执行完成")
    logger.info("="*70)
    logger.info(f"成功模块数: {success_count}/{len(MODULES)}")

    if failed_modules:
        logger.warning("失败模块:")
        for module_name in failed_modules:
            logger.warning(f"  - {module_name}")
    else:
        logger.info("✓ 所有模块执行成功！")

    logger.info(f"结束时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"日志文件: {log_file}")
    logger.info("="*70)

# ==================== 命令行接口 ====================
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Wang (2025) 分析流程')
    parser.add_argument('--skip', type=int, nargs='+',
                       help='跳过指定模块（如 --skip 4 7）')
    parser.add_argument('--only', type=int, nargs='+',
                       help='仅执行指定模块（如 --only 3 5）')

    args = parser.parse_args()

    main(skip_modules=args.skip, only_modules=args.only)
