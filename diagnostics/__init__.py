"""
诊断工具模块

用于验证Wang2025复现项目中各步骤输出数据的正确性

可用工具:
- check_decomposition: 分解数据诊断（TR_fixed_window值分布、enhancement_factor等）
"""

__version__ = "1.0.0"

from pathlib import Path

# 导出模块路径
DIAGNOSTICS_DIR = Path(__file__).parent
PROJECT_ROOT = DIAGNOSTICS_DIR.parent
