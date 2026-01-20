#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module 03C: Fixed-Window decomposition (Method 2)

方法2：实际窗口变化 + 固定窗口强度分解

Decompose:
  TRc_y = TRc_av + TR_window_change + TR_fixed_window

Where:
  TRc_y   = sum_{t=SOS_y..POS_y} B_y(t)  # 当年实际累积
  TRc_av  = sum_{t=SOSav..POSav} A(t)    # 气候态基线（固定窗口）

  TR_window_change = 窗口变化带来的实际额外累积
    = TRc_y_estimated - TRc_y_fixed_window
    （窗口外部分按当年强度估算）

  TR_fixed_window = 固定窗口内的强度差异
    = sum_{t=SOSav..POSav} [B_y(t) - A(t)]
    （在固定窗口内，当年与气候态的累积差异）

Optional split:
  TR_window_change = TR_sos_change + TR_pos_change
  TR_sos_change = SOS变化带来的实际累积（如果SOS_y < SOSav，为正）
  TR_pos_change = POS变化带来的实际累积（如果POS_y > POSav，为正）

核心思想：
  - 固定窗口[SOSav, POSav]内的强度变化可以剥离窗口选择效应
  - TR_fixed_window/(POSav-SOSav+1) = 真实的平均速率差异
  - 这是评估"蒸腾速率是否真的变高变低"的最佳指标

实现方法（使用当年每日TR数据）：
  直接使用当年日尺度TR（B_y(t)）进行分解，无比例缩放假设：
  1. TRc_y = sum_{SOS_y..POS_y} B_y(t)
  2. TR_fixed_window = sum_{SOSav..POSav} [B_y(t) - A(t)]
  3. TR_window_change = TRc_y - sum_{SOSav..POSav} B_y(t)

Notes:
  - A(t) is the climatology TR_daily_av (365 bands).
  - B_y(t) is the year-specific daily TR (observed daily TR, no scaling).
  - 固定窗口为[SOSav, POSav]，所有年份一致
"""

import numpy as np
import rasterio
from rasterio.windows import Window
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
from datetime import datetime, timedelta
import os

# 导入配置
from _config import (
    ROOT, PHENO_DIR, TRC_ANNUAL_DIR, CLIMATOLOGY_DIR, DECOMPOSITION_FIXED_DIR,
    YEAR_START, YEAR_END, NODATA_OUT, TR_DAILY_DIR, TR_FILE_FORMAT, BLOCK_SIZE,
    MAX_WORKERS, TEMPLATE_RASTER, MASK_FILE
)

# GDAL缓存（MB），用于加速连续读取
GDAL_CACHE_MAX_MB = 512
PARALLEL_BY_YEAR = True
YEAR_WORKERS = MAX_WORKERS
MIN_VALID_FRAC = 0.60  # 有效天数阈值（窗口内）

# 确保输出目录存在
DECOMPOSITION_FIXED_DIR.mkdir(parents=True, exist_ok=True)

# 向后兼容：保留旧变量名
TRC_DIR = TRC_ANNUAL_DIR
CLIM_DIR = CLIMATOLOGY_DIR
OUTPUT_DIR = DECOMPOSITION_FIXED_DIR

# 运行模式控制
RUN_MODE = os.getenv("WANG_RUN_MODE", "skip").strip().lower()
OVERWRITE = RUN_MODE == "overwrite"

def should_write(path):
    path = Path(path)
    return OVERWRITE or not path.exists()


def _is_valid_value(value, nodata):
    """
    判断值是否有效（非NODATA）
    兼容各种NODATA类型：-9999, NaN, 3.4e38等
    当nodata元信息缺失时，使用绝对值阈值过滤极大填充值
    """
    # 基础检查：有限性
    valid = np.isfinite(value)

    # 如果有明确的nodata值，检查是否等于nodata
    if nodata is not None and not np.isnan(nodata):
        valid = valid & (value != nodata)

    # 额外保护：过滤常见的极大填充值（1e20, 3.4e38等）
    valid = valid & (np.abs(value) < 1e10)
    # nodata缺失时的兜底：常见无效值为-9999
    if nodata is None or np.isnan(nodata):
        valid = valid & (value > -9000)

    return valid


_GLOBAL_CTX = {}


def _init_worker(trc_av, sos_av, pos_av, nodata_sos, nodata_pos, mask, tr_valid, profile):
    global _GLOBAL_CTX
    _GLOBAL_CTX = {
        "trc_av": trc_av,
        "sos_av": sos_av,
        "pos_av": pos_av,
        "nodata_sos": nodata_sos,
        "nodata_pos": nodata_pos,
        "mask": mask,
        "tr_valid": tr_valid,
        "profile": profile,
    }


def _process_year(year):
    if not OVERWRITE:
        outputs = [
            OUTPUT_DIR / f"TR_window_change_{year}.tif",
            OUTPUT_DIR / f"TR_fixed_window_{year}.tif",
            OUTPUT_DIR / f"TR_sos_change_{year}.tif",
            OUTPUT_DIR / f"TR_pos_change_{year}.tif",
        ]
        if all(p.exists() for p in outputs):
            return year
    ctx = _GLOBAL_CTX
    tr_window_change, tr_fixed_window, tr_sos_change, tr_pos_change = decompose_year(
        year,
        ctx["trc_av"],
        ctx["sos_av"],
        ctx["pos_av"],
        ctx["nodata_sos"],
        ctx["nodata_pos"],
        ctx["mask"],
        ctx["tr_valid"],
        ctx["profile"]
    )
    write_geotiff(OUTPUT_DIR / f"TR_window_change_{year}.tif", tr_window_change, ctx["profile"])
    write_geotiff(OUTPUT_DIR / f"TR_fixed_window_{year}.tif", tr_fixed_window, ctx["profile"])
    write_geotiff(OUTPUT_DIR / f"TR_sos_change_{year}.tif", tr_sos_change, ctx["profile"])
    write_geotiff(OUTPUT_DIR / f"TR_pos_change_{year}.tif", tr_pos_change, ctx["profile"])
    return year


def read_geotiff(file_path):
    with rasterio.open(file_path) as src:
        data = src.read(1)
        profile = src.profile.copy()
        nodata = src.nodata
    return data, profile, nodata


def load_template_profile():
    if not TEMPLATE_RASTER.exists():
        raise FileNotFoundError(f"模板栅格不存在: {TEMPLATE_RASTER}")
    with rasterio.open(TEMPLATE_RASTER) as src:
        return src.profile.copy()


def write_geotiff(file_path, data, profile):
    """写入单波段GeoTIFF（不修改原profile）"""
    if not should_write(file_path):
        print(f"  ⚠ 跳过已存在: {Path(file_path).name}")
        return
    out_profile = profile.copy()
    out_profile.update(dtype=rasterio.float32, count=1, compress="lzw", nodata=NODATA_OUT)
    with rasterio.open(file_path, "w", **out_profile) as dst:
        dst.write(data.astype(np.float32), 1)


def _sum_cum_range(cum_arr, start_doy, end_doy):
    """
    Sum over inclusive [start, end], DOY is 1-based.
    cum_arr: shape (365, H, W) - cumulative sum along DOY axis
    start_doy, end_doy: shape (H, W) - DOY values (1-365)
    """
    flat = cum_arr.reshape(cum_arr.shape[0], -1)
    start = start_doy.ravel()
    end = end_doy.ravel()
    idx = np.arange(start.size)
    end_vals = flat[end - 1, idx]
    start_vals = np.where(start > 1, flat[start - 2, idx], 0.0)
    return (end_vals - start_vals).reshape(start_doy.shape)


def is_leap_year(year):
    """判断是否为闰年"""
    return (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)


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


def accumulate_daily_tr_sums(year, ref_profile, valid_base,
                             sos_y_i, pos_y_i, sos_av_i, pos_av_i):
    """
    按日读取TR并累积所需窗口总量，避免365天三维数组占用内存。
    返回:
      trc_y_actual, trc_y_fixed_actual, tr_sos_advance, tr_sos_delay,
      tr_pos_extend, tr_pos_shorten, tr_valid_y, trc_y_count, trc_y_fixed_count
    """
    height = ref_profile["height"]
    width = ref_profile["width"]

    trc_y_actual = np.zeros((height, width), dtype=np.float32)
    trc_y_fixed_actual = np.zeros((height, width), dtype=np.float32)
    tr_sos_advance = np.zeros((height, width), dtype=np.float32)
    tr_sos_delay = np.zeros((height, width), dtype=np.float32)
    tr_pos_extend = np.zeros((height, width), dtype=np.float32)
    tr_pos_shorten = np.zeros((height, width), dtype=np.float32)
    tr_valid_y = np.zeros((height, width), dtype=bool)
    trc_y_count = np.zeros((height, width), dtype=np.int16)
    trc_y_fixed_count = np.zeros((height, width), dtype=np.int16)

    block_windows = [
        Window(col_off, row_off,
               min(BLOCK_SIZE, width - col_off),
               min(BLOCK_SIZE, height - row_off))
        for row_off in range(0, height, BLOCK_SIZE)
        for col_off in range(0, width, BLOCK_SIZE)
    ]
    block_has_valid = []
    for win in block_windows:
        r0, r1 = win.row_off, win.row_off + win.height
        c0, c1 = win.col_off, win.col_off + win.width
        block_has_valid.append(np.any(valid_base[r0:r1, c0:c1]))

    days_in_year = 366 if is_leap_year(year) else 365
    checked_ref = False

    with rasterio.Env(GDAL_CACHEMAX=GDAL_CACHE_MAX_MB):
        for day in range(1, days_in_year + 1):
            date_obj = datetime(year, 1, 1) + timedelta(days=day - 1)
            doy_noleap = noleap_doy_from_date(date_obj)
            if doy_noleap is None:
                continue

            file_path = TR_DAILY_DIR / TR_FILE_FORMAT.format(date=date_obj.strftime("%Y%m%d"))
            if not file_path.exists():
                continue

            with rasterio.open(file_path) as src:
                profile = src.profile.copy()
                nodata = src.nodata

                if not checked_ref:
                    check_spatial_consistency(ref_profile, file_path, profile,
                                              f"TR_daily_{date_obj.strftime('%Y%m%d')}")
                    checked_ref = True

                for win, has_valid in zip(block_windows, block_has_valid):
                    if not has_valid:
                        continue

                    r0, r1 = win.row_off, win.row_off + win.height
                    c0, c1 = win.col_off, win.col_off + win.width
                    base_block = valid_base[r0:r1, c0:c1]
                    if not np.any(base_block):
                        continue

                    data = src.read(1, window=win).astype(np.float32)
                    if nodata is None or np.isnan(nodata):
                        data[~np.isfinite(data)] = np.nan
                    else:
                        data[data == nodata] = np.nan
                    data[data < 0] = np.nan

                    valid_day = np.isfinite(data)
                    if not np.any(valid_day):
                        continue

                    valid_day = valid_day & base_block
                    if not np.any(valid_day):
                        continue

                    tr_valid_y[r0:r1, c0:c1] |= valid_day
                    doy = doy_noleap

                    sos_y_block = sos_y_i[r0:r1, c0:c1]
                    pos_y_block = pos_y_i[r0:r1, c0:c1]
                    sos_av_block = sos_av_i[r0:r1, c0:c1]
                    pos_av_block = pos_av_i[r0:r1, c0:c1]

                    m_trc = valid_day & (sos_y_block <= doy) & (pos_y_block >= doy)
                    if np.any(m_trc):
                        trc_block = trc_y_actual[r0:r1, c0:c1]
                        trc_block[m_trc] += data[m_trc]
                        trc_count_block = trc_y_count[r0:r1, c0:c1]
                        trc_count_block[m_trc] += 1

                    m_fixed = valid_day & (sos_av_block <= doy) & (pos_av_block >= doy)
                    if np.any(m_fixed):
                        fixed_block = trc_y_fixed_actual[r0:r1, c0:c1]
                        fixed_block[m_fixed] += data[m_fixed]
                        fixed_count_block = trc_y_fixed_count[r0:r1, c0:c1]
                        fixed_count_block[m_fixed] += 1

                    m_sos_adv = valid_day & (sos_y_block < sos_av_block) & (doy >= sos_y_block) & (doy < sos_av_block)
                    if np.any(m_sos_adv):
                        sos_adv_block = tr_sos_advance[r0:r1, c0:c1]
                        sos_adv_block[m_sos_adv] += data[m_sos_adv]

                    m_sos_delay = valid_day & (sos_y_block > sos_av_block) & (doy >= sos_av_block) & (doy < sos_y_block)
                    if np.any(m_sos_delay):
                        sos_delay_block = tr_sos_delay[r0:r1, c0:c1]
                        sos_delay_block[m_sos_delay] += data[m_sos_delay]

                    m_pos_extend = valid_day & (pos_y_block > pos_av_block) & (doy > pos_av_block) & (doy <= pos_y_block)
                    if np.any(m_pos_extend):
                        pos_ext_block = tr_pos_extend[r0:r1, c0:c1]
                        pos_ext_block[m_pos_extend] += data[m_pos_extend]

                    m_pos_shorten = valid_day & (pos_y_block < pos_av_block) & (doy > pos_y_block) & (doy <= pos_av_block)
                    if np.any(m_pos_shorten):
                        pos_short_block = tr_pos_shorten[r0:r1, c0:c1]
                        pos_short_block[m_pos_shorten] += data[m_pos_shorten]

    return (trc_y_actual, trc_y_fixed_actual, tr_sos_advance, tr_sos_delay,
            tr_pos_extend, tr_pos_shorten, tr_valid_y, trc_y_count, trc_y_fixed_count)


def check_spatial_consistency(ref_profile, file_path, data_profile, var_name="data"):
    """检查栅格空间一致性（shape/crs/transform）"""
    if ref_profile['width'] != data_profile['width'] or \
       ref_profile['height'] != data_profile['height']:
        raise ValueError(
            f"Shape mismatch for {var_name}: {file_path}\n"
            f"  Expected: {ref_profile['height']}x{ref_profile['width']}\n"
            f"  Got: {data_profile['height']}x{data_profile['width']}"
        )
    if ref_profile['crs'] != data_profile['crs']:
        raise ValueError(
            f"CRS mismatch for {var_name}: {file_path}\n"
            f"  Expected: {ref_profile['crs']}\n"
            f"  Got: {data_profile['crs']}"
        )

    # Transform比较使用容差（避免浮点精度误报）
    ref_transform = ref_profile['transform']
    data_transform = data_profile['transform']
    if ref_transform is not None and data_transform is not None:
        transform_match = all(
            abs(ref_transform[i] - data_transform[i]) < 1e-6
            for i in range(6)
        )
        if not transform_match:
            raise ValueError(
                f"Transform mismatch for {var_name}: {file_path}\n"
                f"  Expected: {ref_transform}\n"
                f"  Got: {data_transform}"
            )


def load_climatology():
    template_profile = load_template_profile()
    tr_clim_file = CLIM_DIR / "TR_daily_climatology.tif"
    sos_av_file = CLIM_DIR / "SOSav.tif"
    pos_av_file = CLIM_DIR / "POSav.tif"

    if not tr_clim_file.exists():
        raise FileNotFoundError(f"Missing climatology: {tr_clim_file}")
    if not sos_av_file.exists():
        raise FileNotFoundError(f"Missing SOSav: {sos_av_file}")
    if not pos_av_file.exists():
        raise FileNotFoundError(f"Missing POSav: {pos_av_file}")

    with rasterio.open(tr_clim_file) as src:
        tr_daily_av = src.read().astype(np.float32)
        profile = template_profile.copy()
        check_spatial_consistency(profile, tr_clim_file, src.profile, "TR_daily_climatology")
        nodata_tr = src.nodata

    if nodata_tr is None or np.isnan(nodata_tr):
        tr_daily_av[~np.isfinite(tr_daily_av)] = np.nan
    else:
        tr_daily_av[tr_daily_av == nodata_tr] = np.nan
    tr_daily_av[tr_daily_av < 0] = np.nan

    sos_av, sos_profile, nodata_sos = read_geotiff(sos_av_file)
    check_spatial_consistency(profile, sos_av_file, sos_profile, "SOSav")

    pos_av, pos_profile, nodata_pos = read_geotiff(pos_av_file)
    check_spatial_consistency(profile, pos_av_file, pos_profile, "POSav")

    return tr_daily_av, sos_av, pos_av, nodata_sos, nodata_pos, profile


def decompose_year(year, trc_av, sos_av, pos_av,
                   nodata_sos, nodata_pos, mask, tr_valid, ref_profile):
    """
    方法2分解：固定窗口强度 + 实际窗口变化

    Returns:
    --------
    tr_window_change : 窗口变化带来的实际累积
    tr_fixed_window : 固定窗口内的强度差异
    tr_sos_change : SOS变化的贡献
    tr_pos_change : POS变化的贡献
    """
    # 读取当年数据
    trc_file = TRC_DIR / f"TRc_{year}.tif"
    sos_file = PHENO_DIR / f"sos_gpp_{year}.tif"
    pos_file = PHENO_DIR / f"pos_doy_gpp_{year}.tif"

    trc_y, trc_profile, nodata_trc = read_geotiff(trc_file)
    sos_y, sos_profile, nodata_sos_y = read_geotiff(sos_file)
    pos_y, pos_profile, nodata_pos_y = read_geotiff(pos_file)

    # 空间一致性检查
    check_spatial_consistency(ref_profile, trc_file, trc_profile, f"TRc_{year}")
    check_spatial_consistency(ref_profile, sos_file, sos_profile, f"sos_gpp_{year}")
    check_spatial_consistency(ref_profile, pos_file, pos_profile, f"pos_doy_gpp_{year}")

    height, width = trc_y.shape

    # 有效性检查（与03a/03b一致）
    # 注意：TRc_{year}.tif 仅用于一致性校验与有效像元筛选，
    #       分解计算使用当年日尺度TR逐日累积得到的 trc_y_actual。
    valid_base = (
        _is_valid_value(trc_y, nodata_trc) &
        (trc_av != NODATA_OUT) &
        _is_valid_value(sos_y, nodata_sos_y) &
        _is_valid_value(pos_y, nodata_pos_y) &
        _is_valid_value(sos_av, nodata_sos) &
        _is_valid_value(pos_av, nodata_pos) &
        (sos_y > 0) & (sos_y <= 365) &
        (pos_y > 0) & (pos_y <= 365) &
        (sos_av > 0) & (sos_av <= 365) &
        (pos_av > 0) & (pos_av <= 365) &
        (pos_y >= sos_y) &
        (pos_av >= sos_av) &
        mask &
        tr_valid
    )

    # 初始化输出
    tr_window_change = np.full((height, width), NODATA_OUT, dtype=np.float32)
    tr_fixed_window = np.full((height, width), NODATA_OUT, dtype=np.float32)
    tr_sos_change = np.full((height, width), NODATA_OUT, dtype=np.float32)
    tr_pos_change = np.full((height, width), NODATA_OUT, dtype=np.float32)

    if not np.any(valid_base):
        return tr_window_change, tr_fixed_window, tr_sos_change, tr_pos_change

    # Clip DOY to 1-365
    sos_y_i = np.clip(np.rint(sos_y).astype(np.int32), 1, 365)
    pos_y_i = np.clip(np.rint(pos_y).astype(np.int32), 1, 365)
    sos_av_i = np.clip(np.rint(sos_av).astype(np.int32), 1, 365)
    pos_av_i = np.clip(np.rint(pos_av).astype(np.int32), 1, 365)

    # 进一步检查窗口有效性
    valid_sos = valid_base & (pos_y_i >= sos_y_i) & (pos_av_i >= sos_av_i)

    if not np.any(valid_sos):
        return tr_window_change, tr_fixed_window, tr_sos_change, tr_pos_change

    # ========== 核心计算（逐日累积，避免365天三维数组） ==========
    (trc_y_actual, trc_y_fixed_actual, tr_sos_advance, tr_sos_delay,
     tr_pos_extend, tr_pos_shorten, tr_valid_y, trc_y_count, trc_y_fixed_count) = accumulate_daily_tr_sums(
        year, ref_profile, valid_base, sos_y_i, pos_y_i, sos_av_i, pos_av_i
    )

    valid_sos = valid_sos & tr_valid_y
    # 有效天数比例过滤（窗口内有效天数需 >= MIN_VALID_FRAC * 窗口长度）
    window_len_y = (pos_y_i - sos_y_i + 1).astype(np.int16)
    window_len_av = (pos_av_i - sos_av_i + 1).astype(np.int16)
    min_required_y = np.ceil(MIN_VALID_FRAC * window_len_y).astype(np.int16)
    min_required_av = np.ceil(MIN_VALID_FRAC * window_len_av).astype(np.int16)
    valid_sos = valid_sos & (trc_y_count >= min_required_y) & (trc_y_fixed_count >= min_required_av)
    if not np.any(valid_sos):
        return tr_window_change, tr_fixed_window, tr_sos_change, tr_pos_change

    # 3. 固定窗口强度差异
    #    TR_fixed_window = sum_{SOSav..POSav} [B_y(t) - A(t)]
    tr_fixed_window[valid_sos] = (trc_y_fixed_actual - trc_av)[valid_sos]

    # 4. 窗口变化部分（残差）
    #    TR_window_change = TRc_y - sum_{SOSav..POSav} B_y(t)
    tr_window_change[valid_sos] = (trc_y_actual - trc_y_fixed_actual)[valid_sos]

    # 5. 进一步分解窗口变化为SOS和POS部分
    #    这部分需要更详细的估算
    #    注意：避免np.where包裹_sum_cum_range，防止边界越界

    # SOS变化效应估算
    # 如果SOS_y < SOSav（提前），计算[SOS_y, SOSav-1]的贡献
    # 如果SOS_y > SOSav（推迟），计算[SOSav, SOS_y-1]的贡献（为负）
    sos_before = sos_y_i < sos_av_i  # SOS提前

    # 根据提前/推迟设置正负号
    tr_sos_change[valid_sos] = np.where(
        sos_before,
        tr_sos_advance,
        -tr_sos_delay
    )[valid_sos]

    # POS变化效应估算（类似SOS）
    pos_after = pos_y_i > pos_av_i   # POS延后

    # 根据延后/提前设置正负号
    tr_pos_change[valid_sos] = np.where(
        pos_after,
        tr_pos_extend,
        -tr_pos_shorten
    )[valid_sos]

    return tr_window_change, tr_fixed_window, tr_sos_change, tr_pos_change


def main():
    print("\n" + "=" * 80)
    print("Fixed-Window Decomposition (Method 2)")
    print("=" * 80)
    print("\n核心思想：")
    print("  - 固定窗口[SOSav, POSav]内的强度变化 → 剥离窗口选择效应")
    print("  - TR_fixed_window/(POSav-SOSav+1) → 真实的平均速率差异")
    print("  - 窗口变化部分 → 窗口延长/缩短带来的实际累积\n")

    years = list(range(YEAR_START, YEAR_END + 1))

    if RUN_MODE == "skip":
        base_outputs = [
            OUTPUT_DIR / "TRc_av.tif",
            OUTPUT_DIR / "Fixed_Window_Length.tif",
        ]
        year_outputs = []
        for year in years:
            year_outputs.extend([
                OUTPUT_DIR / f"TR_window_change_{year}.tif",
                OUTPUT_DIR / f"TR_fixed_window_{year}.tif",
                OUTPUT_DIR / f"TR_sos_change_{year}.tif",
                OUTPUT_DIR / f"TR_pos_change_{year}.tif",
            ])
        if all(p.exists() for p in base_outputs + year_outputs):
            print("  ✓ 输出齐全，跳过 Module 03c")
            return

    # Load mask
    mask_file = MASK_FILE
    if not mask_file.exists():
        raise FileNotFoundError(f"Missing mask: {mask_file}")
    mask, mask_profile, mask_nodata = read_geotiff(mask_file)
    mask = _is_valid_value(mask, mask_nodata) & (mask > 0)

    # Load climatology
    tr_daily_av, sos_av, pos_av, nodata_sos, nodata_pos, profile = load_climatology()

    # 空间一致性检查
    check_spatial_consistency(profile, mask_file, mask_profile, "combined_mask")
    tr_valid = np.any(np.isfinite(tr_daily_av), axis=0)
    tr_cum = np.nancumsum(tr_daily_av, axis=0).astype(np.float32)

    # Compute TRc_av (固定窗口基线)
    sos_av_i = np.clip(np.rint(sos_av).astype(np.int32), 1, 365)
    pos_av_i = np.clip(np.rint(pos_av).astype(np.int32), 1, 365)
    trc_av = np.full(sos_av.shape, NODATA_OUT, dtype=np.float32)
    valid_av = (
        _is_valid_value(sos_av, nodata_sos) &
        _is_valid_value(pos_av, nodata_pos) &
        (sos_av > 0) & (sos_av <= 366) &  # 允许闰年DOY=366，排除异常值
        (pos_av > 0) & (pos_av <= 366) &  # 允许闰年DOY=366，排除异常值
        (pos_av_i >= sos_av_i) &
        mask & tr_valid
    )
    trc_av[valid_av] = _sum_cum_range(tr_cum, sos_av_i, pos_av_i)[valid_av]
    write_geotiff(OUTPUT_DIR / "TRc_av.tif", trc_av, profile)

    # 保存固定窗口长度（用于后续速率计算）
    fixed_window_length = np.full(sos_av.shape, NODATA_OUT, dtype=np.float32)
    fixed_window_length[valid_av] = (pos_av_i - sos_av_i + 1)[valid_av]
    write_geotiff(OUTPUT_DIR / "Fixed_Window_Length.tif", fixed_window_length, profile)

    # Decompose each year
    if PARALLEL_BY_YEAR and YEAR_WORKERS > 1:
        print(f"\n[Parallel] 按年并行分解: {YEAR_WORKERS} 进程")
        with ProcessPoolExecutor(
            max_workers=YEAR_WORKERS,
            initializer=_init_worker,
            initargs=(trc_av, sos_av, pos_av, nodata_sos, nodata_pos, mask, tr_valid, profile)
        ) as executor:
            futures = {executor.submit(_process_year, year): year for year in years}
            for future in tqdm(as_completed(futures), total=len(futures),
                               desc="Fixed-Window Decomposition (parallel)"):
                year = futures[future]
                try:
                    _ = future.result()
                except Exception as exc:
                    raise RuntimeError(f"Year {year} failed: {exc}") from exc
    else:
        for year in tqdm(years, desc="Fixed-Window Decomposition"):
            if not OVERWRITE:
                outputs = [
                    OUTPUT_DIR / f"TR_window_change_{year}.tif",
                    OUTPUT_DIR / f"TR_fixed_window_{year}.tif",
                    OUTPUT_DIR / f"TR_sos_change_{year}.tif",
                    OUTPUT_DIR / f"TR_pos_change_{year}.tif",
                ]
                if all(p.exists() for p in outputs):
                    continue
            tr_window_change, tr_fixed_window, tr_sos_change, tr_pos_change = decompose_year(
                year, trc_av, sos_av, pos_av,
                nodata_sos, nodata_pos, mask, tr_valid, profile
            )
            write_geotiff(OUTPUT_DIR / f"TR_window_change_{year}.tif", tr_window_change, profile)
            write_geotiff(OUTPUT_DIR / f"TR_fixed_window_{year}.tif", tr_fixed_window, profile)
            write_geotiff(OUTPUT_DIR / f"TR_sos_change_{year}.tif", tr_sos_change, profile)
            write_geotiff(OUTPUT_DIR / f"TR_pos_change_{year}.tif", tr_pos_change, profile)

    # Quality check
    print("\n[Quality Check]")
    test_year = 2000
    trc_test, _, nodata_trc_test = read_geotiff(TRC_DIR / f"TRc_{test_year}.tif")
    tr_window_test, _, _ = read_geotiff(OUTPUT_DIR / f"TR_window_change_{test_year}.tif")
    tr_fixed_test, _, _ = read_geotiff(OUTPUT_DIR / f"TR_fixed_window_{test_year}.tif")
    tr_sos_test, _, _ = read_geotiff(OUTPUT_DIR / f"TR_sos_change_{test_year}.tif")
    tr_pos_test, _, _ = read_geotiff(OUTPUT_DIR / f"TR_pos_change_{test_year}.tif")

    sos_y, _, nodata_sos_y = read_geotiff(PHENO_DIR / f"sos_gpp_{test_year}.tif")
    pos_y, _, nodata_pos_y = read_geotiff(PHENO_DIR / f"pos_doy_gpp_{test_year}.tif")

    sos_y_i = np.clip(np.rint(sos_y).astype(np.int32), 1, 365)
    pos_y_i = np.clip(np.rint(pos_y).astype(np.int32), 1, 365)
    sos_av_i = np.clip(np.rint(sos_av).astype(np.int32), 1, 365)
    pos_av_i = np.clip(np.rint(pos_av).astype(np.int32), 1, 365)

    valid_base = (
        _is_valid_value(trc_test, nodata_trc_test) &
        (trc_av != NODATA_OUT) &
        _is_valid_value(sos_y, nodata_sos_y) &
        _is_valid_value(pos_y, nodata_pos_y) &
        _is_valid_value(sos_av, nodata_sos) &
        _is_valid_value(pos_av, nodata_pos) &
        (sos_y > 0) & (sos_y <= 365) &
        (pos_y > 0) & (pos_y <= 365) &
        (sos_av > 0) & (sos_av <= 365) &
        (pos_av > 0) & (pos_av <= 365) &
        (pos_y >= sos_y) &
        (pos_av >= sos_av) &
        mask &
        tr_valid
    )

    (trc_y_actual, _, _, _, _, _, tr_valid_y, trc_y_count, trc_y_fixed_count) = accumulate_daily_tr_sums(
        test_year, profile, valid_base, sos_y_i, pos_y_i, sos_av_i, pos_av_i
    )

    window_len_y = (pos_y_i - sos_y_i + 1).astype(np.int16)
    window_len_av = (pos_av_i - sos_av_i + 1).astype(np.int16)
    min_required_y = np.ceil(MIN_VALID_FRAC * window_len_y).astype(np.int16)
    min_required_av = np.ceil(MIN_VALID_FRAC * window_len_av).astype(np.int16)
    valid_count = (trc_y_count >= min_required_y) & (trc_y_fixed_count >= min_required_av)

    # 额外校验：TRc_{year}.tif 与 当年日序列累计的一致性
    valid_trc_compare = (
        _is_valid_value(trc_test, nodata_trc_test) &
        tr_valid_y &
        valid_count &
        np.isfinite(trc_y_actual) &
        (trc_test != NODATA_OUT)
    )
    if np.any(valid_trc_compare):
        trc_diff = trc_test[valid_trc_compare] - trc_y_actual[valid_trc_compare]
        print("\n  TRc一致性检查 (TRc_file - TRc_daily_actual):")
        print(f"    Mean: {np.mean(trc_diff):.6f} mm")
        print(f"    Std:  {np.std(trc_diff):.6f} mm")
        print(f"    Max:  {np.max(np.abs(trc_diff)):.6f} mm")
        print(f"    95th percentile: {np.percentile(np.abs(trc_diff), 95):.6f} mm")

    valid_check = (
        mask &
        _is_valid_value(sos_y, nodata_sos_y) &
        _is_valid_value(pos_y, nodata_pos_y) &
        (sos_y > 0) & (sos_y <= 366) &
        (pos_y > 0) & (pos_y <= 366) &
        (pos_y >= sos_y) &
        tr_valid_y &
        valid_count &
        np.isfinite(trc_y_actual) &
        (trc_av != NODATA_OUT) &
        (tr_window_test != NODATA_OUT) &
        (tr_fixed_test != NODATA_OUT) &
        (tr_sos_test != NODATA_OUT) &
        (tr_pos_test != NODATA_OUT)
    )

    # 残差检验: TRc_y_actual - TRc_av - TR_window_change - TR_fixed_window
    residual = (trc_y_actual[valid_check] - trc_av[valid_check] -
                tr_window_test[valid_check] - tr_fixed_test[valid_check])
    window_balance = (tr_window_test[valid_check] -
                      (tr_sos_test[valid_check] + tr_pos_test[valid_check]))
    residual_split = (trc_y_actual[valid_check] - trc_av[valid_check] -
                      tr_fixed_test[valid_check] -
                      tr_sos_test[valid_check] - tr_pos_test[valid_check])

    print(f"\n  Test year: {test_year}")
    print(f"  Residual (TRc - TRc_av - TR_window_change - TR_fixed_window):")
    print(f"    Mean: {np.mean(residual):.6f} mm")
    print(f"    Std:  {np.std(residual):.6f} mm")
    print(f"    Max:  {np.max(np.abs(residual)):.6f} mm")
    print(f"    95th percentile: {np.percentile(np.abs(residual), 95):.6f} mm")
    print(f"    99th percentile: {np.percentile(np.abs(residual), 99):.6f} mm")
    print(f"  Window balance (TR_window_change - (TR_sos_change + TR_pos_change)):")
    print(f"    Mean: {np.mean(window_balance):.6f} mm")
    print(f"    Std:  {np.std(window_balance):.6f} mm")
    print(f"    Max:  {np.max(np.abs(window_balance)):.6f} mm")
    print(f"    95th percentile: {np.percentile(np.abs(window_balance), 95):.6f} mm")
    print(f"    99th percentile: {np.percentile(np.abs(window_balance), 99):.6f} mm")
    print(f"  Residual (TRc - TRc_av - TR_fixed_window - TR_sos_change - TR_pos_change):")
    print(f"    Mean: {np.mean(residual_split):.6f} mm")
    print(f"    Std:  {np.std(residual_split):.6f} mm")
    print(f"    Max:  {np.max(np.abs(residual_split)):.6f} mm")
    print(f"    95th percentile: {np.percentile(np.abs(residual_split), 95):.6f} mm")
    print(f"    99th percentile: {np.percentile(np.abs(residual_split), 99):.6f} mm")

    if np.max(np.abs(residual)) < 1e-2:
        print("  ✓ Quality check passed! Residual negligible.")
    else:
        print("  ⚠ Warning: Large residual detected, please check calculation.")

    print("\n✓ Fixed-Window decomposition complete.")
    print(f"Output: {OUTPUT_DIR}")
    print("\n输出文件说明：")
    print("  - TRc_av.tif: 固定窗口基线")
    print("  - Fixed_Window_Length.tif: 固定窗口长度（天数）")
    print("  - TR_fixed_window_YYYY.tif: 固定窗口内的强度差异（核心指标）")
    print("  - TR_window_change_YYYY.tif: 窗口变化带来的实际累积")
    print("  - TR_sos_change_YYYY.tif: SOS变化的贡献")
    print("  - TR_pos_change_YYYY.tif: POS变化的贡献")
    print("\n速率计算：")
    print("  Fixed_Trate = TR_fixed_window / Fixed_Window_Length")


if __name__ == "__main__":
    main()
