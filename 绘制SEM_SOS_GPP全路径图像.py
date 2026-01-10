import csv
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from matplotlib.patches import FancyArrowPatch

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['axes.unicode_minus'] = False

# === 读取05b结果数据（固定窗口，双时间尺度）===
RESULT_DIR_WIN = Path(r"D:\claude-project\结果分析")
RESULT_DIR_WSL = Path("/mnt/d/claude-project/结果分析")
RESULT_DIR = RESULT_DIR_WIN if RESULT_DIR_WIN.exists() else RESULT_DIR_WSL

PARAM_FILE = RESULT_DIR / "SEM_dual_timescale_parameters.csv"
R2_FILE = RESULT_DIR / "SEM_dual_timescale_R2.csv"

def load_sem_params(param_file):
    coef_map = {}
    if not param_file.exists():
        raise FileNotFoundError(f"参数文件不存在: {param_file}")
    with param_file.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("op") != "~":
                continue
            lhs = row.get("lhs")
            rhs = row.get("rhs")
            std_all = row.get("std.all")
            p_val = row.get("pvalue")
            if lhs is None or rhs is None:
                continue
            try:
                coef = float(std_all)
                pval = float(p_val)
            except (TypeError, ValueError):
                continue
            coef_map[(lhs, rhs)] = {"coef": coef, "p": pval}
    return coef_map

def load_r2(r2_file):
    r2_map = {}
    if not r2_file.exists():
        raise FileNotFoundError(f"R2文件不存在: {r2_file}")
    with r2_file.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f)
        _ = next(reader, None)
        for row in reader:
            if len(row) < 2:
                continue
            var = row[0].strip().strip('"')
            try:
                val = float(row[1])
            except ValueError:
                continue
            r2_map[var] = val
    return r2_map

def build_path(coef_map, lhs, rhs):
    item = coef_map.get((lhs, rhs))
    if not item:
        return None
    coef = item["coef"]
    pval = item["p"]
    if not np.isfinite(coef) or not np.isfinite(pval):
        return None
    return {"coef": coef, "sig": pval < 0.05}

coef_map = load_sem_params(PARAM_FILE)
r2_map = load_r2(R2_FILE)

# R² 值（变量概念调整：GPP_season / Fixed_Trate）
r2_values = {
    "SOS": r2_map.get("SOS", np.nan),
    "GPP": r2_map.get("GPP_season", np.nan),
    "TRATE": r2_map.get("Fixed_Trate", np.nan),
}

display_labels = {
    "TRATE": r"T$_{rate}$",
}

# 第一层：季前气候 → SOS
pre_season_to_sos = {}
for key, (lhs, rhs) in {
    "sos~p_pre": ("SOS", "P_pre"),
    "sos~t_pre": ("SOS", "T_pre"),
    "sos~sw_pre": ("SOS", "SW_pre"),
}.items():
    path = build_path(coef_map, lhs, rhs)
    if path:
        pre_season_to_sos[key] = path

# 第二层：SOS + 季前气候 + 季内气候 → GPP
sos_and_pre_and_season_to_gpp = {}
for key, (lhs, rhs) in {
    "gpp~sos": ("GPP_season", "SOS"),
    "gpp~p_pre": ("GPP_season", "P_pre"),
    "gpp~t_pre": ("GPP_season", "T_pre"),
    "gpp~sw_pre": ("GPP_season", "SW_pre"),
    "gpp~p_season": ("GPP_season", "P_season"),
    "gpp~t_season": ("GPP_season", "T_season"),
    "gpp~sw_season": ("GPP_season", "SW_season"),
}.items():
    path = build_path(coef_map, lhs, rhs)
    if path:
        sos_and_pre_and_season_to_gpp[key] = path

# 第三层：SOS + GPP + 季前气候 + 季内气候 → Fixed_Trate
sos_gpp_pre_and_season_to_trate = {}
for key, (lhs, rhs) in {
    "trate~sos": ("Fixed_Trate", "SOS"),
    "trate~gpp": ("Fixed_Trate", "GPP_season"),
    "trate~p_pre": ("Fixed_Trate", "P_pre"),
    "trate~t_pre": ("Fixed_Trate", "T_pre"),
    "trate~sw_pre": ("Fixed_Trate", "SW_pre"),
    "trate~p_season": ("Fixed_Trate", "P_season"),
    "trate~t_season": ("Fixed_Trate", "T_season"),
    "trate~sw_season": ("Fixed_Trate", "SW_season"),
}.items():
    path = build_path(coef_map, lhs, rhs)
    if path:
        sos_gpp_pre_and_season_to_trate[key] = path

# 合并所有直接路径
direct_paths = {
    **pre_season_to_sos,
    **sos_and_pre_and_season_to_gpp,
    **sos_gpp_pre_and_season_to_trate,
}

# === 创建图形 ===
fig = plt.figure(figsize=(11, 12))
fig.patch.set_facecolor('white')

ax1 = plt.subplot2grid((4, 1), (0, 0), rowspan=3)
ax2 = plt.subplot2grid((4, 1), (3, 0))

ax1.set_facecolor('white')
ax2.set_facecolor('#F8F8F8')

ax1.set_xlim(0, 11)
ax1.set_ylim(0, 12)

# 节点位置
positions = {
    'P_PRE': (1.5, 9.5),
    'T_PRE': (5.5, 9.5),
    'SW_PRE': (9.5, 9.5),

    'SOS': (1.5, 5.5),
    'GPP': (5.5, 5.5),
    'TRATE': (9.5, 5.5),

    'P_SEASON': (1.5, 1.5),
    'T_SEASON': (5.5, 1.5),
    'SW_SEASON': (9.5, 1.5),
}

var_to_pos = {
    'p_pre': 'P_PRE', 't_pre': 'T_PRE', 'sw_pre': 'SW_PRE',
    'sos': 'SOS', 'gpp': 'GPP', 'trate': 'TRATE',
    'p_season': 'P_SEASON', 't_season': 'T_SEASON', 'sw_season': 'SW_SEASON'
}

def draw_box(ax, name, pos, is_outcome=False):
    x, y = pos
    if is_outcome:
        width, height = 2.2, 1.4
        fontsize, r2_fontsize = 24, 18
    else:
        width, height = 2.0, 1.2
        fontsize, r2_fontsize = 20, 16

    shadow = patches.Rectangle(
        (x - width/2 + 0.03, y - height/2 - 0.03), width, height,
        facecolor='lightgray', alpha=0.3, zorder=1
    )
    ax.add_patch(shadow)

    box = patches.Rectangle(
        (x - width/2, y - height/2), width, height,
        facecolor='#E8F5E9' if is_outcome else '#F5F5F5',
        edgecolor='#4CAF50' if is_outcome else '#757575',
        linewidth=2.5 if is_outcome else 2, zorder=2
    )
    ax.add_patch(box)

    label = display_labels.get(name, name)
    if name in r2_values and np.isfinite(r2_values[name]):
        ax.text(x, y + 0.15, label, ha='center', va='center',
                fontsize=fontsize, fontweight='bold', color='black')
        ax.text(x, y - 0.25, f'R² = {r2_values[name]:.3f}',
                ha='center', va='center',
                fontsize=r2_fontsize, color='black', fontweight='bold')
    else:
        ax.text(x, y, label, ha='center', va='center',
                fontsize=fontsize, fontweight='bold', color='black')

def draw_arrow(ax, start_pos, end_pos, coef, significant, start_name, end_name, is_sos_to_et=False):
    x1, y1 = start_pos
    x2, y2 = end_pos

    if start_name in ['SOS', 'GPP', 'TRATE']:
        start_width, start_height = 2.2, 1.4
    else:
        start_width, start_height = 2.0, 1.2

    if end_name in ['SOS', 'GPP', 'TRATE']:
        end_width, end_height = 2.2, 1.4
    else:
        end_width, end_height = 2.0, 1.2

    dx = x2 - x1
    dy = y2 - y1

    # 起点调整
    if start_name in ['P_PRE', 'T_PRE', 'SW_PRE']:
        start_x, start_y = x1, y1 - start_height/2
    elif start_name in ['P_SEASON', 'T_SEASON', 'SW_SEASON']:
        start_x, start_y = x1, y1 + start_height/2
    else:
        if end_name in ['P_SEASON', 'T_SEASON', 'SW_SEASON']:
            start_x, start_y = x1, y1 - start_height/2
        elif end_name in ['P_PRE', 'T_PRE', 'SW_PRE']:
            start_x, start_y = x1, y1 + start_height/2
        else:
            if dx > 0:
                start_x, start_y = x1 + start_width/2, y1
            else:
                start_x, start_y = x1 - start_width/2, y1

    # 终点调整
    if end_name in ['P_PRE', 'T_PRE', 'SW_PRE']:
        end_x, end_y = x2, y2 - end_height/2
    elif end_name in ['P_SEASON', 'T_SEASON', 'SW_SEASON']:
        end_x, end_y = x2, y2 + end_height/2
    else:
        if start_name in ['P_PRE', 'T_PRE', 'SW_PRE']:
            end_x, end_y = x2, y2 + end_height/2
        elif start_name in ['P_SEASON', 'T_SEASON', 'SW_SEASON']:
            end_x, end_y = x2, y2 - end_height/2
        else:
            if dx > 0:
                end_x, end_y = x2 - end_width/2, y2
            else:
                end_x, end_y = x2 + end_width/2, y2

    # ✅ 统一颜色：只区分正负，不区分新旧
    if significant:
        color = '#D32F2F' if coef > 0 else '#1976D2'
        linestyle = '-'
        abs_coef = abs(coef)
        # ✅ 重新调整线宽：拉大差异，最粗的更细
        if abs_coef >= 0.4: linewidth = 6
        elif abs_coef >= 0.3: linewidth = 5
        elif abs_coef >= 0.2: linewidth = 4
        elif abs_coef >= 0.1: linewidth = 2.5
        else: linewidth = 1.5
        alpha = 1.0
    else:
        color = '#9E9E9E'
        linestyle = '--'
        linewidth = 2
        alpha = 0.6

    # ✅ 判断路径类型
    is_pre_to_gpp = (end_name == 'GPP' and start_name in ['P_PRE', 'T_PRE', 'SW_PRE'])
    is_season_to_gpp = (end_name == 'GPP' and start_name in ['P_SEASON', 'T_SEASON', 'SW_SEASON'])
    is_sos_to_gpp = (start_name == 'SOS' and end_name == 'GPP')
    is_gpp_to_trate = (start_name == 'GPP' and end_name == 'TRATE')

    # ✅ SOS→ET使用对称向下弯曲的圆滑弧线
    if is_sos_to_et:
        arrow = FancyArrowPatch(
            (start_x, start_y), (end_x, end_y),
            connectionstyle="arc3,rad=0.3",  # ✅ 圆滑对称弧线
            arrowstyle='->',  # 默认箭头样式
            mutation_scale=25,  # ✅ 与直线箭头大小一致
            color=color,
            linewidth=linewidth,
            linestyle=linestyle,
            alpha=alpha,
            zorder=2  # 箭头在底层
        )
        ax.add_patch(arrow)

        # ✅ 标签位置：弧线1/5处，精确计算弧线上的点和切线角度
        t = 0.2  # 1/5位置

        # 直线插值位置
        straight_x = start_x + (end_x - start_x) * t
        straight_y = start_y + (end_y - start_y) * t

        # 弧线向下弯曲的偏移量（更精确的计算）
        # arc3的rad参数定义了弧线的弯曲程度
        # 在参数t处，弧线偏移 ≈ rad * sin(π*t) * 弦长 / 2
        chord_length = ((end_x - start_x)**2 + (end_y - start_y)**2)**0.5
        arc_offset_y = -0.3 * np.sin(np.pi * t) * chord_length * 0.6  # 调整系数使标注在线上

        label_x = straight_x
        label_y = straight_y + arc_offset_y

        # 计算弧线在t处的切线角度
        dx_dt = (end_x - start_x)
        dy_dt = (end_y - start_y) - 0.3 * np.pi * np.cos(np.pi * t) * chord_length * 0.6
        angle = np.arctan2(dy_dt, dx_dt) * 180 / np.pi
        if angle > 90:
            angle = angle - 180
        elif angle < -90:
            angle = angle + 180
    else:
        # 普通直线箭头
        ax.annotate('', xy=(end_x, end_y), xytext=(start_x, start_y),
                    arrowprops=dict(arrowstyle='->', color=color,
                                  linewidth=linewidth, mutation_scale=25,
                                  linestyle=linestyle, alpha=alpha, zorder=2))

        # ✅ 标签位置规则
        if is_pre_to_gpp or is_season_to_gpp:
            # 季前气候→GPP、春季气候→GPP：1/3位置
            label_x = start_x + (end_x - start_x) * 0.33
            label_y = start_y + (end_y - start_y) * 0.33
        elif is_sos_to_gpp or is_gpp_to_trate:
            # SOS→GPP、GPP→TRate：1/2位置
            label_x = start_x + (end_x - start_x) * 0.5
            label_y = start_y + (end_y - start_y) * 0.5
        else:
            # 其他路径（季前→SOS、季前→ET、春季→ET）：2/3位置
            label_x = start_x + (end_x - start_x) * 0.67
            label_y = start_y + (end_y - start_y) * 0.67

        angle = np.arctan2(end_y - start_y, end_x - start_x) * 180 / np.pi
        if angle > 90:
            angle = angle - 180
        elif angle < -90:
            angle = angle + 180

    # ✅ 绘制文本（置于最上层）
    ax.text(label_x, label_y, f'{coef:.2f}' if abs(coef) >= 0.1 else f'{coef:.3f}',
            ha='center', va='center', fontsize=18,
            fontweight='bold', color='black',
            rotation=angle, rotation_mode='anchor',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                     edgecolor='none', alpha=0.95),
            zorder=100)  # ✅ 文本在最上层

# 绘制变量框
for name, pos in positions.items():
    is_outcome = name in ['SOS', 'GPP', 'TRATE']
    draw_box(ax1, name, pos, is_outcome)

# 绘制路径
for path_key, path_data in direct_paths.items():
    parts = path_key.split('~')
    if len(parts) == 2:
        end_var, start_var = parts
        if start_var in var_to_pos and end_var in var_to_pos:
            start_name = var_to_pos[start_var]
            end_name = var_to_pos[end_var]
            is_sos_trate = (start_name == 'SOS' and end_name == 'TRATE')
            draw_arrow(ax1, positions[start_name], positions[end_name],
                      path_data['coef'], path_data['sig'], start_name, end_name, is_sos_trate)

# 标题
ax1.text(5.5, 11.2, 'SOS - Dual-Timescale SEM (Fixed Window, Full Paths)',
        fontsize=22, fontweight='bold', ha='center',
        bbox=dict(boxstyle='round,pad=0.8', facecolor='white',
                 edgecolor='#666666', linewidth=3))

# ✅ 简化图例：不区分新旧路径
legend_elements = [
    plt.Line2D([0], [0], color='#D32F2F', linewidth=6, label='Positive effect'),
    plt.Line2D([0], [0], color='#1976D2', linewidth=6, label='Negative effect'),
    plt.Line2D([0], [0], color='#9E9E9E', linewidth=4, linestyle='--', label='Non-significant')
]
ax1.legend(handles=legend_elements, loc='upper left', fontsize=14, frameon=True,
          fancybox=True, shadow=True, bbox_to_anchor=(0.01, 0.99))

ax1.set_xticks([])
ax1.set_yticks([])
for spine in ax1.spines.values():
    spine.set_visible(False)

# === 下方效应分解图 ===
variables = ['P', 'T', 'SW']
var_names = ['Precipitation', 'Temperature', 'Solar Radiation']

def get_coef(lhs, rhs):
    item = coef_map.get((lhs, rhs))
    return item["coef"] if item else np.nan

e1 = get_coef("Fixed_Trate", "P_season")
e2 = get_coef("Fixed_Trate", "T_season")
e3 = get_coef("Fixed_Trate", "SW_season")
c1 = get_coef("GPP_season", "P_season")
c2 = get_coef("GPP_season", "T_season")
c3 = get_coef("GPP_season", "SW_season")
d_coef = get_coef("Fixed_Trate", "GPP_season")

direct_effects = [e1, e2, e3]
indirect_via_gpp = [c1 * d_coef, c2 * d_coef, c3 * d_coef]
total_effects = [d + i for d, i in zip(direct_effects, indirect_via_gpp)]
mediation_ratios = [
    (i / t) if np.isfinite(t) and t != 0 else np.nan
    for i, t in zip(indirect_via_gpp, total_effects)
]

x_pos = np.arange(len(variables))
width = 0.5

bars1 = ax2.bar(x_pos, direct_effects, width,
               label='Direct (Season Climate → TRate)',
               color='#9C27B0', alpha=0.85, edgecolor='#7B1FA2', linewidth=2)

bars2 = ax2.bar(x_pos, indirect_via_gpp, width, bottom=direct_effects,
               label='Indirect (Season Climate → GPP → TRate)',
               color='#FF9800', alpha=0.85, edgecolor='#F57C00', linewidth=2)

for i, (total, mediation_pct) in enumerate(zip(total_effects, mediation_ratios)):
    pct_text = f"{mediation_pct*100:.1f}%" if np.isfinite(mediation_pct) else "NA"
    ax2.text(x_pos[i], total + 0.02, f'{total:.3f}\n({pct_text})',
            ha='center', va='bottom', fontsize=16, fontweight='bold')

ax2.set_xlim(-0.5, 2.5)
ymin = min(0, np.nanmin(total_effects + direct_effects + indirect_via_gpp)) - 0.05
ymax = max(0, np.nanmax(total_effects + direct_effects + indirect_via_gpp)) + 0.05
ax2.set_ylim(ymin, ymax)
ax2.axhline(y=0, color='black', linewidth=2.5)
ax2.set_ylabel('Standardized Effect on TRate', fontsize=18, fontweight='bold')
ax2.set_xlabel('Season Climate Variables (Fixed Window)', fontsize=18, fontweight='bold')
ax2.set_xticks(x_pos)
ax2.set_xticklabels(var_names, fontsize=16, fontweight='bold')
ax2.tick_params(axis='both', which='major', labelsize=16)

legend = ax2.legend(loc='upper left', fontsize=16, frameon=True, fancybox=True, shadow=True)
legend.get_frame().set_facecolor('white')
legend.get_frame().set_alpha(0.9)

ax2.grid(True, alpha=0.3, axis='y', linestyle='-', linewidth=0.5)
ax2.set_axisbelow(True)
for spine in ['top', 'right']:
    ax2.spines[spine].set_visible(False)
ax2.spines['left'].set_linewidth(2.5)
ax2.spines['bottom'].set_linewidth(2.5)

plt.tight_layout()
output_png = RESULT_DIR / "SEM_model_SOS_GPP_dual_timescale_05b.png"
plt.savefig(output_png, dpi=300, bbox_inches='tight', facecolor='white')

print("✅ SOS双时间尺度SEM可视化已保存!")
print(f"\n模型统计:")
print(f"   • R² 值: SOS={r2_values['SOS']:.3f}, GPP={r2_values['GPP']:.3f}, TRate={r2_values['TRATE']:.3f}")
if ("Fixed_Trate", "SOS") in coef_map and ("Fixed_Trate", "GPP_season") in coef_map and ("GPP_season", "SOS") in coef_map:
    g_coef = get_coef("Fixed_Trate", "SOS")
    b_coef = get_coef("GPP_season", "SOS")
    d_coef = get_coef("Fixed_Trate", "GPP_season")
    total_eff = g_coef + b_coef * d_coef
    print(f"   • SOS → TRate 总效应: {total_eff:.3f} (直接{g_coef:.3f} + 间接{(b_coef*d_coef):.3f})")
