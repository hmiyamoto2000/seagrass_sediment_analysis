#python 3.10 コード.py NMDS_scores.csv

#name,NMDS1,NMDS2,cnd,env
#sample1,-0.5,-0.3,seagrass,Saiki
#sample2,-0.4,-0.2,non,Saiki
#...

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Ellipse
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np

# -----------------------------
# 引数チェック
# -----------------------------
if len(sys.argv) != 2:
    print("Usage: python コード.py <input_file.csv>")
    sys.exit(1)

input_csv = sys.argv[1]
output_dir = "NMDS_output"
os.makedirs(output_dir, exist_ok=True)

# -----------------------------
# データ読み込み
# -----------------------------
nmds_scores = pd.read_csv(input_csv)

# 列名チェックと修正（必要に応じて）
expected_columns = {"name", "NMDS1", "NMDS2", "cnd", "env"}
if not expected_columns.issubset(set(nmds_scores.columns)):
    print(f"Error: CSV must contain the columns: {expected_columns}")
    sys.exit(1)

# `group` 列として env を使う（既存コードの流用のため）
nmds_scores["group"] = nmds_scores["env"]

# -----------------------------
# 楕円描画関数
# -----------------------------
def draw_ellipse(position, covariance, ax, color, alpha=0.25):
    U, s, Vt = np.linalg.svd(covariance)
    orient = np.arctan2(U[1, 0], U[0, 0]) * 180 / np.pi
    width, height = 2 * np.sqrt(s)
    ellipse = Ellipse(position, width, height, angle=orient,
                      facecolor=color, alpha=alpha,
                      edgecolor='black', linewidth=1, zorder=3)
    ax.add_patch(ellipse)

# -----------------------------
# NMDS_plot_cnd（色 + 形状）
# -----------------------------
plt.figure(figsize=(8, 6))
sns.set(style="whitegrid")
plt.rcParams["font.family"] = "Helvetica"

cnd_colors = {"seagrass": "green", "non": "blue"}
cnd_markers = {"seagrass": "o", "non": "s"}

unique_cnd = nmds_scores['cnd'].dropna().unique()
ax = plt.gca()

for cnd_val in unique_cnd:
    subset = nmds_scores[nmds_scores['cnd'] == cnd_val]
    color = cnd_colors.get(cnd_val.lower(), "gray")
    marker = cnd_markers.get(cnd_val.lower(), "o")
    ax.scatter(subset["NMDS1"], subset["NMDS2"],
               color=color,
               edgecolor='black',
               s=100,
               alpha=0.8,
               linewidth=1,
               marker=marker,
               zorder=2)

for cnd_val in unique_cnd:
    subset = nmds_scores[nmds_scores['cnd'] == cnd_val]
    if len(subset) < 3:
        continue
    cov = np.cov(subset[["NMDS1", "NMDS2"]].values.T)
    mean = subset[["NMDS1", "NMDS2"]].mean().values
    color = cnd_colors.get(cnd_val.lower(), "gray")
    draw_ellipse(mean, cov, ax, color=color)

plt.xlabel("NMDS1")
plt.ylabel("NMDS2")
plt.title("NMDS Plot colored by cnd")

# 凡例を右外に配置
handles_cnd = [
    mlines.Line2D([], [], marker=cnd_markers.get(k.lower(), "o"), color='black',
                  markerfacecolor=cnd_colors.get(k.lower(), "gray"), linestyle='None',
                  markersize=8, label=k, markeredgewidth=1)
    for k in unique_cnd
]
plt.legend(handles=handles_cnd, title="CND",
           bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

plt.subplots_adjust(right=0.75)  # 右側スペース確保
plt.savefig(os.path.join(output_dir, "NMDS_plot_cnd.pdf"), dpi=300)
plt.savefig(os.path.join(output_dir, "NMDS_plot_cnd.png"), dpi=300)
plt.close()

# -----------------------------
# NMDS_plot_env_cnd（色: env, 形状: cnd）
# -----------------------------
plt.figure(figsize=(8, 6))
plt.rcParams["font.family"] = "Helvetica"

# envをアルファベット順にソート
unique_env = sorted(nmds_scores['group'].dropna().unique())
palette_env = sns.color_palette("Set2", len(unique_env))
color_dict_env = dict(zip(unique_env, palette_env))

def get_marker(cnd_val):
    if isinstance(cnd_val, str):
        if cnd_val.lower() == 'seagrass':
            return 'o'
        elif cnd_val.lower() == 'non':
            return 's'
    return 'o'

ax = plt.gca()

for env_val in unique_env:
    subset = nmds_scores[nmds_scores['group'] == env_val]
    base_color = color_dict_env[env_val]
    for _, row in subset.iterrows():
        marker = get_marker(row['cnd'])
        ax.scatter(row['NMDS1'], row['NMDS2'],
                   color=base_color,
                   edgecolor='black',
                   s=100,
                   marker=marker,
                   linewidth=1,
                   zorder=2)

    if len(subset) >= 3:
        cov = np.cov(subset[["NMDS1", "NMDS2"]].values.T)
        mean = subset[["NMDS1", "NMDS2"]].mean().values
        draw_ellipse(mean, cov, ax, color=base_color)

plt.xlabel("NMDS1")
plt.ylabel("NMDS2")
plt.title("NMDS Plot colored by env (region) with cnd markers")

# env凡例（色）を右上に配置（ABC順）
handles_env = [mpatches.Patch(color=color_dict_env[env], label=env) for env in unique_env]
legend_env = ax.legend(handles=handles_env, title="Env (region)",
                       bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

# cnd凡例（形状）をその下に配置
handles_cnd = [
    mlines.Line2D([], [], marker=get_marker(k), color='black',
                  markerfacecolor='white', linestyle='None',
                  markersize=8, label=k, markeredgewidth=1)
    for k in cnd_colors.keys()
]

ax.add_artist(legend_env)
ax.legend(handles=handles_cnd, title="CND (marker shape)",
          bbox_to_anchor=(1.05, 0.5), loc='upper left', borderaxespad=0.)

plt.subplots_adjust(right=0.75)
plt.savefig(os.path.join(output_dir, "NMDS_plot_env_cnd.pdf"), dpi=300)
plt.savefig(os.path.join(output_dir, "NMDS_plot_env_cnd.png"), dpi=300)
plt.close()

print("グラフ出力: NMDS_output ")
