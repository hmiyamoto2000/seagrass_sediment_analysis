#完全コード

#このコードを例として shap_compare.py などに保存すれば、コマンドラインで python3.10 shap_compare.py と実行するだけで、Random Forest、XGBoost、LightGBM の SHAP を一つのグラフで比較できます。

#target,feature1,feature2,feature3
#0,10.5,3.2,1
#1,8.1,5.6,0
#0,9.0,2.1,1
#...
#1列目がtarget（0 or 1）

#python3.10 shap_compare.py your_data.csv

# 修正済み shap_compare_fixed.py

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import shap
import re

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
import xgboost as xgb
import lightgbm as lgb
from sklearn.metrics import accuracy_score
import matplotlib

# ======== グローバルフォント設定 ========
if "Helvetica" in matplotlib.font_manager.findSystemFonts(fontpaths=None, fontext='ttf'):
    matplotlib.rcParams['font.family'] = 'Helvetica'
else:
    matplotlib.rcParams['font.family'] = 'Arial'

# ================= SHAP計算 =================
def compute_shap_values(model, X, model_name):
    print(f"Computing SHAP values for {model_name}...")
    explainer = shap.Explainer(model, X)
    shap_values = explainer(X)
    if len(shap_values.values.shape) == 3:
        shap_vals = shap_values.values[:, :, 1]
    else:
        shap_vals = shap_values.values
    return shap_vals


# ================ SHAPサマリープロット =================
def plot_shap_summary(shap_vals, X_test, model_name, output_dir):
    plt.figure(figsize=(10, 6))
    shap.summary_plot(shap_vals, X_test, show=False)
    plt.xlabel("SHAP value")
    plt.grid(True, axis='x', linestyle='--', alpha=0.7)
    png_path = os.path.join(output_dir, f"shap_summary_{model_name}.png")
    pdf_path = os.path.join(output_dir, f"shap_summary_{model_name}.pdf")
    plt.savefig(png_path, bbox_inches='tight', dpi=300)
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(10, 6))
    shap.summary_plot(shap_vals, X_test, show=False, plot_type="bar")
    plt.xlabel("SHAP value")
    plt.grid(True, axis='x', linestyle='--', alpha=0.7)
    no_label_png = os.path.join(output_dir, f"shap_summary_{model_name}_no_label.png")
    no_label_pdf = os.path.join(output_dir, f"shap_summary_{model_name}_no_label.pdf")
    plt.savefig(no_label_png, bbox_inches='tight', dpi=300)
    plt.savefig(no_label_pdf, bbox_inches='tight')
    plt.close()

    print(f"Saved SHAP summary plots for {model_name} at:\n  {png_path}\n  {pdf_path}\n  {no_label_png}\n  {no_label_pdf}")


# ================ SHAP重要度保存 =================
def save_shap_importance(shap_vals, X, model_name, output_dir):
    importance = np.abs(shap_vals).mean(axis=0)
    df = pd.DataFrame({
        'Feature': X.columns,
        'Importance': importance
    }).sort_values(by='Importance', ascending=False)
    csv_path = os.path.join(output_dir, f"{model_name}_SHAP.csv")
    df.to_csv(csv_path, index=False)
    print(f"Saved SHAP importance CSV for {model_name} at:\n  {csv_path}")
    return csv_path


# ================ バブルチャート =================
def plot_shap_bubble_chart(rf_csv, xgb_csv, lgbm_csv, output_dir):
    df_rf = pd.read_csv(rf_csv).rename(columns={'Importance': 'SHAP_RF'})
    df_xgb = pd.read_csv(xgb_csv).rename(columns={'Importance': 'SHAP_XGB'})
    df_lgbm = pd.read_csv(lgbm_csv).rename(columns={'Importance': 'SHAP_LGBM'})

    df_merge = df_rf[['Feature', 'SHAP_RF']].merge(
        df_xgb[['Feature', 'SHAP_XGB']], on='Feature', how='inner').merge(
        df_lgbm[['Feature', 'SHAP_LGBM']], on='Feature', how='inner')

    x_con = df_merge['SHAP_XGB'].values
    y_con = df_merge['SHAP_RF'].values
    z_con = df_merge['SHAP_LGBM'].values
    labels = df_merge['Feature'].values

    sizes = np.sqrt(z_con) * 100
    fixed_size = 100

    fig, ax = plt.subplots(figsize=(12, 9))
    cm = plt.get_cmap('viridis')
    scat = ax.scatter(x_con, y_con, s=sizes, c=z_con,
                      alpha=0.7, cmap=cm,
                      edgecolors='black', linewidth=0.5)
    cbar = fig.colorbar(scat, ax=ax)
    cbar.set_label('SHAP LightGBM')
    ax.set_xlabel('SHAP XGBoost')
    ax.set_ylabel('SHAP Random Forest')
    ax.set_title('SHAP Importance Comparison: XGB vs RF vs LGBM (Bubble Size Varying)')

    for i, label in enumerate(labels):
        ax.annotate(label, (x_con[i], y_con[i]), fontsize=8, alpha=0.7)

    plt.tight_layout()
    png_path = os.path.join(output_dir, "SHAP_bubble_chart.png")
    pdf_path = os.path.join(output_dir, "SHAP_bubble_chart.pdf")
    fig.savefig(png_path, bbox_inches='tight', dpi=300)
    fig.savefig(pdf_path, bbox_inches='tight')
    plt.close()

    # ===== 固定サイズバブルチャート =====
    fig2, ax2 = plt.subplots(figsize=(12, 9))
    scat2 = ax2.scatter(x_con, y_con, s=fixed_size, c=z_con,
                        alpha=0.7, cmap=cm,
                        edgecolors='black', linewidth=0.5)
    cbar2 = fig2.colorbar(scat2, ax=ax2)
    cbar2.set_label('SHAP LightGBM')
    ax2.set_xlabel('SHAP XGBoost')
    ax2.set_ylabel('SHAP Random Forest')
    ax2.set_title('SHAP Importance Comparison: XGB vs RF vs LGBM (Fixed Bubble Size)')

    for i, label in enumerate(labels):
        ax2.annotate(label, (x_con[i], y_con[i]), fontsize=8, alpha=0.7)

    plt.tight_layout()
    png_path_fixed = os.path.join(output_dir, "SHAP_bubble_chart_fixed_size.png")
    pdf_path_fixed = os.path.join(output_dir, "SHAP_bubble_chart_fixed_size.pdf")
    fig2.savefig(png_path_fixed, bbox_inches='tight', dpi=300)
    fig2.savefig(pdf_path_fixed, bbox_inches='tight')
    plt.close()

    return df_merge


# ============ LightGBMサイズ・色変化バブルチャート ============
def plot_lightgbm_variable_bubble_chart(merged_df, output_dir):
    def normalize_colname(name):
        name = name.lower()
        name = re.sub(r'[^a-z0-9]+', '_', name)
        name = re.sub(r'_+', '_', name)
        return name.strip('_')

    def find_column(df, candidates):
        norm_cols = {normalize_colname(col): col for col in df.columns}
        for c in candidates:
            nc = normalize_colname(c)
            if nc in norm_cols:
                return norm_cols[nc]
        raise ValueError(f"None of the candidate columns found: {candidates}")

    xg_col = find_column(merged_df, ['shap_xgb'])
    rf_col = find_column(merged_df, ['shap_rf'])
    lgb_col = find_column(merged_df, ['shap_lgbm'])

    x = merged_df[xg_col]
    y = merged_df[rf_col]
    labels = merged_df['Feature']
    lightgbm_raw = merged_df[lgb_col]

    if (lightgbm_raw < 0).any():
        print("Warning: negative values in LightGBM SHAP detected. Setting them to 0.")
        lightgbm_raw = lightgbm_raw.clip(lower=0)

    min_size = 50
    max_size = 1000
    if lightgbm_raw.max() == lightgbm_raw.min():
        size = np.full_like(lightgbm_raw, (min_size + max_size) / 2)
    else:
        norm_sizes = (lightgbm_raw - lightgbm_raw.min()) / (lightgbm_raw.max() - lightgbm_raw.min())
        size = norm_sizes * (max_size - min_size) + min_size

    sorted_idx = lightgbm_raw.argsort()[::-1]
    x_sorted = x.iloc[sorted_idx]
    y_sorted = y.iloc[sorted_idx]
    labels_sorted = labels.iloc[sorted_idx]
    lightgbm_sorted = lightgbm_raw.iloc[sorted_idx]
    size_sorted = pd.Series(size).iloc[sorted_idx]

    colormap = plt.colormaps['rainbow']

    # ========== ラベル付き ==========
    fig, ax = plt.subplots(figsize=(12, 9))
    scatter = ax.scatter(x_sorted, y_sorted, s=size_sorted,
                         c=lightgbm_sorted, cmap=colormap, alpha=0.85,
                         edgecolors='black', linewidths=0.7)

    for i, txt in enumerate(labels_sorted):
        ax.annotate(
            txt,
            (x_sorted.iloc[i], y_sorted.iloc[i]),
            fontsize=8, alpha=0.7,
            textcoords="offset points", xytext=(5, 3),
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="gray", lw=0.5, alpha=0.6)
        )

    cbar = fig.colorbar(scatter)
    cbar.set_label("LightGBM SHAP Importance (color/size)")
    ax.set_xlabel("XGBoost SHAP Importance")
    ax.set_ylabel("Random Forest SHAP Importance")

    twin_ax = ax.twinx()
    twin_ax.tick_params(
        labelbottom=False, labelleft=False,
        labelright=False, labeltop=False,
        bottom=False, left=False, right=False, top=False
    )
    twin_ax.set_ylabel("LightGBM SHAP", labelpad=110)

    plt.title("SHAP Bubble Chart (LightGBM variable size/color + feature names)")
    plt.tight_layout()

    output_path_png = os.path.join(output_dir, "SHAP_bubble_chart_lightgbm_variable.png")
    output_path_pdf = os.path.join(output_dir, "SHAP_bubble_chart_lightgbm_variable.pdf")
    fig.savefig(output_path_png, bbox_inches="tight", dpi=300)
    fig.savefig(output_path_pdf, bbox_inches="tight")
    plt.close()

    # ========== ラベル無し ==========
    fig2, ax2 = plt.subplots(figsize=(12, 9))
    scatter2 = ax2.scatter(x_sorted, y_sorted, s=size_sorted,
                           c=lightgbm_sorted, cmap=colormap, alpha=0.85,
                           edgecolors='black', linewidths=0.7)

    cbar2 = fig2.colorbar(scatter2)
    cbar2.set_label("LightGBM SHAP Importance (color/size)")
    ax2.set_xlabel("XGBoost SHAP Importance")
    ax2.set_ylabel("Random Forest SHAP Importance")

    twin_ax2 = ax2.twinx()
    twin_ax2.tick_params(labelbottom=False, labelleft=False, labelright=False,
                         labeltop=False, bottom=False, left=False, right=False, top=False)
    twin_ax2.set_ylabel("LightGBM SHAP", labelpad=110)

    plt.title("SHAP Bubble Chart (No Labels)")
    plt.tight_layout()

    no_label_png = os.path.join(output_dir, "SHAP_bubble_chart_lightgbm_variable_no_labels.png")
    no_label_pdf = os.path.join(output_dir, "SHAP_bubble_chart_lightgbm_variable_no_labels.pdf")
    fig2.savefig(no_label_png, bbox_inches="tight", dpi=300)
    fig2.savefig(no_label_pdf, bbox_inches="tight")
    plt.close()

    print(f"Saved labeled and no-label bubble charts to:\n  {output_path_png}\n  {no_label_png}")


# ============== メイン関数 ==============
def main():
    if len(sys.argv) < 2:
        print("Usage: python shap_compare_fixed.py path_to_your_data.csv")
        sys.exit(1)

    csv_path = sys.argv[1]
    df = pd.read_csv(csv_path)

    print("Data info:")
    print(df.info())
    print("Head:")
    print(df.head())

    target_col = "target"
    if target_col not in df.columns:
        print(f"Error: target column '{target_col}' not found in data")
        sys.exit(1)

    X = df.drop(columns=[target_col])
    y = df[target_col]

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.25, random_state=42, stratify=y)

    output_dir = "Result_SHAP"
    os.makedirs(output_dir, exist_ok=True)

    # === Random Forest ===
    rf = RandomForestClassifier(random_state=42)
    rf.fit(X_train, y_train)
    print(f"RF Test accuracy: {accuracy_score(y_test, rf.predict(X_test)):.4f}")
    shap_vals_rf = compute_shap_values(rf, X_test, "Random_Forest")
    plot_shap_summary(shap_vals_rf, X_test, "Random_Forest", output_dir)
    rf_csv = save_shap_importance(shap_vals_rf, X_test, "Random_Forest", output_dir)

    # === XGBoost ===
    xgb_model = xgb.XGBClassifier(use_label_encoder=False, eval_metric="logloss", random_state=42)
    xgb_model.fit(X_train, y_train)
    print(f"XGB Test accuracy: {accuracy_score(y_test, xgb_model.predict(X_test)):.4f}")
    shap_vals_xgb = compute_shap_values(xgb_model, X_test, "XGBoost")
    plot_shap_summary(shap_vals_xgb, X_test, "XGBoost", output_dir)
    xgb_csv = save_shap_importance(shap_vals_xgb, X_test, "XGBoost", output_dir)

    # === LightGBM ===
    lgb_model = lgb.LGBMClassifier(random_state=42)
    lgb_model.fit(X_train, y_train)
    print(f"LGBM Test accuracy: {accuracy_score(y_test, lgb_model.predict(X_test)):.4f}")
    shap_vals_lgb = compute_shap_values(lgb_model, X_test, "LightGBM")
    plot_shap_summary(shap_vals_lgb, X_test, "LightGBM", output_dir)
    lgbm_csv = save_shap_importance(shap_vals_lgb, X_test, "LightGBM", output_dir)

    # === バブルチャート作成 ===
    merged_df = plot_shap_bubble_chart(rf_csv, xgb_csv, lgbm_csv, output_dir)
    plot_lightgbm_variable_bubble_chart(merged_df, output_dir)


if __name__ == "__main__":
    main()
