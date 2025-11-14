# causal_analysis_X_with_tikz.py
# causal_analysis_X_with_tikz.py
import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from causallearn.search.ConstraintBased.PC import pc
from scipy.stats import shapiro
from causallearn.utils.cit import kci

def detect_discrete_columns(df, threshold=5):
    return [col for col in df.columns if df[col].nunique() <= threshold]

def save_tikz(graph: nx.DiGraph, var_names, filepath):
    import numpy as np

    lines = []
    lines.append("\\begin{tikzpicture}[->,>=Stealth, thick, main node/.style={circle,draw,fill=yellow!20, font=\\sffamily\\bfseries, minimum size=1cm, inner sep=2pt}]")

    # NetworkXのspring_layoutで座標取得
    pos = nx.spring_layout(graph, seed=42)

    # posは{node: (x,y)}、x,yはだいたい[-1,1]の範囲なので
    # 横幅を約10cmにスケーリングしつつyは上下反転（TikZの座標系合わせ）
    xs = np.array([pos[node][0] for node in var_names])
    ys = np.array([pos[node][1] for node in var_names])

    x_min, x_max = xs.min(), xs.max()
    if x_max != x_min:
        xs_scaled = 10 * (xs - x_min) / (x_max - x_min)
    else:
        xs_scaled = np.zeros_like(xs)

    y_min, y_max = ys.min(), ys.max()
    if y_max != y_min:
        ys_scaled = -10 * (ys - y_min) / (y_max - y_min)
    else:
        ys_scaled = np.zeros_like(ys)

    # ノード定義
    for i, v in enumerate(var_names):
        x = xs_scaled[i]
        y = ys_scaled[i]
        lines.append(f"\\node[main node] ({v}) at ({x:.2f},{y:.2f}) {{{v}}};")

    # エッジ定義（bend left=20で矢印を少し曲げる）
    for source, target in graph.edges():
        lines.append(f"\\path[->, bend left=20] ({source}) edge ({target});")

    lines.append("\\end{tikzpicture}")

    with open(filepath, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"TikZコードを保存しました: {filepath}")

def main(csv_path):
    print(f"\nCSVファイルを読み込みます: {csv_path}")
    data = pd.read_csv(csv_path)
    n_vars = data.shape[1]

    # 元の変数名と置換
    original_var_names = data.columns.tolist()
    var_names = [f"X{i+1}" for i in range(n_vars)]
    data.columns = var_names

    base_name = os.path.splitext(os.path.basename(csv_path))[0]
    output_dir = "causal_DAGs_X_Robust"
    os.makedirs(output_dir, exist_ok=True)

    # CSVと変数対応表を保存
    mapping_df = pd.DataFrame({"NewName": var_names, "OriginalName": original_var_names})
    mapping_df.to_csv(os.path.join(output_dir, f"{base_name}_variable_mapping.csv"), index=False)
    data.to_csv(os.path.join(output_dir, f"{base_name}_X.csv"), index=False)

    # summary用のログ
    summary_lines = []
    summary_lines.append(f"CSVファイル: {csv_path}\n\n")
    summary_lines.append("変数名対応表:\n")
    for new, orig in zip(var_names, original_var_names):
        summary_lines.append(f"{new} -> {orig}\n")
    summary_lines.append("\n")

    # 離散変数チェック
    discrete_cols = detect_discrete_columns(data, threshold=5)
    summary_lines.append(f"離散変数検出: {discrete_cols if discrete_cols else 'なし'}\n")

    # 正規性検定は連続変数に対してのみ実行
    continuous_cols = [col for col in var_names if col not in discrete_cols]
    shapiro_results = []
    non_normal_count = 0

    if continuous_cols:
        summary_lines.append("\nシャピロ・ウィルク検定結果:\n")
        for col in continuous_cols:
            stat, p = shapiro(data[col])
            shapiro_results.append((col, p))
            summary_lines.append(f"{col}: p = {p:.4f}\n")
            if p < 0.05:
                non_normal_count += 1
        summary_lines.append(f"\np < 0.05 の変数数（非正規と判定）: {non_normal_count}\n")
    else:
        summary_lines.append("連続変数が存在しないため、シャピロ・ウィルク検定はスキップされました。\n")

    # 独立性検定の選択
    if discrete_cols:
        indep_test = 'chisq'
        summary_lines.append("独立性検定方法: カイ二乗検定（離散変数を含むため）\n")
    elif non_normal_count > 0:
        indep_test = kci
        summary_lines.append("独立性検定方法: KCIT（ノンパラメトリック）\n")
    else:
        indep_test = 'fisherz'
        summary_lines.append("独立性検定方法: Fisher-z検定（全連続変数が正規性を満たすため）\n")

    # 【論文記載用の説明】を追記
    paper_description_ja = """
本解析では、変数間の因果構造を明らかにするために、制約ベースのPCアルゴリズムを使用して因果有向非巡回グラフ（DAG）を推定した。
連続変数については、まずシャピロ・ウィルク検定により正規性を評価し、非正規分布が疑われる変数が存在する場合は、ノンパラメトリックな独立性検定（Kernel-based Conditional Independence Test: KCIT）を用いた。
離散変数が含まれる場合はカイ二乗検定を、すべての変数が正規分布に従うと判断された場合はFisher-z検定を適用した。
因果グラフは NetworkX および matplotlib により可視化された。
"""

    paper_description_en = """
To infer the causal structure among variables, we applied the constraint-based PC algorithm to estimate a causal Directed Acyclic Graph (DAG).
For continuous variables, normality was assessed using the Shapiro–Wilk test. When one or more variables were identified as non-normally distributed (p < 0.05), a nonparametric independence test (Kernel-based Conditional Independence Test: KCIT) was employed.
If discrete variables were detected, the chi-square test was used; otherwise, the Fisher-z test was applied for normally distributed data.
The resulting causal graphs were visualized using NetworkX and matplotlib.
"""

    summary_lines.append("\n\n=== 論文記載用説明（日本語） ===\n")
    summary_lines.append(paper_description_ja.strip() + "\n")
    summary_lines.append("\n=== Description for paper (English) ===\n")
    summary_lines.append(paper_description_en.strip() + "\n")

    # summaryファイル出力
    summary_path = os.path.join(output_dir, f"{base_name}_analysis_summary.txt")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.writelines(summary_lines)
    print(f"\n解析サマリーを保存しました: {summary_path}")

    # PCアルゴリズム実行
    print("\nPCアルゴリズムを実行中...")
    cg = pc(data.values, alpha=0.05, indep_test_method=indep_test)

    G = cg.G
    nodes = list(G.get_nodes())

    nx_graph = nx.DiGraph()
    nx_graph.add_nodes_from(var_names)

    for i, n1 in enumerate(nodes):
        for j, n2 in enumerate(nodes):
            if i != j and G.get_edge(n1, n2) is not None:
                nx_graph.add_edge(var_names[i], var_names[j])

    # 因果グラフ描画
    import matplotlib as mpl
    mpl.rcParams['font.family'] = 'Helvetica'  # ★ Helvetica に変更
    print("因果グラフを描画し保存します...")
    plt.figure(figsize=(8, 8))
    pos = nx.spring_layout(nx_graph, seed=42)
    nx.draw(nx_graph, pos, with_labels=True, node_color="gold",
            node_size=500, font_size=10, arrows=True, arrowsize=7)
    plt.title("Causal DAG (PC Algorithm)", fontsize=5)
    plt.tight_layout()

    pdf_path = os.path.join(output_dir, f"{base_name}_X_causal_DAG.pdf")
    png_path = os.path.join(output_dir, f"{base_name}_X_causal_DAG.png")
    plt.savefig(pdf_path,dpi=300)
    plt.savefig(png_path,dpi=300)
    print(f"グラフを保存しました: {pdf_path}")
    print(f"グラフを保存しました: {png_path}")
    plt.show()

    # TikZファイル出力
    tikz_path = os.path.join(output_dir, f"{base_name}_causal_DAG.tex")
    save_tikz(nx_graph, var_names, tikz_path)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("使い方: python3.10 causal_analysis_X_with_tikz.py データファイル.csv")
        sys.exit(1)
    main(sys.argv[1])
