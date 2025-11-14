import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

# 強調表示用の色指定（赤系以外も取り入れたバランス版）
highlight_colors = {
    'B_D4_Desulfobulbaceae': 'red',         # 赤系濃いめ'#d73027',
    'E_D4_Diatomea': 'red',                  # オレンジ寄り赤'#fc8d59'
    'B_D4_Marinilabiliaceae': 'yellow',         # 明るいオレンジイエロー'#fdae61'
    'B_D4_Desulfobacteraceae': 'darkgreen',        # 緑系明るめ#a6d854
    'E_D3_Fungi': 'darkgreen',        # 緑系明るめ#a6d854
    'B_D4_Cyclobacteriaceae': '#66c2a5',         # 青緑系
    'B_D4_Sva1033': '#1f78b4',                    # 青系濃いめ
    'B_D4_Rhodobacteracae': 'darkbrown',
    'B_D4_Sulfurovaceae': '#756bb1',              # 紫寄り青
    'B_D4_Halieaceae': 'orange',              # 紫寄り青
    'E_D2_Alveolata_other': 'orange',              # 紫寄り青
    'E_D4_Pirsonia': 'darkgray',              # # 1%以上の時
    'E_D4_Nucleariidae': 'darkgreen',              # # 1%以上の時
    'E_D4_Raphidophyceae': 'orange',      # #1%以上の時
    'E_D4_Colpodellida': 'pink',                  # 1%以上の時
    'B_D4_Sandaracinaceae': 'darkgray',              # 明るいオレンジイエロー'#fdae61'
    'B_D4_Cyclobacteriaceae': 'lightblue',              # 1%以上の時
    'E_D4_Sorodiplophrys': 'darkgray',              # 明るいオレンジイエロー'#fdae61'
    'E_D4_Colpodellida': 'pink',                  # 1%以上の時
    'B_D4_Syntrophobacteraceae': '#e7298a',       # ピンク系
    'E_D3_Cercozoa_Novel Clade Gran-1': 'yellow',
    'Others': 'lightgray',
    'B_D3_Gammaproteobacteria.Incertae.Sedis.D4_Unknown.Family': '#fdae6b'  # オレンジ系
}

# 26色カスタムカラーパレット（赤・緑・青バランスよく）
custom_color_list = [
    '#d73027', '#fc8d59', '#fee08b', '#d9ef8b', '#a6d854', '#66c2a5', '#3288bd',
    '#5e4fa2', '#e7298a', '#e6f598', '#a1d76a', '#1b9e77', '#7570b3', '#e78ac3',
    '#b2abd2', '#fdae61', '#f46d43', '#a6761d', '#666666', '#99d594', '#ffffbf',
    '#fdae61', '#d53e4f', '#4d9221', '#276419', '#9e0142'
]

# --- 色割り当て関数（highlight_colors優先） ---
def assign_colors(taxa, others_gray=False):
    colors = []
    available_colors = custom_color_list.copy()

    for taxon in taxa:
        if taxon in highlight_colors:
            colors.append(highlight_colors[taxon])
        elif others_gray and taxon.lower() == 'others':
            colors.append('lightgray')
        else:
            if available_colors:
                colors.append(available_colors.pop(0))
            else:
                colors.append("gray")  # 予備色
    return colors

# --- 地域別2群の100%積み上げ棒グラフ ---
def plot_100stacked_bars(df, output_folder, base_filename, fig_width=10, fig_height=5, dpi=300, font_size=8, others_gray=False):
    num_plots = int(df.shape[1] / 2)
    taxa = df.index.tolist()

    plt.rcParams['font.family'] = ['Helvetica', 'Arial', 'sans-serif']
    plt.rcParams['font.size'] = font_size

    fig, axes = plt.subplots(1, num_plots, figsize=(fig_width, fig_height), sharey=True)
    plt.subplots_adjust(wspace=0.05)

    if num_plots == 1:
        axes = [axes]

    colors = assign_colors(taxa, others_gray=others_gray)

    bar_width = 0.35
    x_positions = [0, 0.4]

    for i in range(num_plots):
        col1 = df.columns[2*i]
        col2 = df.columns[2*i + 1]

        subset = df[[col1, col2]]
        if subset.shape[1] != 2:
            raise ValueError(f"期待される2列ではなく、{subset.shape[1]}列あります: {col1}, {col2}")

        subset_percent = subset.divide(subset.sum(axis=0), axis=1) * 100

        bottom = [0, 0]
        for idx, taxon in enumerate(taxa):
            values = subset_percent.loc[taxon]
            axes[i].bar(x_positions, values, bottom=bottom, label=taxon if i == 0 else "", color=colors[idx], width=bar_width)
            bottom = [bottom[j] + values.iloc[j] for j in range(2)]

        axes[i].set_xticks([pos + bar_width / 2 for pos in x_positions])
        region_name = col1.split('_')[0]
        axes[i].set_xticklabels([f"{region_name} (seagrass)", f"{region_name} (non)"], rotation=45, ha='right', fontsize=font_size)
        axes[i].set_title(region_name, fontsize=font_size+2)
        axes[i].set_ylim(0, 100)
        axes[i].tick_params(axis='y', labelsize=font_size)
        axes[i].tick_params(axis='x', labelsize=font_size)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=font_size, title='Taxa')

    plt.tight_layout()

    os.makedirs(output_folder, exist_ok=True)
    suffix = 'others_gray' if others_gray else 'normal'
    pdf_path = os.path.join(output_folder, f'{base_filename}_stacked_bar_plot_custom_{suffix}.pdf')
    png_path = os.path.join(output_folder, f'{base_filename}_stacked_bar_plot_custom_{suffix}.png')

    plt.savefig(pdf_path, dpi=dpi)
    plt.savefig(png_path, dpi=dpi)
    plt.show()

    print(f"保存完了: {pdf_path}, {png_path}")

# --- 全体まとめた100%積み上げ棒グラフ ---
def plot_combined_100stacked_bar(df, output_folder, base_filename, fig_width=12, fig_height=6, dpi=300, font_size=8):
    plt.rcParams['font.family'] = ['Helvetica', 'Arial', 'sans-serif']
    plt.rcParams['font.size'] = font_size

    taxa = df.index.tolist()
    num_bars = df.shape[1]
    colors = assign_colors(taxa, others_gray=False)

    df_percent = df.divide(df.sum(axis=0), axis=1) * 100

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    bottom = [0] * num_bars
    for idx, taxon in enumerate(taxa):
        values = df_percent.loc[taxon]
        ax.bar(range(num_bars), values, bottom=bottom, label=taxon, color=colors[idx])
        bottom = [bottom[i] + values.iloc[i] for i in range(num_bars)]

    xtick_labels = []
    for col in df.columns:
        if '_' in col:
            region, group = col.split('_', 1)
            xtick_labels.append(f"{region} ({group})")
        else:
            xtick_labels.append(col)

    ax.set_xticks(range(num_bars))
    ax.set_xticklabels(xtick_labels, rotation=45, ha='right', fontsize=font_size)
    ax.set_ylim(0, 100)
    ax.set_ylabel("Relative abundance (%)", fontsize=font_size)
    ax.tick_params(axis='y', labelsize=font_size)
    ax.set_title("100% Stacked Bar Chart by Region and Group", fontsize=font_size + 2)

    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=font_size - 1, title='Taxa')

    plt.tight_layout()

    os.makedirs(output_folder, exist_ok=True)
    pdf_path = os.path.join(output_folder, f'{base_filename}_combined_stacked_bar_custom.pdf')
    png_path = os.path.join(output_folder, f'{base_filename}_combined_stacked_bar_custom.png')

    plt.savefig(pdf_path, dpi=dpi)
    plt.savefig(png_path, dpi=dpi)
    plt.show()

    print(f"保存完了: {pdf_path}, {png_path}")

# --- メイン処理 ---
def main():
    if len(sys.argv) != 2:
        print("使い方: python コード.py データファイル.csv")
        sys.exit(1)

    filename = sys.argv[1]
    base_filename = os.path.splitext(os.path.basename(filename))[0]
    output_folder = "folder_100graph"

    try:
        # マルチインデックスヘッダー対応で読み込み
        df = pd.read_csv(filename, header=[0,1], index_col=0)

        # 列名を「地域_群名」形式にまとめる
        df.columns = [f"{col[0]}_{col[1]}" for col in df.columns]

        # 数値変換（エラーはNaN）
        df = df.apply(pd.to_numeric, errors='coerce')

        # 地域ごとの2群グラフ（通常色）
        plot_100stacked_bars(df, output_folder, base_filename, font_size=8, others_gray=False)
        
        # 地域ごとの2群グラフ（通常色）半分
        plot_100stacked_bars(df, output_folder, base_filename, fig_width=4, font_size=8, others_gray=False)
        
        # Othersを灰色にしたバージョン
        plot_100stacked_bars(df, output_folder, base_filename, font_size=8, others_gray=True)
        
        # Othersを灰色にしたバージョン 半分
        plot_100stacked_bars(df, output_folder, base_filename, fig_width=4, font_size=8, others_gray=True)
        
        # 全体まとめたグラフ
        plot_combined_100stacked_bar(df, output_folder, base_filename, font_size=9)
        
        # 全体まとめたグラフ　半分
        plot_combined_100stacked_bar(df, output_folder, base_filename, fig_width=7, font_size=9)

    except Exception as e:
        print(f"エラーが発生しました: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()
