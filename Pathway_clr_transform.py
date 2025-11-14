#FileA.csv
#sample    K00001    K00002    K00003
#S1    10    30    60
#S2    0    20    80
#S3    5    5    90

#python3.10 clr_transform.py FileA.csv

import pandas as pd
import numpy as np
import sys
import os

def clr_transform(df, pseudocount=1e-6):
    """
    CLR変換（Centered Log-Ratio）を行う
    df: サンプル×機能のDataFrame（数値）
    """
    # pseudocount を加えて log(0) を回避
    df = df + pseudocount

    # 幾何平均を各行（サンプル）ごとに計算
    geometric_means = df.apply(lambda row: np.exp(np.mean(np.log(row))), axis=1)

    # 幾何平均で割って log を取る
    clr_df = df.div(geometric_means, axis=0).apply(np.log)

    return clr_df

# === コマンドライン引数処理 ===

if len(sys.argv) != 2:
    print("❌ 使用方法: python3.10 clr_transform.py input_file.csv")
    sys.exit(1)

input_file = sys.argv[1]

if not os.path.isfile(input_file):
    print(f"❌ エラー: ファイル '{input_file}' が見つかりません。")
    sys.exit(1)

# 出力ファイル名の作成
base, ext = os.path.splitext(input_file)
output_file = f"{base}_clr.csv"

# === ファイル読み込み ===

try:
    df_raw = pd.read_csv(input_file)
except Exception as e:
    print(f"❌ CSV読み込み中にエラー: {e}")
    sys.exit(1)

# 最初の列はサンプル名として使う
if 'sample' not in df_raw.columns[0].lower():
    print("❌ エラー: 最初の列にサンプル名（例: sample）が含まれている必要があります。")
    sys.exit(1)

df_raw = df_raw.set_index(df_raw.columns[0])

# 数値列のみ抽出
try:
    df_numeric = df_raw.astype(float)
except:
    print("❌ エラー: 数値以外のデータが含まれています。全ての機能値は数値である必要があります。")
    sys.exit(1)

# === CLR変換 ===
clr_df = clr_transform(df_numeric)

# インデックスを列に戻して出力
clr_df.insert(0, 'sample', clr_df.index)
clr_df.to_csv(output_file, index=False)

print(f"✅ CLR変換完了: 出力ファイル → {output_file}")
