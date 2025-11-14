#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Raw data https://www.data.jma.go.jp/kaiyou/data/db/kaikyo/series/engan/eg_areano.html
#使用データ　area306/319/508/509/518/602
#python3.10 txtデータから時系列グラフ.py *.csv

#pip install pandas matplotlib

#Time course/
#├─ merged_data.csv             ← すべてのファイルを統合したCSV
#├─ yearly_area508_flagP.png
#├─ monthly_area508_flagP.png
#├─ yearly_area509_flagP.png
#├─ monthly_area509_flagP.png

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams

# フォント設定
rcParams['font.family'] = 'Helvetica'
rcParams['font.size'] = 16  # サイズ変更可

def detect_date_column(df):
    df.columns = [c.lower().replace(".", "").replace(" ", "").replace("　","") for c in df.columns]
    if all(x in df.columns for x in ["yyyy","mm","dd"]):
        for col in ["yyyy","mm","dd"]:
            df[col] = pd.to_numeric(df[col].astype(str).str.strip(), errors="coerce")
        df = df.dropna(subset=["yyyy","mm","dd"])
        df["date"] = pd.to_datetime(df.apply(lambda row: f"{int(row['yyyy'])}-{int(row['mm']):02d}-{int(row['dd']):02d}", axis=1))
        return df
    elif all(x in df.columns for x in ["year","month","day"]):
        for col in ["year","month","day"]:
            df[col] = pd.to_numeric(df[col].astype(str).str.strip(), errors="coerce")
        df = df.dropna(subset=["year","month","day"])
        df["date"] = pd.to_datetime(df.apply(lambda row: f"{int(row['year'])}-{int(row['month']):02d}-{int(row['day']):02d}", axis=1))
        return df
    elif "date" in df.columns:
        df["date"] = pd.to_datetime(df["date"].astype(str).str.strip(), errors="coerce")
        return df
    print("❌ 日付列なし")
    sys.exit(1)

def save_both_formats(fig, filepath_base):
    """PNG と PDF 両方で保存"""
    fig.savefig(filepath_base + ".png", dpi=300)
    fig.savefig(filepath_base + ".pdf", dpi=300)
    plt.close(fig)

def main():
    if len(sys.argv) < 2:
        print("使用法: python コード.py file1.csv file2.csv ...")
        sys.exit(1)

    files = sys.argv[1:]
    output_dir = "Time course"
    os.makedirs(output_dir, exist_ok=True)
    all_data = []

    for f in files:
        try:
            df = pd.read_csv(f, sep=",", header=0, encoding="utf-8")
            df = detect_date_column(df)
            df.columns = [c.lower().replace(".", "").replace(" ", "").replace("　","") for c in df.columns]
            all_data.append(df)
            print(f"読込完了: {f}")
        except Exception as e:
            print(f"読込失敗: {f} ({e})")

    if not all_data:
        print("読込データなし")
        sys.exit(1)

    df_all = pd.concat(all_data, ignore_index=True)
    df_all = df_all[df_all["date"].dt.year >= 2016]

    if df_all.empty:
        print("2016年以降のデータなし")
        sys.exit(1)

    df_all["year"] = df_all["date"].dt.year
    df_all["month"] = df_all["date"].dt.month

    df_yearly = df_all.groupby(["areano","year"])["temp"].mean().reset_index()
    df_monthly = df_all.groupby(["areano","year","month"])["temp"].mean().reset_index()
    df_monthly_10yr = df_monthly.groupby(["areano","month"])["temp"].mean().reset_index()  # 10年平均

    df_yearly.to_csv(os.path.join(output_dir, "yearly_average.csv"), index=False)
    df_monthly.to_csv(os.path.join(output_dir, "monthly_average.csv"), index=False)
    df_monthly_10yr.to_csv(os.path.join(output_dir, "monthly_10yr_average.csv"), index=False)

    area_list = sorted(df_all["areano"].unique())

    # --- 個別グラフ ---
    for area in area_list:
        df_area_year = df_yearly[df_yearly["areano"] == area]
        df_area_month = df_monthly[df_monthly["areano"] == area].copy()
        df_area_month["year_month"] = pd.to_datetime(df_area_month[["year","month"]].assign(day=1))

        # 年平均
        fig = plt.figure(figsize=(10,8))　#サイズ変更可
        plt.plot(df_area_year["year"], df_area_year["temp"], marker="o")
        plt.title(f"Area {area} Yearly Average Temp")
        plt.xlabel("Year")
        plt.ylabel("Temperature")
        plt.grid(True)
        plt.tight_layout()
        save_both_formats(fig, os.path.join(output_dir, f"Area{area}_Yearly"))

        # 月平均
        fig = plt.figure(figsize=(10,8))
        plt.plot(df_area_month["year_month"], df_area_month["temp"], marker="o")
        plt.title(f"Area {area} Monthly Average Temp")
        plt.xlabel("Year-Month")
        plt.ylabel("Temperature")
        plt.grid(True)
        plt.tight_layout()
        save_both_formats(fig, os.path.join(output_dir, f"Area{area}_Monthly"))

        # 月平均（年度別折れ線グラフ）
        import numpy as np
        from matplotlib import cm

        df_area_month_copy = df_area_month.copy()
        df_area_month_copy["month_num"] = df_area_month_copy["month"]

        fig = plt.figure(figsize=(10,8))
        years = sorted(df_area_month_copy["year"].unique())
        colors = cm.viridis(np.linspace(0,1,len(years)))  # カラーマップ

        for i, year in enumerate(years):
            df_year = df_area_month_copy[df_area_month_copy["year"] == year]
            plt.plot(df_year["month_num"], df_year["temp"], marker="o", label=str(year), color=colors[i])

        plt.title(f"Area {area} Monthly Temp by Year")
        plt.xlabel("Month")
        plt.ylabel("Temperature")
        plt.xticks(range(1,13))
        plt.grid(True)
        plt.legend(title="Year")
        plt.tight_layout()
        save_both_formats(fig, os.path.join(output_dir, f"Area{area}_Monthly_ByYear"))

    # --- まとめグラフ ---
    # 年平均まとめ
    fig = plt.figure(figsize=(10,8))
    for area in area_list:
        df_area_year = df_yearly[df_yearly["areano"] == area]
        plt.plot(df_area_year["year"], df_area_year["temp"], marker="o", label=f"Area {area}")
    plt.title("Yearly Average Temp (All Areas)")
    plt.xlabel("Year")
    plt.ylabel("Temperature")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    years = range(df_yearly["year"].min(), df_yearly["year"].max()+1)
    plt.xticks(years)
    save_both_formats(fig, os.path.join(output_dir, "AllAreas_Yearly"))

    # 月別10年平均まとめ（折れ線）
    fig = plt.figure(figsize=(10,8))
    for area in area_list:
        df_area_10yr = df_monthly_10yr[df_monthly_10yr["areano"] == area]
        plt.plot(df_area_10yr["month"], df_area_10yr["temp"], marker="o", label=f"Area {area}")
    plt.title("Monthly 10-Year Average Temp (All Areas)")
    plt.xlabel("Month")
    plt.ylabel("Temperature")
    plt.xticks(range(1,13))
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    save_both_formats(fig, os.path.join(output_dir, "AllAreas_Monthly_10yrAvg_Line"))

    # 月別10年平均まとめ（面グラフ）
    fig = plt.figure(figsize=(10,8))
    for area in area_list:
        df_area_10yr = df_monthly_10yr[df_monthly_10yr["areano"] == area]
        plt.fill_between(df_area_10yr["month"], df_area_10yr["temp"], alpha=0.3)
        plt.plot(df_area_10yr["month"], df_area_10yr["temp"], marker="o", label=f"Area {area}")
    plt.title("Monthly 10-Year Average Temp (All Areas, Area Plot)")
    plt.xlabel("Month")
    plt.ylabel("Temperature")
    plt.xticks(range(1,13))
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    save_both_formats(fig, os.path.join(output_dir, "AllAreas_Monthly_10yrAvg_Area"))

    # 月別10年平均まとめ（棒グラフ）
    fig = plt.figure(figsize=(10,8))
    width = 0.8 / len(area_list)
    for i, area in enumerate(area_list):
        df_area_10yr = df_monthly_10yr[df_monthly_10yr["areano"] == area]
        plt.bar(df_area_10yr["month"] + i*width, df_area_10yr["temp"], width=width, label=f"Area {area}")
    plt.title("Monthly 10-Year Average Temp (All Areas, Bar)")
    plt.xlabel("Month")
    plt.ylabel("Temperature")
    plt.xticks(range(1,13))
    plt.grid(True, axis="y")
    plt.legend()
    plt.tight_layout()
    save_both_formats(fig, os.path.join(output_dir, "AllAreas_Monthly_10yrAvg_Bar"))

    # 月別10年平均まとめ（箱ひげ）
    fig = plt.figure(figsize=(10,8))
    for idx, area in enumerate(area_list):
        df_area_10yr = df_monthly_10yr[df_monthly_10yr["areano"] == area]
        plt.boxplot(df_area_10yr["temp"], positions=[idx+1], widths=0.6)
    plt.title("Monthly 10-Year Average Temp (All Areas, Boxplot)")
    plt.xlabel("Area")
    plt.ylabel("Temperature")
    plt.xticks(range(1, len(area_list)+1), area_list)
    plt.grid(True, axis="y")
    plt.tight_layout()
    save_both_formats(fig, os.path.join(output_dir, "AllAreas_Monthly_10yrAvg_Boxplot"))

    print(f"個別グラフとまとめグラフ格納 '{output_dir}' へ PNG & PDF (dpi=300, Helvetica) 出力")

if __name__ == "__main__":
    main()
