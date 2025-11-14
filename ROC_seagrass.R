#ROC
#https://rpubs.com/kaz_yos/1713 #二つのROCカーブの比較　DeLong's test
#https://note.com/maatan_223220/n/n54ae12ae28cc
#https://qiita.com/insilicomab/items/57902d8838acd795f298
#https://qiita.com/sz_dr/items/96e9306979cb1832d120

#ヒストグラム
#https://hira-labo.com/archives/112
#https://data-science.gr.jp/implementation/ida_r_histogram.html
#http://tips-r.blogspot.com/2017/02/rcol.html
#https://stats.biopapyrus.jp/r/graph/hist.html


# パッケージ
library(Epi)

# 出力フォルダ作成（必要なら）
dir.create("./ROC_seagrass", showWarnings = TRUE, recursive = TRUE)

# データ読み込み
Data <- read.csv("causal_learning_data.csv")

# 菌種名の取得（1列目 group を除く）
bacteria_names <- colnames(Data)[-1]

# AUC 結果を保存するデータフレームを用意
auc_results <- data.frame(
  Feature = character(),
  AUC = numeric(),
  stringsAsFactors = FALSE
)

# ROC解析ループ
for (bacteria in bacteria_names) {
  
  # ROC解析（plot=FALSE で描画しない）
  roc_result <- ROC(test = Data[[bacteria]], stat = Data$Seagrass, plot = FALSE)
  
  # AUC値を取得
  auc_value <- roc_result$AUC
  
  # AUC情報を結果に追加
  auc_results <- rbind(auc_results, data.frame(Feature = bacteria, AUC = auc_value))
  
  # PNG保存
  png_filename <- paste0("./ROC_seagrass/ROC_graph_", bacteria, ".png")
  png(filename = png_filename, width = 6 * 300, height = 6 * 300, res = 300)
  ROC(test = Data[[bacteria]], stat = Data$group, plot = "ROC")
  dev.off()
  
  # PDF保存
  pdf_filename <- paste0("./ROC_seagrass/ROC_graph_", bacteria, ".pdf")
  pdf(file = pdf_filename, width = 6, height = 6)
  ROC(test = Data[[bacteria]], stat = Data$group, plot = "ROC")
  dev.off()
  
  # テキスト保存
  txt_filename <- paste0("./ROC_seagrass/ROC_graph_raw_data_", bacteria, ".txt")
  sink(txt_filename, append = TRUE)
  cat("ROC result for:", bacteria, "\n\n")
  print(roc_result)
  sink()
}

# AUCフィルタ（AUC >= 0.6 または AUC <= 0.4）
filtered_auc <- subset(auc_results, AUC >= 0.6 | AUC <= 0.4)

# 結果をCSVファイルに保存
write.csv(filtered_auc, "./ROC_seagrass/Filtered_AUC_Features.csv", row.names = FALSE)

# whole AUC
whole_auc <- subset(auc_results)

# 結果をCSVファイルに保存
write.csv(whole_auc, "./ROC_seagrass/whole_AUC_Features.csv", row.names = FALSE)


