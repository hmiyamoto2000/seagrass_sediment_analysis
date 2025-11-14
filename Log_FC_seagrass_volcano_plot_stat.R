#正規化していない生データ
#your_counts.csv の例（形式）
#GeneID    Sample1    Sample2    Sample3    Sample4
#GeneA    523    489    512    530
#GeneB    34    28    30    32

#group.csv
#Sample    Group
#Sample1    Cont
#Sample2    Cont
#>>>>
#Sample4    EPS

# === 必要なパッケージ ===
#if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
#BiocManager::install("EnhancedVolcano", ask=FALSE)
#install.packages(c("vegan","ggplot2","pheatmap","RColorBrewer","dplyr"))

# === パッケージ読み込み ===
library(edgeR)
library(limma)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)
library(vegan)
library(dplyr)

# === データ読み込み ===
count <- read.csv("Pathway_whole.csv", row.names=1)  # 生カウントデータyour_counts.csvは自分のデータで
group_df <- read.csv("Group_list.csv", row.names=1)    # サンプルとグループ情報group.csvデータで。
group <- factor(group_df$Group)

# === limma-voom 差次的発現解析パイプライン ===
y <- DGEList(counts = count, group = group)
keep <- filterByExpr(y, group = group)
y <- y[keep,, keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

design <- model.matrix(~group)
v <- voom(y, design, plot=TRUE)

fit <- lmFit(v, design)
fit <- eBayes(fit)
res <- topTable(fit, coef=2, number=Inf, sort.by="P")

# ファイルに保存
write.csv(res, "DEG_results.csv")


# adj.P.Val（FDR） < 0.05 の差次的発現遺伝子を抽出
sig_genes <- res %>% filter(adj.P.Val < 0.05)

# 結果をCSVに保存（任意）
write.csv(sig_genes, file = "DEG_FDR_0.05.csv", row.names = TRUE)

# いくつ検出されたか表示
cat("FDR < 0.05の遺伝子数:", nrow(sig_genes), "\n")

# === 可視化準備 ===
expr <- v$E  # voomで正規化＋log2変換された発現量データ

# === MDSプロット ===
dist_mat <- dist(t(expr))
mds <- cmdscale(dist_mat, k=2)
mds_df <- as.data.frame(mds)
colnames(mds_df) <- c("MDS1","MDS2")
mds_df$Group <- group

p_mds <- ggplot(mds_df, aes(x=MDS1, y=MDS2, color=Group)) +
  geom_point(size=4) + theme_bw() + labs(title="MDS Plot")
ggsave("MDS_plot.png", plot=p_mds, width=6, height=5, dpi=300)

# === NMDSプロット ===
set.seed(123)
nmds <- metaMDS(t(expr), distance="bray", k=2, trymax=100)
nmds_df <- as.data.frame(scores(nmds)$sites)
nmds_df$Group <- group

p_nmds <- ggplot(nmds_df, aes(x=NMDS1, y=NMDS2, color=Group)) +
  geom_point(size=4) + theme_bw() + labs(title="NMDS Plot")
ggsave("NMDS_plot.png", plot=p_nmds, width=6, height=5, dpi=300)

# === Volcano Plot ===
p_volcano <- EnhancedVolcano(res,
                             lab = rownames(res),
                             x = 'logFC',
                             y = 'P.Value',
                             pCutoff = 0.05,
                             FCcutoff = 1.0,
                             title = 'Volcano Plot',
                             legendPosition = 'bottom')
ggsave("Volcano_plot.png", plot=p_volcano, width=8, height=6, dpi=300)

# === ヒートマップ（FDR<0.1 & P.Value<0.05 上位100遺伝子） ===
sig <- res %>% filter(adj.P.Val < 0.1 & P.Value < 0.05) %>% arrange(desc(abs(logFC)))
topN <- head(sig, 100)
genes <- rownames(topN)

heatmat <- expr[genes, , drop=FALSE]

order_idx <- order(group)
heatmat <- heatmat[, order_idx]
group_sorted <- group[order_idx]

annotation_col <- data.frame(Group = group_sorted)
rownames(annotation_col) <- colnames(heatmat)

pheatmap(heatmat,
         annotation_col = annotation_col,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(255),
         show_rownames = FALSE,
         show_colnames = TRUE,
         fontsize_col = 10,
         filename = "Heatmap_FDR0.1_p0.05_top100.png",
         dpi = 300,
         width = 8,
         height = 10)

# === ログFold ChangeやFDRカウントの出力 ===
cat("FDR < 0.05 count:", sum(res$adj.P.Val < 0.05), "\n")
cat("FDR < 0.1 count:", sum(res$adj.P.Val < 0.1), "\n")

## Volcano PlotのpCutoffをFDR（adj.P.Val）に合わせる場合
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'logFC',
                y = 'adj.P.Val',   # FDRを使う場合はこちら
                pCutoff = 0.05,
                FCcutoff = 1.0,
                title = 'Volcano Plot (FDR < 0.05)',
                legendPosition = 'bottom')

# ヒートマップの例（FDR < 0.05 かつ |logFC| > 1 などの条件も可能）
sig_heatmap <- res %>% filter(adj.P.Val < 0.05 & abs(logFC) > 1) %>% arrange(adj.P.Val)
topN <- head(sig_heatmap, 100)  # 上位100遺伝子抽出
genes <- rownames(topN)
heatmat <- expr[genes, , drop=FALSE]


# Volcano Plotを作成し、オブジェクトに保存
p_volcano <- EnhancedVolcano(res,
                             lab = rownames(res),
                             x = 'logFC',
                             y = 'adj.P.Val',   # FDRを使う場合はこちら
                             pCutoff = 0.05,
                             FCcutoff = 1.0,
                             title = 'Volcano Plot (FDR < 0.05)',
                             legendPosition = 'bottom')

# dpi=300でPNG形式で保存
ggsave(filename = "Volcano_plot_FDR0.05.png", plot = p_volcano,
       width = 8, height = 6, dpi = 300)


#統計比較
# すでにv, y, designがある前提で関数を呼び出すだけ

evaluate_voom_stats <- function(v, y, design, alpha = 0.05, output_prefix = "voom_evaluation") {
  # voom.xyが存在しない場合は自作計算
  if (is.null(v$voom.xy)) {
    mean_expr <- rowMeans(v$E)
    var_expr <- apply(v$E, 1, var)
    x <- mean_expr
    y_voom <- log(var_expr + 1e-8)  # 分散にlogと微小値を足して安定化
  } else {
    x <- v$voom.xy$x
    y_voom <- v$voom.xy$y
  }
  
  # Loessフィッティングと残差計算
  loess_fit <- loess(y_voom ~ x)
  loess_pred <- predict(loess_fit)
  residuals <- y_voom - loess_pred
  
  # R2計算
  SSE <- sum(residuals^2)
  TSS <- sum((y_voom - mean(y_voom))^2)
  R2 <- 1 - SSE/TSS
  
  # RMSE (残差の標準誤差)
  RMSE <- sqrt(mean(residuals^2))
  
  # IQR of voom weights
  weights <- as.vector(v$weights)
  IQR_weights <- IQR(weights)
  
  # Spearman相関
  spearman_cor <- cor(x, y_voom, method = "spearman")
  
  # DEG数(voom)
  fit_voom <- lmFit(v, design)
  fit_voom <- eBayes(fit_voom)
  deg_voom <- topTable(fit_voom, number=Inf, adjust.method="BH")
  n_deg_voom <- sum(deg_voom$adj.P.Val < alpha)
  
  # DEG数(raw logCPM)
  logCPM <- cpm(y, log=TRUE, prior.count=2)
  fit_raw <- lmFit(logCPM, design)
  fit_raw <- eBayes(fit_raw)
  deg_raw <- topTable(fit_raw, number=Inf, adjust.method="BH")
  n_deg_raw <- sum(deg_raw$adj.P.Val < alpha)
  
  # 結果まとめ
  result_df <- data.frame(
    Metric = c("Loess_R2", "RMSE_smoothness", "IQR_voom_weights", "Spearman_correlation",
               paste0("DEGs_voom_adjP_", alpha), paste0("DEGs_raw_adjP_", alpha)),
    Value = c(round(R2, 4), round(RMSE, 4), round(IQR_weights, 4), round(spearman_cor, 4),
              n_deg_voom, n_deg_raw)
  )
  
  # ファイル名
  csv_file <- paste0(output_prefix, "_stats.csv")
  txt_file <- paste0(output_prefix, "_stats.txt")
  
  # 保存
  write.csv(result_df, file = csv_file, row.names = FALSE)
  write.table(result_df, file = txt_file, row.names = FALSE, quote = FALSE, sep = "\t")
  
  message("✅ voom統計評価を保存しました: ", csv_file, ", ", txt_file)
  return(result_df)
}

# 呼び出し例（すでにv, y, designが定義されていれば）
stats_res <- evaluate_voom_stats(v, y, design, alpha = 0.05, output_prefix = "voom_evaluation")
print(stats_res)

# stats_res データフレームが既にある前提で

# CSVで保存
write.csv(stats_res, file = "voom_evaluation_results.csv", row.names = FALSE)

# タブ区切りテキストで保存
write.table(stats_res, file = "voom_evaluation_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# 保存完了メッセージ
cat("✅ voom評価結果を以下のファイルに保存しました：\n",
    "CSVファイル: voom_evaluation_results.csv\n",
    "テキストファイル: voom_evaluation_results.txt\n")



# Metric	Value	意味と妥当性
# Loess_R2 (0.9253)	0.9超えで良好	loessでの平均と分散の関係モデルへの適合度。高いR2はトレンドがしっかり捉えられている証拠。
# RMSE_smoothness (0.6616)	小さい値は良好	平滑化後の残差のばらつき。0.66は適度に滑らかなトレンドが得られていることを示唆。
# IQR_voom_weights (47.6624)	数値の大きさはデータ依存	voomの重みの分布の広がり。サンプル数や発現範囲によって変動。大きすぎず小さすぎず。
# Spearman_correlation (-0.9538)	絶対値が1に近い	順位相関が強い（負の相関で期待通りの減少傾向）。voomの平均-分散トレンドを示す指標として理想的。
# DEGs_voom_adjP_0.05 (260)	DEGs数	voom解析でFDR<0.05の遺伝子数。生物学的に意味があれば妥当。
# DEGs_raw_adjP_0.05 (37)	DEGs数	正規化なし生データでのDEG数。voomと比べて少ないのはvoomのパワー向上を示唆。
#

##https://note.com/ash_coding/n/n6ef8afc4494c
library(ggplot2)
library(tidyverse)

# データフレーム res を使用（topTableの出力）
# 必要なら以下で変換しておく：
res <- res %>%
  mutate(log10FDR = -log10(adj.P.Val),
         Significant = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, "Yes", "No"))

# Volcano Plot作成
p_volcano_custom <- ggplot(res, aes(x = logFC, y = log10FDR)) +
  geom_point(aes(color = Significant), size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red", size = 0.4) +
  scale_color_manual(values = c("Yes" = "firebrick", "No" = "grey")) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "Volcano Plot (FDR < 0.05 & |logFC| > 1)",
       x = "log2 Fold Change",
       y = "-log10(FDR)")

# 表示
print(p_volcano_custom)

# 保存（PNG, 300dpi）
ggsave("Volcano_plot_custom_FDR.png", plot = p_volcano_custom,
       width = 8, height = 6, dpi = 300)



#グラフその2
library(ggplot2)
library(tidyverse)

# FDRや有意差の列を追加
res <- res %>%
  mutate(log10FDR = -log10(adj.P.Val),
         Significant = ifelse(adj.P.Val < 0.05 & abs(logFC) > 1, "Yes", "No"))

# プロット作成（フォントは theme() で指定）
p_volcano_custom <- ggplot(res, aes(x = logFC, y = log10FDR)) +
  geom_point(aes(color = Significant), size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red", size = 0.4) +
  scale_color_manual(values = c("Yes" = "firebrick", "No" = "grey")) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "Helvetica")  # フォント指定はここ！
  ) +
  labs(
    title = expression("Volcano Plot (FDR < 0.05 & " * abs(log[2]*FC) * " > 1)"),
    x = expression(log[2]*" Fold Change"),
    y = expression(-log[10]*"(FDR)")
  )

# 表示
print(p_volcano_custom)

# === 保存設定 ===

# PNG保存（300dpi）
ggsave("Volcano_plot_FDR0.05_logFC1.png",
       plot = p_volcano_custom,
       width = 8, height = 6, dpi = 300)

# PDF保存（ベクター）
ggsave("Volcano_plot_FDR0.05_logFC1.pdf",
       plot = p_volcano_custom,
       width = 8, height = 6)


#有意なデータの排出
# 有意な遺伝子（赤で表示された点）だけ抽出
sig_res <- res %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

# 行名（遺伝子ID）を列として追加
sig_res_out <- sig_res %>%
  rownames_to_column(var = "GeneID")  # tibble::rownames_to_column()が必要

# CSVで保存
write.csv(sig_res_out, file = "Significant_genes_volcano.csv", row.names = FALSE)

# txt（タブ区切り）で保存
write.table(sig_res_out, file = "Significant_genes_volcano.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("✅ 有意な遺伝子を保存しました（FDR < 0.05 & |logFC| > 1）\n")


#heatmap選抜　FDR<0.05 LogFC絶対値>1
library(pheatmap)
library(RColorBrewer)

# 前処理は同じ
sig_genes <- res %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1) %>%
  arrange(adj.P.Val)

topN <- head(sig_genes, 100)
genes <- rownames(topN)

heatmat <- expr[genes, , drop = FALSE]

order_idx <- order(group)
heatmat <- heatmat[, order_idx]
group_sorted <- group[order_idx]

annotation_col <- data.frame(Group = group_sorted)
rownames(annotation_col) <- colnames(heatmat)

# --- PDF保存 ---
pdf("Heatmap_FDR0.05_logFC1_top100.pdf", width = 8, height = 10)
pheatmap(heatmat,
         annotation_col = annotation_col,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(255),
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_col = 10)
dev.off()

# --- PNG保存 ---
png("Heatmap_FDR0.05_logFC1_top100.png",
    width = 8*300, height = 10*300, res = 300)  # 8インチ×300dpiなど指定
pheatmap(heatmat,
         annotation_col = annotation_col,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(255),
         show_rownames = FALSE,
         show_colnames = TRUE,
         fontsize_col = 10)
dev.off()

#group分　地域わけ
# group_df には Sample名（行名）、Group（seagrass/non）、env（地域）が含まれている前提

group <- factor(group_df$Group)
group_env <- factor(group_df$env)

# annotation_col を複数の列（地域とグループ）で作成
annotation_col <- data.frame(
  Group = group,
  Environment = group_env
)
rownames(annotation_col) <- rownames(group_df)  # もしくは colnames(heatmat)

# ヒートマップの列の順番をgroupで並び替えたい場合
order_idx <- order(group, group_env)  # 例えばGroupでまず並べ、次に地域で並べる
heatmat <- heatmat[, order_idx]
annotation_col <- annotation_col[order_idx, ]

# --- PDF保存 ---
pdf("Heatmap_FDR0.05_logFC1_top100_region.pdf", width = 8, height = 10)
pheatmap(heatmat,
         annotation_col = annotation_col,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(255),
         show_rownames = TRUE,  # 遺伝子名を表示したいならTRUE
         show_colnames = TRUE,
         fontsize_col = 10)
dev.off()

# pheatmapで表示 png
png("Heatmap_FDR0.05_logFC1_top100_region.png",
    width = 8*300, height = 10*300, res = 300)  # 8インチ×300dpiなど指定

pheatmap(heatmat,
         annotation_col = annotation_col,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(255),
         show_rownames = TRUE,  # 遺伝子名を表示したいならTRUE
         show_colnames = TRUE,
         fontsize_col = 10)

dev.off()


#中央値まとめ
library(pheatmap)
library(RColorBrewer)
library(dplyr)

# データ読み込み例（適宜パスを変更してください）
# expr <- read.csv("expression_data.csv", row.names = 1)  # 遺伝子発現マトリックス
group_df <- read.csv("Group_list.csv", row.names = 1)    # サンプルごとのグループ情報

# 1. サンプル名（列名/行名）の確認と統一（大文字小文字など）
colnames(expr) <- toupper(colnames(expr))
rownames(group_df) <- toupper(rownames(group_df))

# 2. 共通するサンプルだけ抽出し順番を揃える
common_samples <- intersect(colnames(expr), rownames(group_df))

expr <- expr[, common_samples]
group_df <- group_df[common_samples, ]

# 3. 差分チェック（必要に応じて確認）
cat("exprにだけあるサンプル:", setdiff(colnames(expr), rownames(group_df)), "\n")
cat("group_dfにだけあるサンプル:", setdiff(rownames(group_df), colnames(expr)), "\n")

# 4. 差分がなければグループ情報を因子化
group <- factor(group_df$Group)
group_env <- factor(group_df$env)

# 5. DE解析結果が入っているデータフレーム res がある前提でフィルタリング
sig_genes <- res %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1) %>%
  arrange(adj.P.Val)

topN <- head(sig_genes, 100)
genes <- rownames(topN)

# 6. 発現マトリックスを遺伝子で抽出
heatmat <- expr[genes, , drop = FALSE]

# 7. サンプルを環境（env）とグループ（Group）で並び替え
order_idx <- order(group_env, group)
heatmat <- heatmat[, order_idx]
group_sorted <- group[order_idx]
env_sorted <- group_env[order_idx]

# 8. 注釈データフレーム作成
annotation_col <- data.frame(
  Environment = env_sorted,
  Group = group_sorted
)
rownames(annotation_col) <- colnames(heatmat)

# 8.5 因子レベルをカラー指定に合わせる
annotation_col$Environment <- factor(annotation_col$Environment, levels = c("Seagrass", "Non"))
annotation_col$Group <- factor(annotation_col$Group)  # 必要に応じて

# 9. カラーパレット設定
group_levels <- levels(annotation_col$Group)
n_group <- length(group_levels)

if(n_group >= 3){
  group_colors <- brewer.pal(n_group, "Set3")
} else if(n_group == 2){
  group_colors <- c("#66c2a5", "#fc8d62")  # 適宜好きな色を選んでください
} else if(n_group == 1){
  group_colors <- c("#8da0cb")
} else {
  stop("Groupのレベルが空です。")
}

# 10. PDF保存
pdf("Heatmap_FDR0.05_logFC1_top100_region_median.pdf", width = 8, height = 10)
pheatmap(heatmat,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(255),
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_col = 10)
dev.off()

# 11. PNG保存
png("Heatmap_FDR0.05_logFC1_top100_median.png",
    width = 8 * 300, height = 10 * 300, res = 300)
pheatmap(heatmat,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(255),
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_col = 10)
dev.off()

#八つの地域
library(RColorBrewer)
library(scales)  # 色の調整に使う

# 地域のレベル
region_levels <- levels(annotation_col$Region)
n_region <- length(region_levels)  # 8

# 基本色を8色用意
base_colors <- brewer.pal(n_region, "Set3")

# 濃い色: seagrass, 薄い色: non
# 濃淡調整関数
darker <- function(col) { col }  # そのまま（濃い色）
lighter <- function(col) { alpha(col, 0.5) }  # 透明度で薄くする例

# Regionごとにseagrass/nonの色を作る
region_env_colors <- c()
for (i in seq_along(region_levels)) {
  region <- region_levels[i]
  seagrass_col <- darker(base_colors[i])
  non_col <- lighter(base_colors[i])
  
  region_env_colors[paste0(region, "_seagrass")] <- seagrass_col
  region_env_colors[paste0(region, "_non")] <- non_col
}

# annotation_colの新列「Region_Env」を作成
annotation_col$Region_Env <- paste0(annotation_col$Region, "_", annotation_col$Environment)

# annotation_colorsはRegion_Envだけ使う
ann_colors <- list(
  Region_Env = region_env_colors
)

# heatmap描画時の設定例
library(pheatmap)
pheatmap(heatmat,
         annotation_col = data.frame(Region_Env=annotation_col$Region_Env),
         annotation_colors = ann_colors,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255),
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_col = 10)

# ヒートマップ描画

# 11. PNG保存
png("Heatmap_FDR0.05_logFC1_top100_median_each2.png",
    width = 8 * 300, height = 10 * 300, res = 300)
pheatmap(heatmat,
         annotation_col = data.frame(Region_Env=annotation_col$Region_Env),
         annotation_colors = ann_colors,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255),
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_col = 10)


dev.off()

#pdf保存
pdf("Heatmap_FDR0.05_logFC1_top100_region_median_each2.pdf", width = 8, height = 10)

pheatmap(heatmat,
         annotation_col = data.frame(Region_Env=annotation_col$Region_Env),
         annotation_colors = ann_colors,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255),
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_col = 10)

dev.off()

