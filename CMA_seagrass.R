# ============================================================
# Bacteria–Eukaryota SEM対応: 階層的媒介解析 + Chain mediation + 可視化
# ============================================================

# 必要ライブラリ
library(mediation)
library(dplyr)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# ============================================================
# データ読込
# ============================================================
input <- read.csv("EFA_seagrass.csv")
dir.create("./CMA", showWarnings = FALSE)

# ============================================================
# mediation経路定義（SEM構造に対応）
# ============================================================
mediation_paths <- list(
  list(treat = "B_Hyphomonadaceae", mediator = "E_Corallinophycidae", outcome = "B_Desulfobulbaceae"),
  list(treat = "E_Corallinophycidae", mediator = "B_Desulfobulbaceae", outcome = "Seagrass"),
  list(treat = "B_Hyphomonadaceae", mediator = "B_Desulfobulbaceae", outcome = "Seagrass"),
  list(treat = "B_Hyphomonadaceae", mediator = "E_Corallinophycidae", outcome = "Seagrass")
)

all_results <- list()

# ============================================================
# 各媒介解析（Bootstrap & Quasi-Bayes）+ 結果出力
# ============================================================
sink("./CMA/mediation_detailed_results.txt")

for (path in mediation_paths) {
  cat("============================================================\n")
  cat("Treat:", path$treat, " Mediator:", path$mediator, " Outcome:", path$outcome, "\n")
  cat("============================================================\n")
  
  fit_m <- lm(as.formula(paste(path$mediator, "~", path$treat)), data = input)
  fit_y <- lm(as.formula(paste(path$outcome, "~", path$treat, "+", path$mediator)), data = input)
  
  cat("\n--- LM (Mediator ~ Treat) ---\n")
  print(summary(fit_m))
  cat("\n--- LM (Outcome ~ Treat + Mediator) ---\n")
  print(summary(fit_y))
  
  # Quasi-Bayesian
  med_bayes <- tryCatch(
    mediate(fit_m, fit_y, treat = path$treat, mediator = path$mediator, sims = 1000),
    error = function(e) NULL
  )
  if (!is.null(med_bayes)) {
    cat("\n--- Quasi-Bayesian Mediation Analysis ---\n")
    print(summary(med_bayes))
  } else cat("Quasi-Bayes mediation failed\n")
  
  # Bootstrap
  med_boot <- tryCatch(
    mediate(fit_m, fit_y, treat = path$treat, mediator = path$mediator, boot = TRUE, sims = 1000),
    error = function(e) NULL
  )
  if (!is.null(med_boot)) {
    cat("\n--- Bootstrap Mediation Analysis ---\n")
    print(summary(med_boot))
  } else cat("Bootstrap mediation failed\n")
  
  df <- data.frame(
    Treat = path$treat,
    Mediator = path$mediator,
    Outcome = path$outcome,
    Method = c("Bootstrap", "QuasiBayes"),
    ACME = c(ifelse(is.null(med_boot), NA, med_boot$d0),
             ifelse(is.null(med_bayes), NA, med_bayes$d0)),
    ADE = c(ifelse(is.null(med_boot), NA, med_boot$z0),
            ifelse(is.null(med_bayes), NA, med_bayes$z0)),
    Total = c(ifelse(is.null(med_boot), NA, med_boot$tau.coef),
              ifelse(is.null(med_bayes), NA, med_bayes$tau.coef)),
    Prop.Mediated = c(ifelse(is.null(med_boot), NA, med_boot$n0),
                      ifelse(is.null(med_bayes), NA, med_bayes$n0))
  )
  
  all_results <- append(all_results, list(df))
}
sink()

# ============================================================
# 連鎖媒介 (B_Hyphomonadaceae → E_Corallinophycidae → B_Desulfobulbaceae → Seagrass)
# ============================================================
all_df <- bind_rows(all_results)

get_acme_safe <- function(df, treat, med, method) {
  res <- df %>% filter(Treat == treat, Mediator == med, Method == method)
  if (nrow(res) == 0) return(NA)
  return(res$ACME[1])
}

chain_mediation <- function(df, path, method) {
  acmes <- c()
  for (i in 1:(length(path)-1)) {
    acmes <- c(acmes, get_acme_safe(df, path[i], path[i+1], method))
  }
  prod(acmes, na.rm = TRUE)
}

path <- c("B_Hyphomonadaceae", "E_Corallinophycidae", "B_Desulfobulbaceae", "Seagrass")
boot_chain <- chain_mediation(all_df, path, "Bootstrap")
bayes_chain <- chain_mediation(all_df, path, "QuasiBayes")

chain_summary <- data.frame(
  Path = paste(path, collapse = " → "),
  Chain_Indirect_Effect_Boot = boot_chain,
  Chain_Indirect_Effect_Bayes = bayes_chain
)
write.csv(chain_summary, "./CMA/chain_mediation_summary.csv", row.names = FALSE)

sink("./CMA/chain_mediation_comments.txt")
cat("=== Chain Mediation Interpretation ===\n")
cat("Path:", paste(path, collapse = " → "), "\n\n")
cat("Bootstrap estimate:", boot_chain, "\n")
cat("Quasi-Bayesian estimate:", bayes_chain, "\n\n")
cat("Interpretation:\n")
cat("この値は B_Hyphomonadaceae → E_Corallinophycidae → B_Desulfobulbaceae → Seagrass の間接効果\n")
cat("正なら促進的媒介、負なら抑制的媒介を意味する\n")
sink()

cat("mediation解析・chain mediation\n")

# ============================================================
# SEM + CMA可視化（縦長、p値アスタリスク付）
# ============================================================

edges <- data.frame(
  from = c("B_Hyphomonadaceae", "E_Corallinophycidae", "B_Hyphomonadaceae",
           "B_Desulfobulbaceae", "E_Corallinophycidae", "B_Hyphomonadaceae"),
  to   = c("E_Corallinophycidae", "B_Desulfobulbaceae", "B_Desulfobulbaceae",
           "Seagrass", "Seagrass", "Seagrass"),
  est  = c(11.65, -0.057, -6.68, 0.041, 0.026, 0.857),
  p    = c(0.0001, 0.42, 0.0002, 0.005, 0.008, 0.03)
)

get_p_star <- function(p) {
  sapply(p, function(x) {
    if (is.na(x)) return("")
    if (x < 0.001) return("***")
    else if (x < 0.01) return("**")
    else if (x < 0.05) return("*")
    else return("(ns)")
  })
}

edges$label <- paste0(edges$est, get_p_star(edges$p))
edges$color <- ifelse(edges$p < 0.05 & edges$est > 0, "green",
                 ifelse(edges$p < 0.05 & edges$est < 0, "red", "gray"))
edges$style <- ifelse(edges$p < 0.05, "solid", "dashed")

edge_lines <- paste0(edges$from, " -> ", edges$to,
                     " [label='", edges$label, "', color='",
                     edges$color, "', style='", edges$style, "'];")

dot_code <- paste0("
digraph SEM_CMA_BacteriaEukaryota {
  graph [layout = dot, rankdir = TB]  # 縦方向 (Top→Bottom)
  node [fontname = Helvetica, shape = rectangle, style = filled, fillcolor = 'white'];
  ", paste(edge_lines, collapse = "\n  "), "
  label = 'SEM + CMA model (Bacteria–Eukaryota)\\nGreen: positive sig, Red: negative sig, Gray: ns';
}")

graph <- grViz(dot_code)

svg_code <- export_svg(graph)
rsvg_pdf(charToRaw(svg_code), file = "./CMA/SEM_Bacteria_Eukaryota_vertical.pdf")
rsvg_png(charToRaw(svg_code),
         file = "./CMA/SEM_Bacteria_Eukaryota_vertical.png",
         width = 2000, height = 2800)

cat("SEM+CMA統合グラフ（縦長）を ./CMA へ\n")
