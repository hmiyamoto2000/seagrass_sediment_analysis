# ============================================================
# 🧩 階層的媒介解析 + Chain mediation + SEM統合可視化（p値付き）
# ============================================================

# 必要パッケージ
library(mediation)
library(dplyr)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# ============================================================
# データ読み込み
# ============================================================
dir.create("./CMA", showWarnings = FALSE)
input <- read.csv("Seagrass_path_selected.csv")

# ============================================================
# mediation 経路定義（SEM構造に対応）
# ============================================================
mediation_paths <- list(
  list(treat = "Seagrass", mediator = "PWY_2941", outcome = "PWY_4722"),
  list(treat = "PWY_6471", mediator = "P241_PWY", outcome = "PWY_4722"),
  list(treat = "PWY_6471", mediator = "P241_PWY", outcome = "PWY_2941")
)

all_results <- list()

# ============================================================
# mediation解析 + 結果出力
# ============================================================
sink("./CMA/mediation_detailed_results.txt")

for (path in mediation_paths) {
  cat("============================================================\n")
  cat("Treat:", path$treat, " Mediator:", path$mediator, " Outcome:", path$outcome, "\n")
  cat("============================================================\n")
  
  fit_m <- lm(as.formula(paste(path$mediator, "~", path$treat)), data = input)
  fit_y <- lm(as.formula(paste(path$outcome, "~", path$treat, "+", path$mediator)), data = input)
  
  cat("\n--- LM: Mediator ~ Treat ---\n")
  print(summary(fit_m))
  
  cat("\n--- LM: Outcome ~ Treat + Mediator ---\n")
  print(summary(fit_y))
  
  med_bayes <- tryCatch(
    mediate(fit_m, fit_y, treat = path$treat, mediator = path$mediator, sims = 1000),
    error = function(e) NULL
  )
  if (!is.null(med_bayes)) {
    cat("\n--- Quasi-Bayesian Mediation Analysis ---\n")
    print(summary(med_bayes))
  } else cat("⚠️ Quasi-Bayes mediation failed\n")
  
  med_boot <- tryCatch(
    mediate(fit_m, fit_y, treat = path$treat, mediator = path$mediator, boot = TRUE, sims = 1000),
    error = function(e) NULL
  )
  if (!is.null(med_boot)) {
    cat("\n--- Bootstrap Mediation Analysis ---\n")
    print(summary(med_boot))
  } else cat("⚠️ Bootstrap mediation failed\n")
  
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
# 連鎖媒介（chain mediation）
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

# 例: Seagrass → PWY_2941 → PWY_4722
path <- c("Seagrass", "PWY_2941", "PWY_4722")

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
cat("この値は Seagrass → PWY_2941 → PWY_4722 の間接効果を表します。\n")
cat("正なら正の媒介効果、負なら負の媒介効果を示します。\n")
sink()

cat("✅ Mediation 解析完了\n")

# ============================================================
# グラフ出力 (縦長・p値に応じたアスタリスク)
# ============================================================

# 効果値とp値を仮に指定（実際はlm/mediate結果から動的に抽出してもOK）
edges <- data.frame(
  from = c("Seagrass", "Seagrass", "PWY_2941",
           "PWY_6471", "PWY_6471", "P241_PWY"),
  to =   c("PWY_2941", "PWY_4722", "PWY_4722",
           "PWY_4722", "P241_PWY", "PWY_2941"),
  est = c(0.604, 2.752, 0.330, 0.093, -0.100, -0.098),
  p = c(0.003, 0.012, 0.245, 0.310, 0.087, 0.230)
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

# DOTコード生成（縦長）
edge_lines <- paste0(edges$from, " -> ", edges$to,
                     " [label = '", edges$label, "', color='",
                     edges$color, "', style='", edges$style, "'];")

dot_code <- paste0("
digraph SEM_CMA_PWY {
  graph [layout = dot, rankdir = TB]  # ← 縦長（TB=Top-Bottom）
  node [fontname = Helvetica, shape = rectangle, style = filled, fillcolor = 'white'];
  ", paste(edge_lines, collapse = "\n  "), "
  label = 'SEM + CMA integrated model\\n(Green=positive sig, Red=negative sig, Gray=ns)';
}")

graph <- grViz(dot_code)

svg_code <- export_svg(graph)
rsvg_pdf(charToRaw(svg_code), file = "./CMA/SEM_CMA_PWY_vertical.pdf")
rsvg_png(charToRaw(svg_code),
         file = "./CMA/SEM_CMA_PWY_vertical.png",
         width = 2000, height = 3000)

cat("✅ 縦長SEM+CMAグラフを ./CMA/ に出力しました。\n")
