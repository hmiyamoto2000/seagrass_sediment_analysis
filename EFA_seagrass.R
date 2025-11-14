#For EFA
dir.create("./EFA_result", showWarnings = TRUE, recursive = FALSE, mode = "0777")

library( psych )
library( GPArotation )

data=read.csv("EFA_seagrass.csv")

KMO(data) #E_Rhodymeniophycidae 0.27であるため削除

k_d=KMO(data)
sink('./EFA_result/KMO_seagrass.txt', append = TRUE)
print (k_d)
sink()

par("mar"=c(5,5,5,5))
fa.parallel(data, fa = "fa", use = "complete.obs")
abline(h = 0)
parallel=fa.parallel(data, fa = "fa", use = "complete.obs")
print(parallel)
abline(h = 0)
sink('./EFA_result/parallel.txt', append = TRUE)
print (parallel)
sink()

MAPminres <- vss(data, fm="minres")
print(MAPminres)

sink('./EFA_result/VSS_MAPminres.txt', append = TRUE)
print (MAPminres)
sink()

#Calculate as factor=8 based on the data using the function "vss"
#Use rotate = "promax" based on the result using CA

result_4 = fa(data, nfactors = 4, fm = "minres", rotate = "promax", use = "complete.obs" )
print( result_4, digits = 4, sort = T)

sink('./EFA_result/nfactors_4.txt', append = TRUE)
print (result_4)
sink()

print(result_4$loadings, cutoff = 0.3) #こちらを評価して、0.3以下なら除外|loading| < 0.3

sink('./EFA_result/nfactors_4_selected.txt', append = TRUE)
print (result_4)
sink()

fa_result <- fa(data, nfactors = 4, rotate = "promax")
print(fa_result$loadings, cutoff = 0.3)

fa_result_loadings=fa_result$loadings
sink('./EFA_result/fa_result_loadings.txt', append = TRUE)
print (fa_result_loadings)
sink()

result = fa( data, nfactors = 4, fm = "minres", rotate = "promax", use = "complete.obs" )
fa.diagram( result )

library(heatmaply)

heatmaply(fa(data, nfactors = 4, fm = "minres", rotate = "promax")$loadings,grid_gap = 1,subplot_widths = c(0.3, 0.2),subplot_heights = c(0.20, 0.70))
