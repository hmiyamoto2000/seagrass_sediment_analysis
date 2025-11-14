
#install.packages("vegan")
library(vegan)
pc = read.csv("Bacteria_diversity.csv", header= TRUE)
dim(pc) #[1] 132 997
#make community matrix - extract columns with abundance information
com = pc[,2:997]
m_com = as.matrix(com)
set.seed(123)
nmds = metaMDS(m_com, distance = "bray")
nmds
plot(nmds)

#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
data.scores = as.data.frame(scores(nmds)$sites)

write.csv(data.scores,"NMDS.csv") #一行目にTypeと記載し、分類を記載

values=read.csv("NMDS.csv")

library(pairwiseAdonis)
values=read.csv("List.csv")

library(pairwiseAdonis)
p0=pairwise.adonis(pc[,2:997],values$cnd)
sink('NMDS_pairwiseR_cnd.txt', append = TRUE)　
print (p0)
sink()

library(pairwiseAdonis)
p0=pairwise.adonis(pc[,2:997],values$env)
sink('NMDS_pairwiseR_envtxt', append = TRUE)　#Comp ThB別
print (p0)
sink()






