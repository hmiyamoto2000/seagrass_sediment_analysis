# BayesLiNGAM for seagrass

getwd()
setwd("/.../Applications/Bayeslingamtest")
setwd("/.../Applications/Bayeslingamtest")
library('fastICA')
source('bayeslingam/main/loud.R')
loud()

df <- read.csv("Seagrass_path_selected.csv",header=T,sep=",")

#seagrass path
x <- df$Seagrass
y <- df$PWY_2941
w <- df$PWY_6471
w1 <- df$PWY_4722
w2 <- df$P241_PWY
d <- data.frame(x1=x, x2=y,x3=w, x4=w1, x5=w2)
result <- greedybayeslingam(d,model='GL')


#グラフの作成
par("mar"=c(1,1,1,1)) 
library(igraph)
par(mfrow=c(2,3))
prob <- round(result$prob,digits=4) * 100
n.node <- max(result$components$node)
case <- nrow(result$DAGs)
node <- result$DAGs[,1:n.node]
res <- as.matrix(result$DAGs[,-c(1:n.node)],nrow=case)
name <-paste("X",1:n.node,sep="")

for(i in order(prob,decreasing=T)[1:6]){
  amat <- matrix(0,n.node,n.node)
  index <- node[i,]
  amat[lower.tri(amat,diag=FALSE)] <- res[i,]
  amat <- t(amat[index,index])
  g <- graph.adjacency(amat,mode="directed")
  E(g)$label <- NA
  pos <- layout.circle(g)
  rot <- matrix(c(0,1,-1,0),nrow=2,byrow=T)
  la <- pos %*% rot
  if(n.node == 2)la <- matrix(c(-1,0.5,1,0.5),,nrow=2,byrow=T)
  plot(g,layout=la,edge.arrow.size = 0.5,vertex.size=30,vertex.color = "yellow",vertex.label=name) #70を30に変更など
  mtext(paste(prob[i],"%"),line=0)
}


#file取り出し
sink('----------.txt', append = TRUE)
print (result)
sink()

