#For Shapiro-Wilk test and SEM
dir.create("./SEM_result", showWarnings = TRUE, recursive = FALSE, mode = "0777")

library(MVN)
input <- read.csv("EFA_seagrass.csv")

mvn(input[,-1])

mvn_data=mvn(input[,-1])

sink('./SEM_result/mvn_data.txt', append = TRUE)
print (mvn_data)
sink()

#For SEM
createModel10 <- function() {
  +     # variable ∼ dependent variable
    +     # + dependent variable
    return ("
    Seagrass ~ B_Desulfobulbaceae
    B_Desulfobulbaceae + Seagrass ~ B_Hyphomonadaceae + E_Corallinophycidae + E_Diatomea
    B_Hyphomonadaceae ~ E_Corallinophycidae
                ")
}
#Top srmr0.011 gfi 0.999 agfi 0.984 rmsea0 pvalue  0.639

#createModel10 <- function() {
#  +     # variable ∼ dependent variable
#    +     # + dependent variable
#    return ("
#    B_Desulfobulbaceae ~ Seagrass
#    B_Desulfobulbaceae + Seagrass ~ B_Hyphomonadaceae + E_Corallinophycidae + E_Diatomea
#    B_Hyphomonadaceae ~ E_Corallinophycidae
#                ")
#}
#2nd  srmr0.011 gfi 0.999 agfi 0.984 rmsea0 pvalue  0.639

#createModel10 <- function() {
#  +     # variable ∼ dependent variable
#    +     # + dependent variable
#    return ("
#    Seagrass ~ B_Desulfobulbaceae
#    B_Desulfobulbaceae + Seagrass ~ B_Hyphomonadaceae + E_Corallinophycidae + E_Diatomea
#    E_Corallinophycidae ~ B_Hyphomonadaceae
#                ")
#}
#3rd srmr 0.008 gfi 0.999 agfi 0.992 pvalue0.737

model.l10 <- createModel10()

library(lavaan)

res.l10 <- lavaan(model.l10, data = input, estimator = "MLR", auto.var = TRUE)
summary(res.l10, fit.measure=TRUE,　standardized = TRUE)
fitMeasures(res.l10)

p0=model.l10
sink('./SEM_result/lavaanstaticsmodel.l10.txt', append = TRUE)
print (p0)
sink()

p1=summary(res.l10, fit.measure=TRUE,　standardized = TRUE)
sink('./SEM_result/lavaanstaticssummary.txt', append = TRUE)
print (p1)
sink()

p2=fitMeasures(res.l10)
sink('./SEM_result/lavaanstaticsfitMeasures.txt', append = TRUE)
print (p2)
sink()

#optinal evaluation:|SR|<2.58
sr <- residuals(res.l10, type = "standardized")
sink('./SEM_result/residuals.txt', append = TRUE)
print (sr)
sink()

# bootstrap based on ML

res.l12 <- lavaan(model.l10, data = input, auto.var = TRUE,se="bootstrap",bootstrap=1000)
summary(res.l12, fit.measure=TRUE,　standardized = TRUE)
P4=summary(res.l12, fit.measure=TRUE,　standardized = TRUE)
fitMeasures(res.l12)

sink('./SEM_result/lavaanstaticsmodel.lbootres12.txt', append = TRUE)
print (P4)
sink()

#Visualize the SEM
library(OpenMx)
library(psych)
library(semPlot)
library(rgl)

semPaths(res.l10, what="stand", layout="tree2", style="lisrel",
         shapeMan="rectangle", shapeLat="ellipse",
         sizeMan=6, residScale=5, posCol="darkgreen",
         negCol="violet", fade=FALSE, edge.label.cex=0.8) #Type1

semPaths(res.l10, what="stand", layout="tree", style="lisrel",
         shapeMan="rectangle", shapeLat="ellipse",
         sizeMan=6, residScale=5, posCol="darkgreen",
         negCol="violet", fade=FALSE, edge.label.cex=0.8) #Type2

semPaths(res.l10, what="stand", layout="tree3", style="lisrel",
         shapeMan="rectangle", shapeLat="ellipse",
         sizeMan=6, residScale=5, posCol="darkgreen",
         negCol="violet", fade=FALSE, edge.label.cex=0.8) #Type3

semPaths(res.l10, what="stand", layout="circle", style="lisrel",
         shapeMan="rectangle", shapeLat="ellipse",
         sizeMan=6, residScale=5, posCol="darkgreen",
         negCol="violet", fade=FALSE, edge.label.cex=0.8) #Type4

semPaths(res.l10, what="stand", layout="spring", style="lisrel",
         shapeMan="rectangle", shapeLat="ellipse",
         sizeMan=6, residScale=5, posCol="darkgreen",
         negCol="violet", fade=FALSE, edge.label.cex=0.8) #Type5

semPaths(res.l10, what="stand", layout="tree", style="lisrel",
         shapeMan="rectangle", shapeLat="ellipse",
         sizeMan=8, residScale=5, posCol="darkgreen",
         negCol="violet", fade=FALSE, edge.label.cex=1.0) #Type6

semPaths(res.l10, what="stand", layout="tree", style="lisrel",
         shapeMan="rectangle", shapeLat="ellipse",
         sizeMan=10, residScale=5, posCol="darkgreen",
         negCol="violet", fade=FALSE, edge.label.cex=1.0) #gType7
semPaths(res.l10, what="stand", layout="groups", style="lisrel",
         shapeMan="rectangle", shapeLat="ellipse",
         sizeMan=10, residScale=5, posCol="darkgreen",
         negCol="violet", fade=FALSE, edge.label.cex=1.0) #gType8

semPaths(res.l10, what="stand", layout="groups", style="lisrel",
         shapeMan="rectangle", shapeLat="ellipse",
         sizeMan=8, residScale=5, posCol="darkgreen",
         negCol="violet", fade=FALSE, edge.label.cex=1.0) #Type9

