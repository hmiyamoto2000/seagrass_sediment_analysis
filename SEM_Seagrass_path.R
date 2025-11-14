#For Shapiro-Wilk test and SEM
dir.create("./SEM_result", showWarnings = TRUE, recursive = FALSE, mode = "0777")

library(MVN)

input <- read.csv("Seagrass_path_selected.csv")#251009

mvn_data=mvn(input[,-1])

sink('./SEM_result/mvn_data.txt', append = TRUE)
print (mvn_data)
sink()

#path
createModel10 <- function() {
  +     # variable ∼ dependent variable
    +     # + dependent variable
    return ("
    PWY_2941 + P241_PWY + PWY_4722 + PWY_6471 ~ Seagrass
    PWY_2941 ~ PWY_4722
    P241_PWY ~ PWY_6471
                ")
}

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
         negCol="violet", fade=FALSE, edge.label.cex=0.8) #Type1

semPaths(res.l10, what="stand", layout="tree3", style="lisrel",
         shapeMan="rectangle", shapeLat="ellipse",
         sizeMan=6, residScale=5, posCol="darkgreen",
         negCol="violet", fade=FALSE, edge.label.cex=0.8) #Type1

semPaths(res.l10, what="stand", layout="circle", style="lisrel",
         shapeMan="rectangle", shapeLat="ellipse",
         sizeMan=6, residScale=5, posCol="darkgreen",
         negCol="violet", fade=FALSE, edge.label.cex=0.8) #Type1

semPaths(res.l10, what="stand", layout="spring", style="lisrel",
         shapeMan="rectangle", shapeLat="ellipse",
         sizeMan=6, residScale=5, posCol="darkgreen",
         negCol="violet", fade=FALSE, edge.label.cex=0.8) #Type1

semPaths(res.l10, what="stand", layout="tree", style="lisrel",
         shapeMan="rectangle", shapeLat="ellipse",
         sizeMan=8, residScale=5, posCol="darkgreen",
         negCol="violet", fade=FALSE, edge.label.cex=1.0) #Type2

semPaths(res.l10, what="stand", layout="tree", style="lisrel",
         shapeMan="rectangle", shapeLat="ellipse",
         sizeMan=10, residScale=5, posCol="darkgreen",
         negCol="violet", fade=FALSE, edge.label.cex=1.0) #gType3
semPaths(res.l10, what="stand", layout="groups", style="lisrel",
         shapeMan="rectangle", shapeLat="ellipse",
         sizeMan=10, residScale=5, posCol="darkgreen",
         negCol="violet", fade=FALSE, edge.label.cex=1.0) #gType3

semPaths(res.l10, what="stand", layout="groups", style="lisrel",
         shapeMan="rectangle", shapeLat="ellipse",
         sizeMan=8, residScale=5, posCol="darkgreen",
         negCol="violet", fade=FALSE, edge.label.cex=1.0) #Type2






