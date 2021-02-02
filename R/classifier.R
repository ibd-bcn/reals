library(gdata)
library(ggplot2)
library(ggpubr)
library(plyr)
library(ggthemes)
library(Rmisc, quietly =T)
library(gridExtra)
library(grid)
library(e1071)
library("caret")


# Dades
dades <- readxl::read_xlsx('processed/AU_markers.xlsx')


d2 <- pivot_wider(dades, values_from = AU, names_from = Target,
                  values_fn = list(AU = mean))
upa <- d2[d2$Study == "UPA" & d2$Time == "w0", ]


inTrain <- createDataPartition(
  y = upa$remission,
  ## the outcome data are needed
  p = .75,
  ## The percentage of data in the
  ## training set
  list = FALSE
)
training_upa <- upa[ inTrain,]
testing_upa  <- upa[-inTrain,]

# Following https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/
control <- trainControl(method = "repeatedcv", number = 10, repeats = 10,
                        classProbs = TRUE)
plsFit_upa <- train(
  remission~AQP7+CHI3L1+COL3A1+DERL3+GZMH+HTR3E+`IFN GAMMA`+IL17A
  +OSM+PDGFD+PTGDR2+RTNLB+S100A8+SOX6+TBX21+THY1,
  data = training_upa,
  tuneLength = 15,
  method = "rf",
  trControl = control,
  metric = "Accuracy",
  na.action = na.omit
)
plot(plsFit_upa, type = c("g", "o"), main = "UPA")

# estimate variable importance
importance <- varImp(plsFit_upa, scale=FALSE)
plot(importance, main = "UPA")


tnf <- d2[d2$Study == "TNF" & d2$Time == "w0", ]

inTrain <- createDataPartition(
  y = tnf$remission,
  ## the outcome data are needed
  p = .75,
  ## The percentage of data in the
  ## training set
  list = FALSE
)
training_tnf <- tnf[ inTrain,]
testing_tnf  <- tnf[-inTrain,]

# Following https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/
plsFit_tnf <- train(
  remission~AQP7+CHI3L1+COL3A1+DERL3+GZMH+HTR3E+`IFN GAMMA`+IL17A
  +OSM+PDGFD+PTGDR2+RTNLB+S100A8+SOX6+TBX21+THY1,
  data = training_tnf,
  method = "rf",
  trControl = control,
  tuneLength = 15,
  metric = "Accuracy",
  na.action = na.omit
)

# estimate variable importance
importance <- varImp(plsFit_tnf, scale=FALSE)
plot(plsFit_tnf, type = c("g", "o"), main = "TNF")
plot(importance, main = "TNF")
