#RiboMNaseML_learn.R
#Predict P sites using Random Forest trained model
#Reference: https://github.com/mvanins/scRiboSeq_manuscript/blob/main/random_forest/tuneModel_CVpred.R
#Version: Yu Sun, 2021-11

print("Loading packages...")
suppressMessages(library(dplyr))
suppressMessages(library(mlr))
suppressMessages(library(mlr3))
suppressMessages(library(ranger))
suppressMessages(library(ggplot2))
suppressMessages(library(forcats))
suppressMessages(library(parallelMap))
suppressMessages(library(pheatmap))

args <- commandArgs(TRUE)
File <- args[1]
ModelFile <- args[2]
Output <- args[3]

TunedModel <- readRDS(ModelFile)

ColName35Fea <- c("base5m1", "base5m2", "base5m3", "base5m4", "base5m5", "base5m6", "base5m7", "base5m8","base5",
                  "base5p1", "base5p2", "base5p3", "base5p4", "base5p5", "base5p6", "base5p7", "base5p8",
                  "base3m1", "base3m2", "base3m3", "base3m4", "base3m5", "base3m6", "base3m7", "base3m8","base3",
                  "base3p1", "base3p2", "base3p3", "base3p4", "base3p5", "base3p6", "base3p7", "base3p8",
                  "length")

print("Reading data...")
Data <- read.table(File, col.names = ColName35Fea)
Data <- Data %>% mutate_if(is.character, as.factor)
Data$p_offset <- as.factor(1) 

#Create task
PredictionTask <- makeClassifTask(data = Data, target = "p_offset")

#Train using tuned parameters
print("Making predictions...")
tuned.predict.test <- predict(TunedModel, PredictionTask)
Results <- tuned.predict.test$data$response
Results <- as.numeric(as.character(Results))
print("Writing outputs...")
write(as.numeric(as.character(Results)), file = Output, ncolumns = 1)
print("Prediction workflow done.")
