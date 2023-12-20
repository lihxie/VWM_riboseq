#RiboMNaseML_learn.R
#Train Random Forest model using reads spanning STOP codon
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
Prefix <- args[2]

#Default value:
mtry_tuned <- 13
numtrees_tuned <- 500

ColName35Fea <- c("base5m1", "base5m2", "base5m3", "base5m4", "base5m5", "base5m6", "base5m7", "base5m8","base5",
                  "base5p1", "base5p2", "base5p3", "base5p4", "base5p5", "base5p6", "base5p7", "base5p8",
                  "base3m1", "base3m2", "base3m3", "base3m4", "base3m5", "base3m6", "base3m7", "base3m8","base3",
                  "base3p1", "base3p2", "base3p3", "base3p4", "base3p5", "base3p6", "base3p7", "base3p8",
                  "length","p_offset")

Data <- read.table(File, col.names = ColName35Fea)
Data <- Data %>% mutate_if(is.character, as.factor)    #Change all character to factor, then rewrite Data. makeClassifTask doesn't allow characters
Data$p_offset <- as.factor(Data$p_offset)              #Make the outcome to be factor also

train_set <- sample(nrow(Data), size = nrow(Data)*.8)
test_set <- (1:nrow(Data))[-train_set]
trainReads <- Data[train_set,]

#Create tasks
TrainTask <- makeClassifTask(data = trainReads, target = "p_offset")  #This is for hyperparameter tuning
AssessTask <- makeClassifTask(data = Data, target = "p_offset")       #For assessment

#Basic Learner
RangerLearner <- makeLearner("classif.ranger", predict.type = "response")

#Search all hyperparameters using the given range.
#This is time-consuming. Better to use multiple cores.
#Hyperparameters: num.trees: Number of trees
#                 mtry: Number of variables randomly sampled as candidates at each split.
#                 min.node.size: Minimal data in a terminal node, usually default for classification (1) and regression (5).
#Use the following default if the training step fails
#[Tune] Result: num.trees=500; mtry=13; min.node.size=1 : acc.test.mean=0.8116545
parallelStartSocket(parallel::detectCores())
res <- tuneParams(RangerLearner,
                  task = TrainTask,
                  resampling = makeResampleDesc("CV", iters = 5L),
                  par.set = makeParamSet(makeIntegerParam("num.trees", lower = 500, upper = 500),
                                         makeIntegerParam("mtry", lower = 13, upper = 27),
                                         makeIntegerParam("min.node.size", lower = 1, upper = 5)),
                  control = makeTuneControlRandom(maxit=100),
                  measure = acc)
parallelStop()

mtry_tuned <- res$x$mtry
numtrees_tuned <- res$x$num.trees

print(paste0("Use tuned mtry value: ",mtry_tuned))
print(paste0("Use tuned numtrees value: ",numtrees_tuned))

RangerLearner_Tuned <- makeLearner("classif.ranger",
                                    predict.type = "response",
                                    importance = 'permutation',
                                    respect.unordered.factors = 'partition',
                                    mtry = mtry_tuned,
                                    num.trees = numtrees_tuned)

#Train using tuned parameters
TunedModel <- train(RangerLearner_Tuned, AssessTask, subset = train_set)
TestResults <- predict(TunedModel, AssessTask, subset = test_set)
print("Accuracy on test dataset: ")
performance(TestResults, measures=list(acc))

modelInfo <- data.frame(importance = TunedModel$learner.model$variable.importance,
                        names = as.factor(names(TunedModel$learner.model$variable.importance)))

print("Plotting figures...")
pdf(paste0(Prefix,".RFModel.35features.importance.pdf"),width = 30, height = 9)
barplot(modelInfo$importance, names.arg = modelInfo$names,ylab = "Permutation importance",
        main = "Base importance using sequences surrounding 5end",
        cex.names = 0.8)
dev.off()

pdf(paste0(Prefix,".RFModel.35features.importance.sorted.pdf"),width = 30, height = 9)
modelInfo_order <- modelInfo[order(modelInfo$importance, decreasing = T),]
barplot(modelInfo_order$importance, names.arg = modelInfo_order$names,ylab = "Permutation importance",
        main = "Base importance using sequences surrounding 5end, sorted",
        cex.names = 0.8)
dev.off()

ConfusionMatrix <- calculateConfusionMatrix(TestResults)
matrix <- ConfusionMatrix$result
matrix <- matrix[1:ncol(matrix)-1,1:ncol(matrix)-1]

pdf(paste0(Prefix,".RFModel.35features.heatmap.pdf"),width = 10, height = 10)
pheatmap(matrix, cluster_rows=FALSE, cluster_cols=FALSE)
dev.off()

saveRDS(TunedModel, file = paste0(Prefix,".TunedModel.rds"))