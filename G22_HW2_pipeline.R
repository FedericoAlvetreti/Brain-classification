# Libraries ---------------------------------------------------------------
library(igraph)
library(readr)
library(foreach)
library(iterators)
library(doParallel)
library(ggplot2)
library(e1071)
library(tidyverse)
library(MNITemplate)
library(rgl)
library(vars)
library(caret)
library(graphkernels)
library(Hmisc)
library(sgof)
library(fclust)
library(cluster)
library(stepPlr)
library(tseries)
library(stats)


# setting doParallel parameter
registerDoParallel(cores = detectCores())

source('library_HW2.R')
train <- read.csv('sutff/train_HW03.csv')
test <- read.csv('stuff/test_HW03.csv')

# Preprocessing -----------------------------------------------------------


train$age <- (train$age - min(train$age))/ (max(train$age)-min(train$age))
train$sex <- as.numeric(train$sex=="male")
train$y <- as.numeric(train$y == "autism")
train$male <- train$sex
train$female <- 1-train$sex 
train$sex <- NULL

test$age <- (test$age - min(test$age))/ (max(test$age)-min(test$age))
test$sex <- as.numeric(test$sex=="male")
test$male <- test$sex
test$female <- 1-test$sex 
test$sex <- NULL

train <- find_null_series(train)
test <- correct_null_series(test)

train_prepro <- ROI_train_preprocessing(train)
test_prepro <- ROI_test_preprocessing(test)

# Skip: just load

train <- read_csv("stuff/train_prepro.csv")
test <- read_csv("stuff/test_prepro.csv")

train$...1 <- NULL
test$...1 <- NULL


train <- train %>% drop_na()


# Feature extraction ------------------------------------------------------

all <- rbind(train[,-3],test)
k_opt <- median(best_k_distro(all))                   #Best K for fuzzy cluster
train_features <- get_patients_features(train, k_opt)
test_features <- get_patients_features(test, k_opt)

#Skip: just load 
train_features <- read_csv("stuff/final_train_features.csv")
test_features <- read_csv("stuff/final_test_features.csv")

train$...1 <- NULL
test$...1 <- NULL

train <- train %>% drop_na()
#  Feature assessment and selection----------------------------------------

aut <- which(train_features[,362] == 1)
contr <- which(train_features[,362] == 0)

for (i in 1:361) {
  aut_ <- train_features[aut,i]
  contr_ <- train_features[contr,i]
  test <- wilcox.test(as.numeric(aut_),as.numeric(contr_))
  pvalues <- append(pvalues,test$p.value)
}


relevant <- which(pvalues <=0.05)

# Model validation and training -------------------------------------------

cvfit <-  cv.glmnet(x=data.matrix(train[,relevant]),
                    y=data.matrix(train[,362]), 
                    family = "binomial", 
                    type.measure = "class", 
                    alpha = 1,
                    nfolds = 10)

log <- glmnet(x=data.matrix(train[,relevant]),
              y=data.matrix(train[,362]), 
              family = "binomial",
              alpha = 1,
              lambda=cvfit$lambda.min)


# LOCO assessment ---------------------------------------------------------

train_relevant <- train[,c(relevant,362)]
loco_train <- LOCO(train_relevant)

killed <- which(log$beta==0)
loco_train

h_0 <- intersect(which(loco_train$upper > 0),which(loco_train$lower < 0))


plot(loco_train$feature_importance, ylim=c(-0.15,0.15),xlab='Feature index',ylab='LOCO inference')
segments(x0 = loco_train$feature, y0 = loco_train$lower, y1=loco_train$upper)
segments(x0 = h_0, y0 = loco_train$lower[h_0], y1 = loco_train$upper[h_0], col='pink', lwd=2)
points(killed, loco_train$feature_importance[killed], col='purple',cex=1.2, pch=16)
abline(h=0)
legend(x=c(20,50),
       y=c(-0.14,-0.1),
       legend=c('LOCO-wise indifferent','Killed by LASSO'), 
       col=c('pink','purple'),
       y.intersp = 0.1,
       horiz=T,
       text.font=3,
       lty=c(1,NA),
       pch=c(NA,16),
       bty = "n",
       lwd=2)


