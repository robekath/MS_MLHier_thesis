#________________________________________________________________________
# 
#  THESIS: Random Forest (RF)
#
#  Author: Katie Roberts
#  Last edited: 8/28/2017
# 
# 
#  Goals: 
#  1. Create function to sample one level2 per level1 (giving unique/uncorrelated records)
#     Steps:
#         i) Randomly sample 1 record per level1 data
#        ii) Grow a forest of (500) trees
#       iii) Predict outcome using forest
#        iv) Capture the predicted responses along with both hierarchy levels
#         v) Repeat steps i-v 100 times - Will now have multiple predicted responses per unique level1 and level2
#        vi) Take majority vote per unique level1 and level2 to get final outcome
#  2. Create testing iris data hierarchy
#  3. Test function with manipulated iris data
#         i) Look at results and confusion matrix
#  4. Additional ROC Curve code
#
# Packages used:
# "randomForest", "caret", "plyr" 
# Note from randomForest() documentation:
# "Any ties are broken at random, so if this is undesirable, avoid it by using odd number ntree in randomForest()."
#
# Citations for packages used and modified here:
#
# randomForest:
#   A. Liaw and M. Wiener (2002). Classification and Regression by randomForest. R News
#     2(3), 18--22.
#
# plyr:
#  Hadley Wickham (2011). The Split-Apply-Combine Strategy for Data Analysis. Journal
#     of Statistical Software, 40(1), 1-29. URL http://www.jstatsoft.org/v40/i01/.
#
# caret:
#  Max Kuhn. Contributions from Jed Wing, Steve Weston, Andre Williams, Chris Keefer,
#     Allan Engelhardt, Tony Cooper, Zachary Mayer, Brenton Kenkel, the R Core Team,
#     Michael Benesty, Reynald Lescarbeau, Andrew Ziem, Luca Scrucca, Yuan Tang, Can
#     Candan and Tyler Hunt. (2016). caret: Classification and Regression Training. R
#     package version 6.0-73. https://CRAN.R-project.org/package=caret
#
#________________________________________________________________________
#
#Packages needed
list.of.packages <- c("randomForest", "plyr", "caret")
#Install packaged if they are not already installed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
#Load packages
lapply(list.of.packages, require, character.only = TRUE)

#----------------------------------------------------------------
# 1. Function to sample one level2 per level1
#----------------------------------------------------------------
#formula: this is the model to evaluate. Use formula info from randomForest()
#data: this is the data you intend to use for the function. Must include the outcome variable, first level hierarchical variable, ...
#x.rep: this numerical value indicates how many different de-correlated samples you want to create from the original dataset. Default is 100.
#L1hier: this is the first level hierarchical variable from the dataset (patient_num)
#L2hier: this is the second level hierarchical variable from the dataset (node)

mult.rf <- function(formula = formula(data), data = sys.frame(sys.parent()), L1hier,L2hier, x.rep=100,...){
  # #the following code chunk meant for testing
  # formula = Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width
  # data=my.iris
  # L1hier="level1"
  # L2hier="level2"
  # x.rep=25
  # outcome <- formula[[2]]
  # outcome.loc <- which(names(data)==outcome) #column location for L1Hier
  
  l1h.loc <- which(names(data)==L1hier) #column location for L1Hier
  l2h.loc <- which(names(data)==L2hier) #column location for L2Hier
  outcome <- formula[[2]] #capture the response variable name
  
  #create placeholder lists to store data for each iteration
  outdata = list()
  outrfmodel = list()
  #outplot = list()
  outimportance = list()
  outrfconf = list()
  outconfbyClass = list()
  
  for(i in 1:x.rep){
    #sample one level2 per level1
    require(plyr)
    newdata <- ddply(data,.(data[,l1h.loc]),function(x) x[sample(nrow(x),1),])

    #Grow a (500) tree forest
    library(randomForest)
    samp.rf <- randomForest(formula, data=newdata,ntrees=500,importance=TRUE)
    outrfmodel[[i]] <- samp.rf #save the samp.rf for each of the x iterations
    #outplot[[i]] <- plot(samp.rf) #save plot for each x iteration for OOB error
    outimportance[[i]] <- importance(samp.rf) #save var importance for each iteration
    
    # Predicting response variable on testing data
    newdata$predicted.response <- predict(samp.rf)
    newdata$predicted.probs <- predict(samp.rf,type='prob')
    
    # Capture the confusion matrix and the matrix info by class
    outcome.loc <- which(names(newdata)==outcome) #column location for outcome/response
    library(caret)
    outrfconf[[i]] <- confusionMatrix(data=newdata$predicted.response,
                                      reference=newdata[,outcome.loc])
    outconfbyClass[[i]] <- outrfconf[[i]]$byClass
    
    #subset into data with only L1hier, L2hier, and rf.outcomes
    tousedata <- newdata[,c(L1hier, L2hier, "predicted.response", "predicted.probs")]
    tousedata$repnum <- i
    outdata[[i]] <- tousedata
  } #end for loop
  
  simple_data = do.call(rbind, outdata) #should contain index, level1 info, level2 info, and predicted.response
  l1h.loc.b <- which(names(simple_data)==L1hier)
  l2h.loc.b <- which(names(simple_data)==L2hier)
  sorted.simple_data <- simple_data[with(simple_data, order(simple_data[,l1h.loc.b], simple_data[,l2h.loc.b])), ]
  
  #loop through the unique combinations of level1 and node to come up with the majority vote and output a dataset with the majority vote predictions from all the forests
  final.data <- data.frame(matrix(ncol = 4, nrow = 0))
  x <- c("l1h", "l2h", "maj.vote.pred", "avg.prob")
  colnames(final.data) <- x
  
  l1h.loc.c <- which(names(sorted.simple_data)==L1hier)
  l2h.loc.c <- which(names(sorted.simple_data)==L2hier)
  
  m=1 #set starting index for rows in final data
  
  for (i in unique(sorted.simple_data[,l1h.loc.c])) {
    #i=1 #for testing
    uniq.data <- sorted.simple_data[(sorted.simple_data[,l1h.loc.c]==i),]
    for (j in unique(uniq.data[,l2h.loc.c])){
      uniq.data.b <- uniq.data[(uniq.data[,l2h.loc.c]==j),]
      #majority vote
      maj.vote.pred <- names(which.max(table(uniq.data.b$predicted.response))) 
      
      #average of predicted probs Metastatic
      avg.prob <- names(which.max(table(t(uniq.data.b$predicted.probs[,1]))))
      avg.prob <- mean(uniq.data.b$predicted.probs[,1])
      
      l1h.loc.d <- which(names(uniq.data.b)==L1hier)
      l2h.loc.d <- which(names(uniq.data.b)==L2hier)
      
      final.data[m,]$l1h <- i
      final.data[m,]$l2h <- j
      final.data[m,]$maj.vote.pred <- maj.vote.pred
      final.data[m,]$avg.prob <- avg.prob
      
      final.data <- rbind(final.data, final.data[m,])
      
      m=m+1
    }# end level 2 loop
    
  }#end level 1 loop
  
  rownames(final.data) <- NULL
  final.data.YES <- unique(final.data)
  final <- as.data.frame(final.data.YES)
  
  final.output <- merge(x=data, y=final, by.x=c(L1hier,L2hier), by.y=c("l1h", "l2h"), all=T)
  
  # return list of lists
  # outList <- list("Final.dat" = final.output, "RF.models" = outrfmodel, "pred.dats"=outdata, "oob.plot"=outplot, "var.imp"=outimportance, "confusion"=outrfconf, "conf.byClass"=outconfbyClass)
  
  outList <- list("Final.dat" = final.output, "RF.models" = outrfmodel, "pred.dats"=outdata, "var.imp"=outimportance, "confusion"=outrfconf, "conf.byClass"=outconfbyClass)
  
  return(outList)
  
}#end function


#-----------------------------------------------------------------
# 2. Create testing iris data hierarchy
#-----------------------------------------------------------------
# NOTE: this hierarchy is generated to test the functionality of the code ONLY. There is no correlation here so the results should be less than interesting
data(iris)
my.iris <- iris
table(my.iris$Species)
my.iris$l1 <- c(rep(1:50,each=3))
my.iris <- my.iris[with(my.iris, order(l1)), ]
my.iris
set.seed(165)
my.iris <- my.iris[sample(nrow(my.iris), 100), ]
#sort data by level1 and level2
my.iris <- my.iris[order(my.iris$l1),]
my.iris$level2 <- ave(my.iris$l1, my.iris$l1, FUN = seq_along)
my.iris
#need to randomize the order of l1 and node so that the outcome (species) will be spread over l1's
l1l2 <- my.iris[6:7]
l1l2
names(l1l2) <- c("l1","l2")
set.seed(165)
ran <- l1l2[sample(nrow(l1l2)),]
my.iris <- cbind(my.iris, ran)
my.iris <- my.iris[,-c(6:7)]
names(my.iris)[names(my.iris)=="l1"] <- "level1"
names(my.iris)[names(my.iris)=="l2"] <- "level2"
my.iris <- my.iris[with(my.iris, order(level1, level2)), ]
my.iris


#----------------------------------------------------------------
# 3. Test function with manipulated iris data
#----------------------------------------------------------------
#formula: this is the model to evaluate. Use formula info from randomForest()
#data: this is the data you intend to use for the function. Must include the outcome variable, first level hierarchical variable, ...
#x.rep: this numerical value indicates how many different de-correlated samples you want to create from the original dataset. Default is 100.
#L1hier: this is the first level hierarchical variable from the dataset (level1)
#L2hier: this is the second level hierarchical variable from the dataset (level2)

head(my.iris)

rf.iris.hier <- mult.rf(Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width, 
                   data=my.iris, 
                   L1hier="level1",
                   L2hier="level2", x.rep=25)
rf.iris.hier[1] #just gets out the final data with majority vote predictions
str(rf.iris.hier[1]$Final.dat) #structure of Final data with majority vote column
str(rf.iris.hier[2]$RF.models[1]) #structure of the first RF model (of x.rep models)
str(rf.iris.hier[3]$pred.dats[1]) #structure of Pred data for first model
str(rf.iris.hier[4]$var.imp[1]) #structure of variable importance for first model
str(rf.iris.hier[5]$confusion[1]) #structure of confusion matrix for first model
str(rf.iris.hier[6]$conf.byClass[1]) #structure of confusion matrix by class for first model

rf.iris.dat <- rf.iris.hier[1]$Final.dat
rf.iris.dat$maj.vote.pred <- as.factor(rf.iris.dat$maj.vote.pred)
#rf.iris.dat$Species <- ordered(rf.iris.dat$Species, levels = c("virginica","setosa", "versicolor"))


#confusion matrix
table(rf.iris.dat$maj.vote.pred, rf.iris.dat$Species)
library(caret)
rf.confusion <- confusionMatrix(data=rf.iris.dat$maj.vote.pred,
                                reference=rf.iris.dat$Species)
rf.confusion
rfmodel.accuracy<-rf.confusion$overall[1]
rfmodel.accuracy

#-------------------------------------------------------------------------
# Build ROC Curve
# NOTE AGAIN: The Iris data is fabricated leading to uninformative results
#-------------------------------------------------------------------------
apost.prob <- rf.iris.hier[1]$Final.dat$avg.prob #for ROC Curve
#ROC curve
library(pROC)
# Plot:
rocobj <- plot.roc(my.iris$Species, apost.prob, print.auc=T, ci=T,
                   print.thres="best", print.thres.best.method="youden")
title(main = "Random Forest ROC Curve - Hier", line = 3)


#Further information:
optimals <- data.frame(t(coords(rocobj,best.method="youden", rocobj$thresholds, "thr", ret=c("thr","se", "sp"))))
rownames(optimals) <- NULL
optimals$youden <- optimals$sensitivity + optimals$specificity -1
#optimals
maxyouden <- optimals[which(optimals$youden == max(optimals$youden)), ]
maxyouden









