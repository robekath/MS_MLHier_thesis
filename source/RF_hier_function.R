#________________________________________________________________________
# 
#  THESIS: Random Forest (RF)
#
#  Author: Katie Roberts
#  Last edited: 1/13/2017
# 
# 
#  Goals: 
#  1. Create function to sample one level2 per level1 (giving unique/uncorrelated records)
#     Steps:
#         i) Randomly sample 1 record per level1 data
#        ii) Use sample to create training(60%) and testing(40%) data
#       iii) Grow a forest of (500) trees on training set
#        iv) Predict outcome of testing set using forest
#         v) Capture the predicted responses on the testing set along with both hierarchy levels
#        vi) Repeat steps i-v 100 times - Will now have multiple predicted responses per unique level1 and level2
#       vii) Take majority vote per unique level1 and level2 to get final outcome
#  2. Create testing iris data hierarchy
#  3. Test function with manipulated iris data
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
#install.packages("randomForest")
#install.packages("plyr")
#install.packages("caret")

library(randomForest)
library(plyr)
library(caret)


#----------------------------------------------------------------
# 1. Function to sample one level2 per level1
#----------------------------------------------------------------
#formula: this is the model to evaluate. Use formula info from randomForest()
#data: this is the data you intend to use for the function. Must include the outcome variable, first level hierarchical variable, ...
#x.rep: this numerical value indicates how many different de-correlated samples you want to create from the original dataset. Default is 100.
#L1hier: this is the first level hierarchical variable from the dataset (patient_num)
#L2hier: this is the second level hierarchical variable from the dataset (node)

mult.rf <- function(formula = formula(data), data = sys.frame(sys.parent()), L1hier,L2hier, x.rep=100,...){
  
  l1h.loc <- which(names(data)==L1hier) #column location for L1Hier
  l2h.loc <- which(names(data)==L2hier) #column location for L2Hier
  outcome <- formula[[2]] #capture the response variable name
  
  #create placeholder lists to store data for each iteration
  outdata = list()
  outrfmodel = list()
  outplot = list()
  outimportance = list()
  outrfconf = list()
  outconfbyClass = list()
  
  for(i in 1:x.rep){
    #sample one level2 per level1
    require(plyr)
    newdata <- ddply(data,.(data[,l1h.loc]),function(x) x[sample(nrow(x),1),])
    #Sample 60/40 training/testing data
    sample.ind <- sample(2, nrow(newdata),replace = T,prob = c(0.6,0.4))
    training <- newdata[sample.ind==1,]
    testing <- newdata[sample.ind==2,]
    
    #Grow a (500) tree forest
    library(randomForest)
    samp.rf <- randomForest(formula, data=training,ntrees=500,importance=TRUE)
    outrfmodel[[i]] <- samp.rf #save the samp.rf for each of the x iterations
    outplot[[i]] <- plot(samp.rf) #save plot for each x iteration for OOB error
    outimportance[[i]] <- importance(samp.rf) #save var importance for each iteration
    
    # Predicting response variable on testing data
    testing$predicted.response <- predict(samp.rf,testing)
    testing$predicted.probs <- predict(samp.rf,testing,type='prob')
    
    # Capture the confusion matrix and the matrix info by class
    outcome.loc <- which(names(testing)==outcome) #column location for outcome/response
    library(caret)
    outrfconf[[i]] <- confusionMatrix(data=testing$predicted.response,
                                      reference=testing[,outcome.loc])
    outconfbyClass[[i]] <- outrfconf[[i]]$byClass
    
    #subset into data with only L1hier, L2hier, and rf.outcomes
    tousedata <- testing[,c(L1hier, L2hier, "predicted.response", "predicted.probs")]
    tousedata$repnum <- i
    outdata[[i]] <- tousedata
  } #end for loop
  
  simple_data = do.call(rbind, outdata) #should contain index, level1 info, level2 info, and predicted.response
  l1h.loc.b <- which(names(simple_data)==L1hier)
  l2h.loc.b <- which(names(simple_data)==L2hier)
  sorted.simple_data <- simple_data[with(simple_data, order(simple_data[,l1h.loc.b], simple_data[,l2h.loc.b])), ]
  
  #loop through the unique combinations of level1 and node to come up with the majority vote and output a dataset with the majority vote predictions from all the forests
  final.data <- data.frame(matrix(ncol = 3, nrow = 0))
  x <- c("l1h", "l2h", "maj.vote.pred")
  colnames(final.data) <- x
  
  l1h.loc.c <- which(names(sorted.simple_data)==L1hier)
  l2h.loc.c <- which(names(sorted.simple_data)==L2hier)
  
  m=1 #set starting index for rows in final data
  
  for (i in unique(sorted.simple_data[,l1h.loc.c])) {
    uniq.data <- sorted.simple_data[(sorted.simple_data[,l1h.loc.c]==i),]
    for (j in unique(uniq.data[,l2h.loc.c])){
      uniq.data.b <- uniq.data[(uniq.data[,l2h.loc.c]==j),]
      #majority vote
      maj.vote.pred <- names(which.max(table(uniq.data.b$predicted.response))) 
      
      l1h.loc.d <- which(names(uniq.data.b)==L1hier)
      l2h.loc.d <- which(names(uniq.data.b)==L2hier)
      
      final.data[m,]$l1h <- i
      final.data[m,]$l2h <- j
      final.data[m,]$maj.vote.pred <- maj.vote.pred
      
      final.data <- rbind(final.data, final.data[m,])
      
      m=m+1
    }# end level 2 loop
    
  }#end level 1 loop
  
  rownames(final.data) <- NULL
  final.data.YES <- unique(final.data)
  final <- as.data.frame(final.data.YES)
  
  final.output <- merge(x=data, y=final, by.x=c(L1hier,L2hier), by.y=c("l1h", "l2h"), all=T)
  
  # return list of lists
  outList <- list("Final.dat" = final.output, "RF.models" = outrfmodel, "pred.dats"=outdata, "oob.plot"=outplot,
                  "var.imp"=outimportance, "confusion"=outrfconf, "conf.byClass"=outconfbyClass)
  
  return(outList)
  
}#end function



#-----------------------------------------------------------------
# 2. Create testing iris data hierarchy
#-----------------------------------------------------------------
# NOTE: this hierarchy is generated to test the functionality of the code ONLY. There is no correlation
# here so the results should be less than interesting
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
str(rf.iris.hier[4]$oob.plot[1]) #structure of OOB error for first model
str(rf.iris.hier[5]$var.imp[1]) #structure of variable importance for first model
str(rf.iris.hier[6]$confusion[1]) #structure of confusion matrix for first model
str(rf.iris.hier[7]$conf.byClass[1]) #structure of confusion matrix by class for first model

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









