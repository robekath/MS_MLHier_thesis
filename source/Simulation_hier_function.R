#________________________________________________________________________
# 
#  THESIS: Simulation_hier_function
#
#  Author: Katie Roberts
#  Last edited: 7/18/2017
# 
# 
#  Steps: 
#  1. Set-up, load, and Initialize packages, data, and functions to use
#  2. Capture Beta's to Use in Simulation
#         i) Pull in all relevant variables from dataset:
#             a) thyroidDataAllNodes
#             b) final model variables (plus age continuous variable) were used
#        ii) Run glmer on model with those variables
#       iii) Collect beta's from model
#        iv) Capture distribution details for variables
#  3. Simulation Function
#         i) Set up lists and vectors
#        ii) Include Function for warnings and errors
#       iii) For each simulation:
#             a) simulate lesions
#             b) expand the dataset to be correct size
#             c) For each person/subject:
#                 1) get their number of nodes (d_i)
#                 2) get variance between subject (B_0i)
#                 3) For each lesion/node
#                       i) simulate covariates using betas and distributions
#                      ii) gather probability (p) for each node
#                     iii) simulate outcome
#             d) Factor certain variables if methods require it
#             e) Partition data sets into 60/40 for trainings/testing
#             f) Fit Models:
#                 1) GLMM
#                         i) subset on variables
#                        ii) Fit model and catch warnings & errors
#                       iii) IF no warnings/errors, then...
#                             a) Fit two-way Naive Bayes model
#                             b) Fit PDA non-hierarchical
#                             c) Fit PDA hierarchical **catch warns/errors
#                             d) Fit RF non-hierarchical
#                             e) Fit RF hierarchical
#                       iii) ELSE IF warnings/errors, then...
#                             a) Capture the warning/error message
#             g) Capture simulation #
#             h) Row bind simulation info (fits, confusion matrices, messages, etc.)
#        iv) Return the list of outputs
#  4. Run a Test of the Simulation and check output
#  5. Run Simulations and save as R datasets
#         i) Sim1: nsim=500, n=100, sd=sd1
#         i) Sim2: nsim=500, n=100, sd=sd2
#         i) Sim3: nsim=500, n=100, sd=sd3
#         i) Sim4: nsim=500, n=200, sd=sd1
#         i) Sim5: nsim=500, n=200, sd=sd2
#         i) Sim6: nsim=500, n=200, sd=sd3
#  6. For cleanup and results see simulation_results.R program
#
# Separate/related programs: 
#     PDA_hier_function.R
#     RF_hier_function.R
#
# Packages used:
# "ggplot2", "LaplacesDemon", "MASS", "lme4", "clusterPower", "caret", "Runuran"
# "pushoverr" used to send message updates to my phone when analysis was completed

#
#________________________________________________________________________
#Packages needed
list.of.packages <- c("ggplot2", "LaplacesDemon", "MASS", "lme4", "clusterPower", "caret", "Runuran", "pushoverr")
#Install packaged if they are not already installed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
#Load packages
lapply(list.of.packages, require, character.only = TRUE)

# Need to source these as well
# https://github.com/robekath/MS_MLHier_thesis/blob/master/source/PDA_hier_function.R
# https://github.com/robekath/MS_MLHier_thesis/blob/master/source/RF_hier_function.R
library(RCurl)
pdafun <- getURL("https://raw.githubusercontent.com/robekath/MS_MLHier_thesis/master/source/PDA_hier_function.R", ssl.verifypeer = FALSE)
invisible(eval(parse(text = pdafun)))
rffun <- getURL("https://raw.githubusercontent.com/robekath/MS_MLHier_thesis/master/source/RF_hier_function.R", ssl.verifypeer = FALSE)
invisible(eval(parse(text = rffun)))


#=================================================================================
# Set up Beta's from GLMM in original analysis
#=================================================================================
# Initialize beta's
  # **these were collected from an original analysis using glmer, hardcoded here
  # glmmfit<-glmer(y~x+(1|patient_num),data=dataToSim,family=binomial(link="logit"),nAGQ=10)
  # (beta_int <- unname(fixef(glmmfit)[1]))
  # (beta_x1 <- unname(fixef(glmmfit)[2]))
  # (beta_x2 <- unname(fixef(glmmfit)[3]))
  # (beta_x3 <- unname(fixef(glmmfit)[4]))
  # (beta_x4 <- unname(fixef(glmmfit)[5]))
  # (sigma_int <- attributes(VarCorr(glmmfit)$patient_num)$stddev)

(beta_int <- -40.4287)
(beta_x1 <- 5.321149)
(beta_x2 <- 12.87692)
(beta_x3 <- 7.989059)
(beta_x4 <- -0.1347207)
(sigma_int <- 12.38476)

#=================================================================================
# Distributions of covariates in the simulation:
#=================================================================================
# The distributions were simulated based on the original dataset variables, but are hardcoded here:
x1_mean <- 5.64672 #x1 Mean
x1_sd <- 1.076254 #x1 sd
x2_prob <- 0.5633803 #x2 success prob
x3_prob <- 0.2394366 #x3 success prob
x4_mean <- 46.5493 #x4 Mean
x4_sd <- 11.11857 #x4 sd

#=================================================================================
# Generic Simulation Function
#=================================================================================
# nsim is number of simulations
# n is the number of SUBJECTS (i.e. Level1 hierarchy) in training AND testing data
# B0, B1, B2, B3, and B4 are the betas for the model
# sd is the between-subject sd for G matrix (bi0)
# kmin is the min amount of lesions a person can have
# kmax is the max amount of lesions a person can have


hiersimfun <- function(nsim, n, B0, B1, B2, B3, B4, sd, kmin=1, kmax=5) {
  # Set up lists and data frames to collect output from simulation
  sim.df = list()
  training.df = list()
  testing.df = list()
  confmat.glmm = list()
  confmat.6040 = list()
  confmat.pda.nh = list()
  confmat.pda.h = list()
  confmat.rf.nh = list()
  confmat.rf.h = list()
  warn.error = data.frame(sim=numeric(), GW=numeric(), Gmess=character(), GE=numeric(), Gerr=character(),
                          PDAHE=numeric(), PDAHerr=character(),
                          NBW=numeric(), NBmess=character(),
                          PDANW=numeric(), PDANmess=character(),
                          PDAHW=numeric(), PDAHmess=character(),
                          RFNW=numeric(), RFNmess=character(),
                          RFHW=numeric(), RFHmess=character())
  
  for(s in 1:nsim){ #for every simulation...
    #simulate number of lesions per person (n)
    k <- c(kmin:kmax)
    d <- sample(x=k, size = n,replace=T) # d is the number of nodes within subject
    #expand a data frame to be the correct size
    sim.df[[s]] <- data.frame(subject = factor(rep.int(x=c(1:n), times=d))) #total rows = sum(d)
    sim.df[[s]]$node <-  ave(sim.df[[s]]$subject, sim.df[[s]]$subject, FUN = seq_along)
    namevector <- c("B0i","x1","x2","x3","x4", "prob","outcome")
    sim.df[[s]][,namevector] <- NA
    
    for(i in 1:n){ #for every subject...
      di <- d[i] #the subject i's number of nodes
      B0i <- rnorm(n=1, mean=B0, sd=sd) #simulates the B0i. Captures the variance between subjects
      sim.df[[s]][sim.df[[s]]$subject %in% i,]$B0i <- rep(B0i, di)

      for(j in 1:di){ #for every lesion...
        #simulate covariates using GLMM distributions.
        # Note: Argument values used in random generation have been created and captured prior to function
        #x1_rv, x2_rv, x2_rv_fact, x3_rv, and x3_rv_fact are not subject-specific
        #x1
        x1_rv <- rnorm(n=1, mean=x1_mean, sd=x1_sd)
        sim.df[[s]][sim.df[[s]]$subject %in% i,][j,]$x1 <- x1_rv
        #x2
        x2_rv <- rbinom(n=1, size = 1, p=x2_prob)
        sim.df[[s]][sim.df[[s]]$subject %in% i,][j,]$x2 <- x2_rv
        #x3
        x3_rv <- rbinom(n=1, size = 1, p=x3_prob)
        sim.df[[s]][sim.df[[s]]$subject %in% i,][j,]$x3 <- x3_rv
        #x4
        x4_rv <- round(urnorm(n=1, mean=x4_mean, sd=x4_sd, lb=25, ub=80), 2)
        sim.df[[s]][sim.df[[s]]$subject %in% i,][j,]$x4 <- x4_rv
        
        #gather probability for subject and node, p
        p <- expit(B0i + (B1*x1_rv) + (B2*x2_rv) + (B3*x3_rv) + (B4*x4_rv))
        sim.df[[s]][sim.df[[s]]$subject %in% i,][j,]$prob <- p
        
        #simulate the outcome y variable for every subject and node
        outcome.cancer <- rbinom(n=1, size = 1, p=p)
        sim.df[[s]][sim.df[[s]]$subject %in% i,][j,]$outcome <- outcome.cancer
        
      } #end for(j in 1:di) lesion loop
    } #end for(i in 1:n) subject loop
    
    #Add factored columns to the data.frame for some methods
    sim.df[[s]]$x2_factor <- factor(sim.df[[s]]$x2, 
                                       levels=c(1,0),
                                       labels=c("Yes", "No"))
    sim.df[[s]]$x3_factor <- factor(sim.df[[s]]$x3, 
                           levels=c(1,0),
                           labels=c("Yes", "No"))
    sim.df[[s]]$outcome_factor <- factor(sim.df[[s]]$outcome, 
                                     levels=c(1,0),
                                     labels=c("Yes", "No"))
    
    #---------------------------------------------
    # data partition: 60% to train, 40% to test
    #---------------------------------------------
    # sample 60% training rows within subject groups
    sub <- unique(sim.df[[s]]$subject)
    samp <- sample(sub, 0.6*length(sub), replace=FALSE)
    training.df[[s]] <- sim.df[[s]][sim.df[[s]]$subject %in% samp , ]
    testing.df[[s]] <- sim.df[[s]][!(sim.df[[s]]$subject %in% samp) , ]
    
    #---------------------------------------------
    # For each of the 6 models:
    #     Fit model to training data
    #     Predict on the testing data
    #     Save confusion matrix
    #---------------------------------------------
    # 1)
    # GLMM, using final model
    glmm.varnames <- c("subject","x1", "x2", "x3", "x4", "outcome")
    tmp.traindf <- training.df[[s]][,glmm.varnames]
    tmp.testdf <- testing.df[[s]][,glmm.varnames]
    
    #catch fit, warnings, and errors
    catchToList <- function(expr) {
      val <- NULL
      ## warnings 
      myWarnings <- NULL # init warning
      wHandler <- function(w) {
        myWarnings <<- c(myWarnings, w$message)
        invokeRestart("muffleWarning")
      }
      val <- tryCatch(withCallingHandlers(expr, warning = wHandler))
      ## return result with warnings and errors 
      list(value = val, warnings = myWarnings)
    }
    
    possibleError <- tryCatch(
      capture <- catchToList(glmer(outcome~ x1 + x2 + x3 + x4  + (1|subject),
                                   data=tmp.traindf,family=binomial(link = "logit"), nAGQ=10)),
      error = function(e) e
    )
    
    #### NOTE: THIS CODE ASSUMES ONLY GLMM WILL PRODUCE AN ERROR
    
    if(!inherits(possibleError, "error")){ #if there is no error
      Gerr <- paste("No error")
      GE <- 0 #1 if there is an error, 0 if no error
      
      tmpGlmmfit <- capture$value #tmpGlmmfit
      Gmess <- ifelse(is.null(capture$warnings), "No warnings",capture$warnings)
      GW <- ifelse(!is.null(capture$warnings),1,0) #1 if there is a warning, 0 if no warning
      # predict the test value
      tmp.pred = predict(tmpGlmmfit, newdata=tmp.testdf, allow.new.levels = TRUE, type="response")
      # Use cutoff from original GLMM work
      # my work = 0.26
      # [note: article = 0.38]
      tmp.pred2 <- ifelse(tmp.pred < 0.26,"No","Yes")
      tmp.test = cbind(predProb=tmp.pred,tmp.pred2, tmp.testdf)   
      tmp.test$outcome_factor <- factor(tmp.test$outcome, 
                                       levels=c(0,1),
                                       labels=c("No","Yes"))
      # results
      confmat.glmm[[s]] <- confusionMatrix(tmp.pred2, tmp.test$outcome_factor, positive='Yes')
      
      # 2)
      # Naive bayes, 60/40
      nb.varnames <- c("outcome_factor", "x1", "x2_factor", "x4", "x3_factor")
      tmp.traindf <- training.df[[s]][,nb.varnames]
      tmp.testdf <- testing.df[[s]][,nb.varnames]
      library(e1071)
      #catch fit and warnings
      capture <- catchToList(naiveBayes(outcome_factor ~ ., data = tmp.traindf, laplace=1))
      tmpNBfit.6040 <- capture$value
      NBmess <- ifelse(is.null(capture$warnings), "No warnings",capture$warnings)
      NBW <- ifelse(!is.null(capture$warnings),1,0) #1 if there is a warning, 0 if no warning
      # predict the test value
      tmp.pred <- predict(tmpNBfit.6040, tmp.testdf[,-1])
      tmp.test = cbind(tmp.pred, tmp.testdf)   
      # results
      confmat.6040[[s]] <- confusionMatrix(tmp.pred, tmp.test$outcome_factor, positive='Yes')
      
      # 3)
      # PDA, ignore hierarchy
      pda.nh.varnames <- c("outcome_factor", "x1", "x2_factor", "x4", "x3_factor")
      tmp.traindf <- training.df[[s]][,pda.nh.varnames]
      tmp.testdf <- testing.df[[s]][,pda.nh.varnames]
      #catch fit and warnings
      capture <- catchToList(pda(outcome_factor ~ ., data = tmp.traindf, method=gen.ridge))
      tmpPDAfit <- capture$value
      PDANmess <- ifelse(is.null(capture$warnings), "No warnings",capture$warnings)
      PDANW <- ifelse(!is.null(capture$warnings),1,0) #1 if there is a warning, 0 if no warning
      # make predictions
      pda.pred <- predict(tmpPDAfit, tmp.testdf[,-1])
      # results
      confmat.pda.nh[[s]] <- confusionMatrix(pda.pred, tmp.testdf$outcome_factor, positive='Yes')
      
      # 4)
      # PDA, account for hierarchy
      pda.h.varnames <- c("outcome_factor","subject", "x1", "x2_factor", "x4", 
                          "x3_factor")
      tmp.traindf <- training.df[[s]][,pda.h.varnames]
      tmp.testdf <- testing.df[[s]][,pda.h.varnames]
      #catch fit and warnings
      capture <- catchToList(pda(outcome_factor ~ x1 + x2_factor + x3_factor + x4, hier=subject,  
                                 data = tmp.traindf, method=gen.ridge))
      tmpPDAfit2 <- capture$value
      PDAHmess <- ifelse(is.null(capture$warnings), "No warnings",capture$warnings)
      PDAHW <- ifelse(!is.null(capture$warnings),1,0) #1 if there is a warning, 0 if no warning
      # make predictions
      possibleError2 <- tryCatch(
        pda.pred <- predict(tmpPDAfit2, tmp.testdf[,-1]),
        error = function(e) e
      )
      #pda.pred <- predict(tmpPDAfit2, tmp.testdf[,-1])
      if(!inherits(possibleError2, "error")){
        PDAHerr <- paste("No error")
        PDAHE <- 0 #1 if there is an error, 0 if no error
        # results
        confmat.pda.h[[s]] <- confusionMatrix(pda.pred, tmp.testdf$outcome_factor, positive='Yes')
      }else{
        PDAHerr <- possibleError2$message #capture the error message from PDA H
        PDAHE <- 1 #1 if there is an error, 0 if no error
      }
      # results
      #confmat.pda.h[[s]] <- confusionMatrix(pda.pred, tmp.testdf$cancer_factor, positive='Metastatic')
      
      # 5)
      # RF, ignore hierarchy
      rf.nh.varnames <- c("outcome_factor", "x1", "x2_factor", "x4", "x3_factor")
      tmp.rf <- sim.df[[s]][,rf.nh.varnames]
      library(randomForest)
      #catch fit and warnings
      capture <- catchToList(randomForest(outcome_factor ~ x1 + x2_factor + x3_factor + x4, 
                                          data=tmp.rf,ntrees=500,importance=TRUE, x.rep=(3*n)))
      tmpRFfit <- capture$value
      RFNmess <- ifelse(is.null(capture$warnings), "No warnings",capture$warnings)
      RFNW <- ifelse(!is.null(capture$warnings),1,0) #1 if there is a warning, 0 if no warning
      # make predictions
      rf.pred <- predict(tmpRFfit)
      # results
      confmat.rf.nh[[s]] <- confusionMatrix(rf.pred, tmp.rf$outcome_factor, positive='Yes')
      
      # 6)
      # RF, account for hierarchy
      rf.h.varnames <- c("outcome_factor","subject","node", "x1", "x2_factor", "x4", "x3_factor")
      tmp.rf <- sim.df[[s]][,rf.h.varnames]
      #catch fit and warnings
      capture <- catchToList(mult.rf(outcome_factor ~ x1 + x2_factor +  x3_factor + x4, 
                                     data=tmp.rf, 
                                     L1hier="subject",
                                     L2hier="node", x.rep=25))
      tmpRFfit2 <- capture$value
      RFHmess <- ifelse(is.null(capture$warnings), "No warnings",capture$warnings)
      RFHW <- ifelse(!is.null(capture$warnings),1,0) #1 if there is a warning, 0 if no warning
      # collect predictions
      rf.hier.dat <- tmpRFfit2[1]$Final.dat
      rf.hier.dat$maj.vote.pred <- as.factor(rf.hier.dat$maj.vote.pred)
      rf.hier.dat$outcome_factor <- ordered(rf.hier.dat$outcome_factor, levels = c("No","Yes"))
      # results
      confmat.rf.h[[s]] <- confusionMatrix(rf.hier.dat$maj.vote.pred, rf.hier.dat$outcome_factor, 
                                           positive='Yes')
  
    }else{ #if there is an error in GLMM...
      Gerr <- possibleError$message #capture the error message from GLMM
      GE <- 1 #1 if there is an error, 0 if no error
    }

    
    #catch all warnings and errors
    sim <- s
    warn.error = rbind(warn.error, 
                       data.frame(sim=sim, GW=GW, Gmess=Gmess, GE=GE, Gerr=Gerr, 
                                  PDAHE=PDAHE, PDAHerr=PDAHerr,
                                  NBW=NBW, NBmess=NBmess,
                                  PDANW=PDANW, PDANmess=PDANmess,
                                  PDAHW=PDAHW, PDAHmess=PDAHmess,
                                  RFNW=RFNW, RFNmess=RFNmess,
                                  RFHW=RFHW, RFHmess=RFHmess))
  } #end for(s in 1:nsim) simulation loop
  #Prepare the results for the entire simulation, return list of lists
  outList <- list("GLMM.ConfMatrix" = confmat.glmm, 
                  "NB.ConfMatrix.6040" = confmat.6040, 
                  "PDA.ConfMatrix.NoHier" = confmat.pda.nh, 
                  "PDA.ConfMatrix.Hier" = confmat.pda.h, 
                  "RF.ConfMatrix.NoHier" = confmat.rf.nh,
                  "RF.ConfMatrix.Hier" = confmat.rf.h,
                  "Warnings.Errors" = warn.error)
  return(outList)
} #end Simulation Function 


#=================================================================================
# Run a Test of the Simulation
#=================================================================================
# REMINDER:
# nsim is number of simulations
# n is the number of SUBJECTS (i.e. Level1 hierarchy) in training & testing data
# B0, B1, B2, B3, and B4 are the betas for the model
# sd is the between-subject sd for G matrix (bi0)
# kmin is the min amount of lesions a person can have
# kmax is the max amount of lesions a person can have

system.time( sim.test <- hiersimfun(nsim=10, n=80, 
              B0=beta_int, B1=beta_x1, B2=beta_x2, B3=beta_x3, B4=beta_x4,
              sd=sigma_int,kmin=1, kmax=5) )
sim.test$Warnings.Errors

# use results
str(sim.test)
sim.test[[1]] #GLMM, or do...
sim.test$GLMM.ConfMatrix
sim.test[[2]] #NB60/40 or do...
sim.test$NB.ConfMatrix.6040
sim.test$Warnings.Errors
#proportion of bad fit GLMM's (warnings)
(glmm.bad.fits.test <- sum(sim.test$Warnings.Errors$GW)/ length(sim.test$GLMM.ConfMatrix))
sim.test$GLMM.ConfMatrix[[1]]$overall


