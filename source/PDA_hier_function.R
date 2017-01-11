#________________________________________________________________________
# 
#  THESIS: Penalized Discriminant Analysis (PDA)
#
#  Author: Katie Roberts
#  Last edited: 12/6/2016
# 
#
#  Goals: 
#  1. Update functions from the 'mda' package to accomodate 
#     a hierarchically structured data set
#         i) Update the gen.ridge function -> gen.ridge.pda
#        ii) Update the contr.fda function -> contr.pda
#       iii) Update the fda function -> pda
#  2. Test functions ignoring hierarchy
#         i) Test original IRIS data fda vs. pda with no hierarchy
#        ii) Test thyroid data fda vs. pda with no hierarchy
#  3. Create IRIS data
#         i) Manipulate the Iris data set to contain a psuedo 
#             hierarchically structured data set 
#  4. Test functions assuming hierarchy
#         i) Test manip. IRIS data with pda using hierarchy
#        ii) Test thyroid data with pda using hierarchy
#  5. Try Thyroid hierarchy data with these models:
#         i) Univariate
#        ii) Final GLMM model used
#             final GLMM model: cancer~ logvol + echoTextureNonHomo + microcalcYes  + (1|patient_num)
#       
#
# Packages used:
# "mda"
#
# Citation of "mda" package used and modified here:
#    S original by Trevor Hastie & Robert Tibshirani. Original R port by Friedrich
#     Leisch, Kurt Hornik and Brian D. Ripley. (2016). mda: Mixture and Flexible
#     Discriminant Analysis. R package version 0.4-9.
#     https://CRAN.R-project.org/package=mda
#________________________________________________________________________

#=================================================================================
#
# PDA EDITED FUNCTIONS FOR HIERARCHICAL DATA
#
#=================================================================================
#install.packages("mda")
library(mda)

#-----------------------------------------------------------------
# GEN RIDGE FUNCTION
#-----------------------------------------------------------------
#gen.ridge
#Perform a penalized regression, as used in penalized discriminant analysis.
gen.ridge.pda <- function (x, y, weights, lambda = 1, h=FALSE,means, omega, df, ...) 
{
  
  
  
  if (h==FALSE){
    if (missing(df) && lambda <= .Machine$double.eps) #.Machine$double.esp is the smallest postitive floating-point number on the machine running R.
      return(polyreg(x, y)) #do a simple polynomial regression on x and y if missing df and lambda is essentially zero.
    d <- dim(x) #dim of the predictor matrix
    mm <- apply(x, 2, mean) #means by column of the pred matrix
    x <- scale(x, mm, FALSE) #center x about the means but don't scale
    
    simple <- if (missing(omega)) 
      TRUE
    else FALSE
    if (!simple) {
      if (!all(match(c("values", "vectors"), names(omega), 
                     FALSE))) 
        stop("You must supply an  eigen-decomposed version of omega")
      vals <- pmax(sqrt(.Machine$double.eps), sqrt(omega$values))
      basis <- scale(omega$vectors, FALSE, vals)
      x <- x %*% basis
    }
    svd.x <- svd(x)
    dd <- svd.x$d
    if (!missing(df)) 
      lambda = df.gold(dd, df)
    df = sum(dd^2/(dd^2 + lambda))
    y <- (t(t(y) %*% svd.x$u) * dd)/(dd^2 + lambda)
    coef <- svd.x$v %*% y
    fitted <- x %*% coef
    if (!simple) 
      coef <- basis %*% coef
    structure(list(fitted.values = fitted, coefficients = coef, 
                   df = df, lambda = lambda, xmeans = mm), class = "gen.ridge")
    
  }else {
    
    if(!missing(means)){
      mm <- means
    }
    
    if (missing(df) && lambda <= .Machine$double.eps)
      return(polyreg(x, y)) 
    
    simple <- if (missing(omega)) 
      TRUE
    else FALSE
    if (!simple) {
      if (!all(match(c("values", "vectors"), names(omega), 
                     FALSE))) 
        stop("You must supply an  eigen-decomposed version of omega")
      vals <- pmax(sqrt(.Machine$double.eps), sqrt(omega$values))
      basis <- scale(omega$vectors, FALSE, vals)
      x <- x %*% basis
    }
    
    svd.x <- svd(x)
    dd <- svd.x$d
    if (!missing(df)) 
      lambda = df.gold(dd, df)
    df = sum(dd^2/(dd^2 + lambda))
    y <- (t(t(y) %*% svd.x$u) * dd)/(dd^2 + lambda)
    coef <- svd.x$v %*% y
    fitted <- x %*% coef
    if (!simple) 
      coef <- basis %*% coef
    structure(list(fitted.values = fitted, coefficients = coef, 
                   df = df, lambda = lambda, xmeans = mm), class = "gen.ridge")
  }
  
  
}
#gen.ridge(x, y, weights, lambda=1, omega, df, ...)
#gen.ridge.pda()

#-----------------------------------------------------------------
# CONTRAST FUNCTION
#-----------------------------------------------------------------
#contr.fda function
#runs matrix of contrasts to compute QR decomp. matrix
contr.pda <- function (p = rep(1, d[1]), contrast.default = contr.helmert(length(p))) 
{
  d <- dim(contrast.default)
  sqp <- sqrt(p/sum(p))
  x <- cbind(1, contrast.default) * outer(sqp, rep(1, d[2] +  1))
  qx <- qr(x)
  J <- qx$rank
  qr.qy(qx, diag(d[1])[, seq(2, J)])/outer(sqp, rep(1, J - 1)) #computes QR decomp. matrix
}


#-----------------------------------------------------------------
# PDA FUNCTION
#-----------------------------------------------------------------
#PDA
#hier = level 1 hierarchy in data if applicable. For the thyroid data this would be subject
pda <- function (formula = formula(data), data = sys.frame(sys.parent()), hier,
                    weights, theta, dimension = J - 1, eps = .Machine$double.eps, 
                    method = gen.ridge.pda, keep.fitted = (n * dimension < 5000), ...) 
{
  
  this.call <- match.call() #will add argument names if they weren't explicit
  m <- match.call(expand.dots = FALSE) #don't expand the ... arguments  
  m[[1]] <- as.name("model.frame") #identify/label the model.frame that's read in
  m <- m[match(names(m), c("", "formula", "data", "hier", "weights"), 0)] #identify/label parts of the function that are read in
  m <- eval(m, parent.frame()) #evaluate m at the parent.frame environment (default)
  
  Terms <- attr(m, "terms") #get term attributes of m
  g <- model.extract(m, "response") #returns the response component of the model frame m
  x <- model.matrix(Terms, m) #creates a design (model) matrix
  
  
  if (attr(Terms, "intercept")) 
    x = x[, -1, drop = FALSE] #if there's an intercept, drop it from x
  dd <- dim(x)
  n <- dd[1] #number of records
  weights <- model.extract(m, weights)
  if (!length(weights)) 
    weights <- rep(1, n) #if no weights, then create numeric list of 1's of length n
  else if (any(weights < 0)) 
    stop("negative weights not allowed")
  if (length(g) != n) 
    stop("g should have length nrow(x)")
  fg <- factor(g)
  prior <- table(fg) #table of factored response variable
  prior <- prior/sum(prior) #converted to percentage (fraction)
  cnames <- levels(fg) #response variable names
  g <- as.numeric(fg) #converts factored levels to numbers
  J <- length(cnames)  #number of levels for response variable
  iswt <- FALSE
  if (missing(weights)) 
    dp <- table(g)/n
  else {
    weights <- (n * weights)/sum(weights)
    dp <- tapply(weights, g, sum)/n
    iswt <- TRUE
  }
  
  if (missing(theta)) 
    theta <- contr.helmert(J) #runs matrix of contrasts
  theta <- contr.pda(dp, theta) #function that creates contrasts to compute QR decomp. matrix
  
  
  
  
  if (missing(hier)) {
    #continue with original function...
    Theta <- theta[g, , drop = FALSE] #applies the theta matrix contrasts to full n matrix 
    fit <- method(x, Theta, weights, ...) # polyreg fit method with x=predictor matrix and Theta=response matrix
    if (iswt) 
      Theta <- Theta * weights
    ssm <- t(Theta) %*% fitted(fit)/n #transpose Theta multiplied by the fitted values of fit/n
    
  }  else {
    #utilize the hierarchy if applicable
    hier <- model.extract(m, hier)
    N <- length(unique(hier)) #number of unique level 1 hierarchical records
    #Theta placeholder
    dimT <- dim(theta) #collect dimensions of theta
    Tcol <- dimT[2] #collect number of columns from theta
    Theta = matrix(ncol=Tcol) #create a placeholder matrix for Theta
    #my.gen.x placeholder
    # my.gen.x = data.matrix(x[0]) #create empty data matrix to fill
    my.gen.x <- matrix(ncol =  ncol(x))
    colnames(my.gen.x) <- colnames(x)
    #my.gen.mm placeholder to collect xmeans
    my.gen.mm <- matrix(ncol =  ncol(x))
    colnames(my.gen.mm) <- colnames(x)
    
    #create if statement to separate those that only have one factor level vs. those with two or more.
    subj <- data.frame(m$`(hier)`)
    names(subj) <- "subj"
    x2 <- cbind(x, subj)
    # unique(x2$subj)
    
    # BEGIN THE LEVEL 1 FOR-LOOP
    for (i in unique(m$`(hier)`)){ #for each level 1 hierarchy (subject...)
      #execute the pda function for nodes within each subject
      hier.data <- m[m$`(hier)`==i,]
      #
      xi.var <- x2[x2$subj==i,]
      xi <- xi.var[,-which(names(xi.var) == "subj")] 
      
      ddi <- dim(xi)
      ni <- ddi[1]
      gi <- model.extract(hier.data, "response")
      gi <- as.numeric(gi) 
      #
      Thetai <- theta[gi, , drop=FALSE]
      Thetai <- Thetai/ni
      # 
      my.gen.d <- dim(xi)
      my.gen.mmi <- apply(xi, 2, mean)
      my.gen.xi <- scale(xi, my.gen.mmi, FALSE)
      my.gen.xi <- my.gen.xi/ni
      
      
      my.gen.mm = rbind(my.gen.mm,my.gen.mmi) #stack xmeans
      Theta = rbind(Theta, Thetai) #stack Thetai's
      my.gen.x = rbind(my.gen.x, my.gen.xi) #stack the my.gen.x's
      
    } #end level 1 for-loop
    
    
    # 
    #remove first row in Thetai's and my.gen.x
    # Theta <- Theta[-1,]
    # my.gen.x <- my.gen.x[-1,]
    # my.gen.mm <- my.gen.mm[-1,]
    Theta <- t(t(Theta[-1,]))
    my.gen.x <- t(t(my.gen.x[-1,]))
    my.gen.mm <- t(t(my.gen.mm[-1,]))
    mm <- apply(my.gen.mm, 2, mean) #means by column of the pred matrix
    
    
    #now use the method=gen.ridge for hierarchy (h!=1)
    fit <- method(my.gen.x, Theta, h=TRUE,means=mm, weights, ...) # polyreg fit method with my.gen.x=predictor matrix, Theta=response matrix, h=2 for hierarchy
    #structure(list(fitted.values = fitted, coefficients = coef, 
    #               df = df, lambda = lambda), class = "gen.ridge")
    if (iswt) 
      Theta <- Theta * weights
    ssm <- t(Theta) %*% fitted(fit)/N #transpose Theta multiplied by the fitted values of fit/n
    
    
    
    
  }
  #get out: Theta, fit, ssm
  
  
  
  
  ed <- svd(ssm, nu = 0) #singular value decomposition of matrix ssm. nu= number of left singular vectors to be computed
  thetan <- ed$v #matrix whose columns contain the right singular vectors of ssm
  lambda <- ed$d #vector containing the singular values of ssm
  # eps = .Machine$double.eps means the smallest positive floating-point number x such that 1+x != 1. Normally 2.220446e-16
  #dimension = J - 1 number of response factors minus 1
  lambda[lambda > 1 - eps] <- 1 - eps #convert value of lambda that are essentially greater than 1 to 1 minus essentially zero =~1
  discr.eigen <- lambda/(1 - lambda)
  pe <- (100 * cumsum(discr.eigen))/sum(discr.eigen)
  dimension <- min(dimension, sum(lambda > eps))
  if (dimension == 0) {
    warning("degenerate problem; no discrimination")
    return(structure(list(dimension = 0, fit = fit, call = this.call), 
                     class = "fda"))
  }
  thetan <- thetan[, seq(dimension), drop = FALSE]
  pe <- pe[seq(dimension)]
  alpha <- sqrt(lambda[seq(dimension)])
  sqima <- sqrt(1 - lambda[seq(dimension)])
  vnames <- paste("v", seq(dimension), sep = "")
  means <- scale(theta %*% thetan, FALSE, sqima/alpha) #scale theta%*%thetan by sqima/alpha
  dimnames(means) <- list(cnames, vnames)
  names(lambda) <- c(vnames, rep("", length(lambda) - dimension))
  names(pe) <- vnames
  
  obj <- structure(list(percent.explained = pe, values = lambda, 
                        means = means, theta.mod = thetan, dimension = dimension, 
                        prior = prior, fit = fit, call = this.call, terms = Terms), 
                   class = "fda")
  obj$confusion <- confusion(predict(obj), fg)
  if (!keep.fitted) 
    obj$fit$fitted.values <- NULL
  obj
}






#=================================================================================
#
# TEST THE FUNCTIONS WITH DATA BELOW
#
#=================================================================================

#-----------------------------------------------------------------
# TEST THE FUNCTION WITH NAIVE APPROACHES COMPARED TO KNOWN MDA
#-----------------------------------------------------------------
# IRIS - they match
data(iris)
head(iris)

#known function with iris data
set.seed(2345)
irisfit <- fda(Species ~ ., data = iris, method=gen.ridge)
irisfit
confusion(irisfit, iris)
irisfit$confusion
plot(irisfit,coords=c(1,2))
coef(irisfit)
# str(irisfit)

#pda function with iris data
set.seed(2345)
irisfit1 <- pda(Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width, data = iris, method=gen.ridge)
irisfit1
confusion(irisfit1, iris)
irisfit1$confusion
plot(irisfit1)
coef(irisfit1)
# str(irisfit1)


#-----------------------------------------------------------------
# CREATE TESTING IRIS DATA TO USE WITH HIERARCHY
#-----------------------------------------------------------------
#create iris data with added hierarchy
# NOTE: this hierarchy is generated to test the functionality of the code ONLY. There is no correlation
# here so the results should be less than interesting
data(iris)
my.iris <- iris
table(my.iris$Species)
my.iris$subject <- c(rep(1:50,each=3))
my.iris <- my.iris[with(my.iris, order(subject)), ]
#my.iris$node <- ave(my.iris$subject, my.iris$subject, FUN = seq_along)
#my.iris$indic <- sample(0:1, 150, replace=T) # random indicator variable
my.iris
set.seed(165)
my.iris <- my.iris[sample(nrow(my.iris), 100), ]
#sort data by subject and node
my.iris <- my.iris[order(my.iris$subject),]
my.iris$node <- ave(my.iris$subject, my.iris$subject, FUN = seq_along)
my.iris
#need to randomize the order of subject and node so that the outcome (species) will be spread over subjects
subnode <- my.iris[6:7]
subnode
names(subnode) <- c("sub1","node1")
set.seed(165)
ran <- subnode[sample(nrow(subnode)),]
my.iris <- cbind(my.iris, ran)
my.iris <- my.iris[,-c(6:7)]
names(my.iris)[names(my.iris)=="sub1"] <- "subject"
names(my.iris)[names(my.iris)=="node1"] <- "node"
my.iris <- my.iris[with(my.iris, order(subject)), ]
my.iris


#-----------------------------------------------------------------
# TEST THE FUNCTION WITH HIERARCHICALLY STRUCTURED DATA
#-----------------------------------------------------------------
# IRIS
head(my.iris)

#my function with hierarchical my.iris data
set.seed(2345)
irisfit2 <- pda(Species ~ Sepal.Length+Sepal.Width+Petal.Length+Petal.Width, hier= subject,data = my.iris, method=gen.ridge)
irisfit2
#confusion(irisfit2, my.iris)
irisfit2$confusion
plot(irisfit2)
coef(irisfit2)





# Compare non-hier to hier IRIS data
head(my.iris)

#my function with iris data
set.seed(2345)
irisfit.NOhier <- pda(Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width, data = my.iris, method=gen.ridge)
irisfit.NOhier
confusion(irisfit.NOhier, my.iris)
irisfit.NOhier$confusion
plot(irisfit.NOhier)
coef(irisfit.NOhier)
# str(irisfit.NOhier)

#my function with hierarchical my.iris data
set.seed(2345)
irisfit.hier <- pda(Species ~ Sepal.Length+Sepal.Width+Petal.Length+Petal.Width, hier= subject,data = my.iris, method=gen.ridge)
irisfit.hier
#confusion(irisfit.hier, my.iris)
irisfit.hier$confusion
plot(irisfit.hier)
coef(irisfit.hier)



