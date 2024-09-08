#' Predictive confidence interval for survival data
#'
#' The main function to generate a predictive conformal confidence interval for a unit's survival time.
#'
#' @param x a vector of the covariate for test point. 
#' @param c the censoring time for the test point.
#' @param Xtrain a n-by-p matrix of the covariate of the training data.
#' @param C a length n vector of the censoring time of the training data.
#' @param event a length n vector of indicators if the time observed is censored. TRUE corresponds to NOT censored, and FALSE censored.
#' @param time  a vevtor of length n, containing the observed survival time.
#' @param alpha a number between 0 and 1, speciifying the miscoverage rate.
#' @param seed an integer random seed (default: 24601).
#' @param model Options include "cox", "randomforest", "Powell", "Portnoy" and "PengHuang". This determines the model used to fit the condditional quantile (default: "cox").
#' @param dist either "weibull", "exponential" or "gaussian" (default: "weibull"). The distribution of T used in the cox model. 
#' @param h the bandwidth for the local confidence interval. Default is 1.
#'
#' @return low_ci a value of the lower bound for the survival time of the test point.
#' @return includeR 0 or 1, indicating if [r,inf) is included in the confidence interval.
#'
#' @examples
#' # Generate data
#' n <- 500
#' X <- runif(n,0,2)
#' T <- exp(X+rnorm(n,0,1))
#' R <- rexp(n,rate = 0.01)
#' event <- T<=R
#' time <- pmin(T,R)
#' data <- data.frame(X=X,R=R,event=event,censored_T=censored_T)
#' # Prediction point
#' x <- seq(0,2,by=.4)
#' r <- 2
#' # Run cfsurv
#' res <- cfsurv(x,r,X,R,event,time,alpha=0.1,model="cox")
#'
#' @export

# function to construct conformal confidence interval
cfsurv_train <- function(p, c_list=NULL,
                   pr_list=NULL,
                   pr_new_list=NULL,
                   Xtrain,C,event,time,
                   alpha=0.05,
                   type="quantile",
                   seed = 24601,
                   model = "distBoost_train",
                   dist= "weibull",
                   I_fit = NULL,
                   ftol=.1,tol=.1,
                   n.tree=100
                   ){
  ## Check if the required packages are installed
  ## Solution found from https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
  #list.of.packages <- c("ggplot2",
  #                      "quantreg",
  #                      "grf",
  #                      "quantregForest",
  #                      "randomForestSRC",
  #                      "survival",
  #                      "tidyverse",
  #                      "fishmethods",
  #                      "foreach",
  #                      "doParallel",
  #                      "GauPro",
  #                      "gbm",
  #                      "np")
  #list.of.packages <- c("ggplot2",
  #                      "survival",
  #                      "gbm")
  #new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  #if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')
  #suppressPackageStartupMessages(res <- lapply(X=list.of.packages,FUN=require,character.only=TRUE))
  ## Process the input
  ## Check the length of x and c: only two cases are supported. length(r)=1, or length(r)=length(x)
  X <- Xtrain
  if(is.null(dim(X)[1])){
    n <- length(X)
    pX <- 1
  }else{
    n <- dim(X)[1]
    pX <- dim(X)[2]
  }
  

  ## Check the type of the model. Only "cox" and "randomforest" are supported
  #if(model %in% c("cox","randomforest","pow","portnoy","PengHuang",
  #                "distBoost_train","gpr", "quantBoost")==0) 
  #  stop("The regression model is not supported.")

  ## Check the type of the confidence inteval
  if(type %in% c("quantile","percentile")==0) stop("The type of confidence interval is not supported.")

  ## Check the value of alpha
  if (alpha>=1 | alpha<=0) stop("The value of alpha is out of bound.")

  ## Check the dimensions of the data 
  xnames <- paste0('X', 1:p)
  if(n != length(C))stop("The number of rows in X does not match the length of R.")
  if(length(C) != length(event))stop("The length of R does not match the length of event.")
  if(length(event) != length(time))stop("The length of event does not match the length of time.")
  if(p != pX) stop("The dimension of the test point does not match the dimension of the training point.")

  data <- as.data.frame(cbind(C,event,time,X))
  colnames(data) <- c("C","event","censored_T",xnames)

  ## set random seed
  set.seed(seed)

  ## Split the data into the training set and the calibration set
  n = dim(data)[1]
  n_train = n/2
  n_calib = n-n_train
  if(is.null(I_fit)){
    I_fit <- sample(1:n,n_train,replace = FALSE)
  }
  data_fit <- data[I_fit,]
  data_calib <- data[-I_fit,]
  
  ## If c is not specified, select c automatically 
  if(is.null(c_list)){
    ref_length <- 10
    c_list <- seq(min(data_fit$C),max(data_fit$C),length=ref_length)
  }
  
  if(length(c_list)==1){
    c <- c_list
    if(is.null(pr_list) | is.null(pr_new_list)){
      res <- censoring_prob_train(data_fit,data_calib, xnames = xnames,c = c ,ftol = ftol,tol = tol)
      pr_calib <- res$pr_calib
      rf_model <- res$rf_model
    }else{
      pr_calib <- pr_list[-I_fit]
    }
  }else{
    if(is.null(pr_list) | is.null(pr_new_list)){
      res <- selection_c(X=data_fit[,colnames(data_fit)%in%xnames],
                         C=data_fit$C,
                         event=data_fit$event,
                         time=data_fit$censored_T,
                         weight_ref=NULL,
                         alpha,c_ref=c_list,
                         type=type,dist=dist)
      c <- res$c_opt
      res <- censoring_prob_train(data_fit,data_calib, xnames = xnames,c = c ,ftol = ftol,tol = tol)
      pr_calib <- res$pr_calib
      rf_model <- res$rf_model
    }else{
      weight_ref <- 1/pr_list[I_fit,]
      res <- selection_c(X=data_fit[,colnames(data_fit)%in%xnames],
                         C=data_fit$C,
                         event=data_fit$event,
                         time=data_fit$censored_T,
                         alpha,c_ref=c_list,
                         weight_ref=weight_ref,
                         type=type,dist=dist)
      c <- res$c_opt
      pr_calib <- pr_list[-I_fit,c_list==c] 
    }
  }
  ## Computing the weight for the calibration data and the test data
  weight_calib <- 1/pr_calib
 
  ## Run the main function and gather resutls
  if(model == "distBoost_train"){
    res = distBoost_based_train(p,c,alpha,
                    data_fit,
                    data_calib,
                    weight_calib,
                    n.tree)
    DB_gbm_mdl = res$gbm_mdl
    res_T = res$res_T
    DB_mdlrb = res$mdlrb
    score = res$score
    weight_calib = res$weight_calib
    ret_lst = list(weight_calib = weight_calib, rf_model = rf_model, DB_gbm_mdl = DB_gbm_mdl, 
                 DB_mdlrb = DB_mdlrb, score = score, res_T = res_T, c = c)
    return(ret_lst)
   }
  return("Error: model not found")

}

cfsurv_predict <- function(x,c,
                   pr_list=NULL,
                   pr_new_list=NULL,
                   alpha=0.05,
                   type="quantile",
                   seed = 24601,
                   model = "distBoost_predict",
                   dist= "weibull",
                   I_fit = NULL,
                   ftol=.1,tol=.1,
                   n.tree=100, weight_calib, rf_model, 
                   DB_gbm_mdl= NULL, DB_mdlrb= NULL, score, res_T= NULL,
                   mdl = NULL, bw = NULL
                   ){
  ## Check if the required packages are installed
  ## Solution found from https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
  #list.of.packages <- c("ggplot2",
  #                      "quantreg",
  #                      "grf",
  #                      "quantregForest",
  #                      "randomForestSRC",
  #                      "survival",
  #                      "tidyverse",
  #                      "fishmethods",
  #                      "foreach",
  #                      "doParallel",
  #                      "GauPro",
  #                      "gbm",
  #                      "np")
  #list.of.packages <- c("ggplot2",
  #                      "survival",
  #                      "gbm")
  #new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  #if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')
  #suppressPackageStartupMessages(res <- lapply(X=list.of.packages,FUN=require,character.only=TRUE))
  ## Process the input
  ## Check the length of x and c: only two cases are supported. length(r)=1, or length(r)=length(x)
  if(is.null(dim(x)[1])){
    len_x <- length(x)
    p <- 1
  }else{
    len_x <- dim(x)[1]
    p <- dim(x)[2]
  }


  ## Check the type of the model. Only "cox" and "randomforest" are supported
  #if(model %in% c("cox","randomforest","pow","portnoy","PengHuang",
  #                "distBoost_predict","gpr", "quantBoost")==0) 
  #  stop("The regression model is not supported.")

  ## Check the type of the confidence inteval
  if(type %in% c("quantile","percentile")==0) stop("The type of confidence interval is not supported.")

  ## Check the value of alpha
  if (alpha>=1 | alpha<=0) stop("The value of alpha is out of bound.")

  ## Check the dimensions of the data 
  xnames <- paste0('X', 1:p)
  #if(n != length(C))stop("The number of rows in X does not match the length of R.")
  #if(length(C) != length(event))stop("The length of R does not match the length of event.")
  #if(length(event) != length(time))stop("The length of event does not match the length of time.")
  #if(p != pX) stop("The dimension of the test point does not match the dimension of the training point.")

  #data <- as.data.frame(cbind(C,event,time,X))
  #colnames(data) <- c("C","event","censored_T",xnames)

  ## set random seed
  set.seed(seed)

  #newdata <- data.frame(x)
  #colnames(newdata) <- xnames
  

  ## Computing the weight for the calibration data and the test data
  pr_new = censoring_prob_test(test = x, xnames, c, rf_model)
  weight_new <- 1/pr_new
 
  ## Run the main function and gather results

  if(model == "distBoost_predict"){
    res = distBoost_based_predict(x,
                      DB_gbm_mdl, 
                      res_T, 
                      DB_mdlrb, weight_new, weight_calib, alpha, score, c)
   }
  return(res)
}

distBoost_based_train <- function(p,c,alpha,
                      data_fit,
                      data_calib,
                      weight_calib,
                      n.tree=100){
  ## Keep only the data points with C>=c
  ## Transform min(T,C) to min(T,c) 
  weight_calib <- weight_calib[data_calib$C>=c]
  data_calib <- data_calib[data_calib$C>=c,]
  data_calib$censored_T <- pmin(data_calib$censored_T,c)
  
  ## Fit the model for S(y)=p(min(T,c)>=y|X)
  xnames <- paste0("X",1:p)
  data_fit <- data_fit[data_fit$C>=c,]
  data_fit$censored_T <- pmin(data_fit$censored_T,c)

  ## Moving on, use surv_data_fit  
  surv_data_fit <- data_fit
  surv_data_fit$censored_T <- -surv_data_fit$censored_T
  fmla <- with(surv_data_fit,as.formula(paste("censored_T~ ",
                                              paste(xnames, collapse= "+"))))

  gbm_mdl <- gbm(fmla,data=surv_data_fit,distribution="gaussian",
                 n.trees=n.tree)
  capture.output(median_fit<- predict(object=gbm_mdl,newdata = surv_data_fit),
                 file=NULL)
  res_T <- surv_data_fit$censored_T-median_fit
  resamp_T <- median_fit + res_T[sample.int(dim(surv_data_fit)[1])]
  if(p==1){
    xdf <- data.frame(X1=surv_data_fit[,colnames(surv_data_fit)%in%xnames],
                    X2=rep(1,dim(data_fit)[1]))
  }else{
    xdf <- surv_data_fit[,colnames(surv_data_fit)%in%xnames]
  }
  mdlrb <- modtrast(xdf, surv_data_fit$censored_T,
                    resamp_T, min.node = 200)

  ## obtain the score of the calibration data
  surv_data_calib <- data_calib
  surv_data_calib$censored_T <- -surv_data_calib$censored_T
  capture.output(median_calib <- predict(object=gbm_mdl,
                                         newdata = surv_data_calib,
                                         n.trees=n.tree),
                 file=NULL)
  if(p==1){
    xdf <- data.frame(X1=surv_data_calib[,colnames(surv_data_calib)%in%xnames],
                      X2=rep(1,dim(data_calib)[1]))
  }else{
    xdf <- surv_data_calib[,colnames(surv_data_calib)%in%xnames]
  }

  score <- rep(NA,dim(data_calib)[1])
  for(i in 1:length(score)){
    score[i] <- distBoost_cdf(mdlrb,xdf[i,],median_calib[i],
                               surv_data_calib$censored_T[i],
                               res_T)
  }
  
  ret = list(gbm_mdl = gbm_mdl, res_T = res_T, mdlrb = mdlrb, score = score, weight_calib = weight_calib)
  return(ret)
}

distBoost_based_predict <- function(x,
                      gbm_mdl, 
                      res_T, 
                      mdlrb, weight_new, weight_calib, alpha, score, c){
  ## Obtain the final confidence interval
  ## Obtain the calibration term
  #weight_calib = rep(1,length(weight_calib))
  #weight_new = rep(1,length(weight_new))
  #print("weight_new:")
  #print(weight_new)
  #print("score:")
  #print(score)
  calib_term <- sapply(X=weight_new,get_calibration,score=score,
                         weight_calib=weight_calib,alpha=alpha)
  ##   calib_term <- pmin(calib_term, 1)
  #print("calib_term:")
  #print(calib_term)

  if(is.null(dim(x)[1])){
    len_x <- length(x)
    p <- 1
  }else{
    len_x <- dim(x)[1]
    p <- dim(x)[2]
  }
  xnames <- paste0("X",1:p)

  lower_bnd <- rep(0,len_x)
  newdata <- data.frame(x)
  colnames(newdata) <- xnames
  if(p==1){newdata$X2 <- rep(1,len_x)}
  median_test <- predict(object=gbm_mdl,newdata = newdata)
  for(i in 1:len_x){
    if(calib_term[i] == Inf){
      lower_bnd[i] <- 0
    }else{
      qres_fit <- as.numeric(quantile(res_T,calib_term[i]))
      bound <- ydist(mdlrb,newdata[i,],median_test[i]+qres_fit)
      lower_bnd[i] <- -bound
    }
  }
  lower_bnd <- pmax(lower_bnd,0)
  lower_bnd <- pmin(lower_bnd,c)
  return(lower_bnd)
}


censoring_prob_train <- function(fit, calib,
                           xnames, c,
                           ftol=.1, tol=.1, n.tree = 40){
  pr_fit = rep(1, dim(fit)[1])
  pr_calib = rep(1, dim(calib)[1])
  ret = list(pr_fit = pr_fit, pr_calib = pr_calib, rf_model = NULL)
  return(ret)
}


censoring_prob_test <- function(test = NULL,
                           xnames, c, rf_model){


  if(!is.null(test)){
    pr_new <- rep(1, dim(test)[1])
  }else{pr_new=NULL}
  return(pr_new = pr_new)
}

get_calibration <- function(score,weight_calib,weight_new,alpha){
  ## Check input format
  if(length(score)!=length(weight_calib)) stop("The length of score is not compatible with the length of weight!")

  if(!is.numeric(alpha)) stop("alpha should be a real number between 0 and 1!")
  if(alpha>1 | alpha<0) stop("alpha should be a real number between 0 and 1!")

  ## Computing the calibration term
  weight <- c(weight_calib,weight_new)
  weight <- weight/sum(weight)
  score_vec <- c(score,Inf)
  sort_score <- sort(score_vec)
  order_score <- order(score_vec)
  sort_w <- weight[order_score]
  idxw <- min(which(cumsum(sort_w)>=1-alpha))
  calib_term <- sort_score[idxw]

  return(calib_term)
}

distBoost_cdf <- function(mdl, x, z, t, resid){
 perc <- (0 : 1000) / 1000
 qres <- as.numeric(quantile(resid, perc))
 yhat <- ydist(mdl, x, z + qres)
 quant_ind <- suppressWarnings(min(which(yhat >= t)))
 if(quant_ind == Inf){
  perc_est <- 1
 }else{
  perc_est <- perc[quant_ind]
 }
 return(perc_est)
}

selection_c <- function(X,C,event,time,alpha,
                        c_ref,weight_ref,
                        model="cox",
                        type="quantile",
                        dist="weibull"){
  
  ## Get the dimension of the input
  if(is.null(dim(X))){
    n <- length(X)
    p <- 1
  }else{
    n <- dim(X)[1]
    p <- dim(X)[2]
  }
  xnames <- paste0("X",1:p)
  data <- cbind(X,C,event,time)
  data <- data.frame(data)
  colnames(data) <- c(xnames,"C","event","censored_T")

  ## Evaluate the average bound for each candidate c
  bnd_ref <- c()
  for(i in 1:length(c_ref)){
  bnd <- evaluate_length(c_ref[i],alpha=alpha,n=n,p=p,model,
                        data=data,weight=weight_ref[,i],xnames=xnames,
                        type=type,dist=dist)
  bnd_ref <- c(bnd_ref,bnd)
  }
  c_opt <- c_ref[which.max(bnd_ref)]
  return(list(c_opt=c_opt,c_ref=c_ref,bnd_ref=bnd_ref))

}

evaluate_length <- function(c,alpha,n,p,
                            model,
                            data,
                            weight,
                            xnames,
                            type = "quantile",
                            dist = "weibull",
                            seed = 2020){
  ## Determine the fitting set, the calibration set and the test set
  set.seed(seed)
  I_fit <- sample(1:n,floor(n/2),replace=FALSE)
  I_calib <- sample((1:n)[-I_fit],floor(n/4),replace=FALSE)
  I_test <- (1:n)[-c(I_fit,I_calib)]
  
  data_fit <- data[I_fit,]
  data_calib <- data[I_calib,]
  data_test <- data[I_test,]
  
  if(is.null(weight)){
    res <- censoring_prob(fit=data_fit,calib=data_calib,test=data_test,
                          xnames=xnames,c)
    pr_calib <- res$pr_calib
    pr_new <- res$pr_new
    weight_calib <- 1/pr_calib
    weight_new <- 1/pr_new
  }else{
    weight_calib <- weight[I_calib]
    weight_new <- weight[I_test]
  }
  x <- data_test[,colnames(data_test)%in%xnames]
  
  if(model == "cox"){
    bnd <- cox_based(x,c,alpha,
                     data_fit,
                     data_calib,
                     type = "quantile",
                     dist,
                     weight_calib,
                     weight_new)
   }
  
  if(model == "randomforest"){
    bnd <- rf_based(x,c,alpha,
                    data_fit,
                    data_calib,
                    weight_calib,
                    weight_new)
  }
  
  if(model == "pow"){
    bnd <- pow_based(x,c,alpha,
                    data_fit,
                    data_calib,
                    weight_calib,
                    weight_new)
  }

  if(model == "portnoy"){
    bnd <- portnoy_based(x,c,alpha,
                        data_fit,
                        data_calib,
                        weight_calib,
                        weight_new)
  }

  if(model == "PengHuang"){
    bnd <- ph_based(x,c,alpha,
                   data_fit,
                   data_calib,
                   weight_calib,
                   weight_new)
  }

  return(mean(bnd))
}
