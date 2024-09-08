#' Fitting the censoring probability P(C>=c|X)
#'
#' @export



censoring_prob_train2 <- function(fit, calib,
                           xnames, c,
                           ftol=.1, tol=.1, n.tree = 40){

  p <- length(xnames)
  ## Fitting P(-C<=-c_0|X) (since P(C>=c_0|X)=P(-C<=-c_0|X))
  fit$C <- -fit$C
  fmla <- with(fit,as.formula(paste("C ~ ", paste(xnames, collapse= "+"))))
  gbm_mdl <- gbm(fmla,data=fit,distribution="gaussian", n.tree = n.tree)
  median_fit<- predict(object=gbm_mdl,newdata = fit)
  res_fit <- fit$C-median_fit
  resamp_fit <- median_fit + res_fit[sample.int(dim(fit)[1])]
  if(p==1){
    xdf <- data.frame(X1=fit[,colnames(fit)%in%xnames],
                  X2=rep(1,dim(fit)[1]))
    }else{
    xdf <- fit[,colnames(fit)%in%xnames]
    }
  mdlrb <- modtrast(xdf,fit$C,resamp_fit,min.node=200)
  
  ## Computing the censoring scores for the fitting data
  pr_fit <- rep(NA, dim(fit)[1])
  for(i in 1:length(pr_fit)){
      pr_fit[i] <- distBoost_cdf(mdlrb,xdf[i,],median_fit[i],-c,res_fit)
  }

  ## Computing the censoring scores for the calibration data
  pr_calib <- rep(NA, dim(calib)[1])
  median_calib<- predict(object=gbm_mdl,newdata = calib)
  if(p==1){
    xdf <- data.frame(X1=calib[,colnames(calib)%in%xnames],
                  X2=rep(1,dim(calib)[1]))
    }else{
    xdf <- calib[,colnames(calib)%in%xnames]
    }

  for(i in 1:length(pr_calib)){
    pr_calib[i] <- distBoost_cdf(mdlrb,xdf[i,],median_calib[i],
                              -c,res_fit)
  }
  ret = list(pr_fit = pr_fit, pr_calib = pr_calib, gbm_mdl = gbm_mdl, mdlrb = mdlrb, res_fit = res_fit)
  return(ret)
}


censoring_prob_test2 <- function(test = NULL,
                           xnames, c, gbm_mdl, mdlrb, res_fit){

  p <- length(xnames)
  ## Computing the censoring scores for the test data
  if(!is.null(test)){
    newdata <- data.frame(test)
    colnames(newdata) <- xnames
    median_test<- predict(object=gbm_mdl,newdata = newdata)
    if(p==1){
      xdf <- data.frame(X1=test,
                  X2=rep(1,length(test)))
      n_new <- length(test)
    }else{
      xdf <- test
      n_new <- dim(test)[1]
    }
    pr_new <- rep(NA, n_new)
    for(i in 1:n_new){
        pr_new[i] <- distBoost_cdf(mdlrb,xdf[i,],median_test[i],-c,res_fit)
    }
  }else{pr_new=NULL}
  return(pr_new = pr_new)
}

##################################################################

censoring_prob_train1 <- function(fit, calib,
                           xnames, c,
                           ftol=.1, tol=.1, n.tree = 40){

  p <- length(xnames)
  X_train <- fit[, colnames(fit) %in% xnames]
  C_train <- as.factor(fit$C > c)
  #print(summary(X_train))
  X_train <- model.matrix(~ 0 + ., X_train)
  print(dim(X_train))
  rf_model <- probability_forest(X_train, C_train)
  pr_fit <- predict(rf_model)$predictions[, "TRUE"]

  X_calib = calib[, colnames(calib) %in% xnames]
  X_calib = model.matrix(~ 0 + ., X_calib)
  pr_calib = predict(rf_model, X_calib)$predictions[, "TRUE"]

  ret = list(pr_fit = pr_fit, pr_calib = pr_calib, rf_model = rf_model)
  return(ret)
}


censoring_prob_test1 <- function(test = NULL,
                           xnames, c, rf_model){

  p <- length(xnames)
  ## Computing the censoring scores for the test data
  if(!is.null(test)){
    #print(summary(test))
    #print(xnames)
    #newdata <- as.matrix(test)
    colnames(test) <- xnames
    #print(summary(newdata))
    #print(summary(rf_model))
    test = model.matrix(~ 0 + ., test)
    print(dim(test))
    pr_new <- predict(rf_model, test)$predictions[, "TRUE"]
  }else{pr_new=NULL}
  return(pr_new = pr_new)
}

##################################################################
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