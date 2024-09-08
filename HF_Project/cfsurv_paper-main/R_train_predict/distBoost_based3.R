#' Confidence interval based on distributional boosting
#'
#' Construct conformal predictive interval based on distributional boosting
#' @export

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
  print("score:")
  print(score)
  calib_term <- sapply(X=weight_new,get_calibration,score=score,
                         weight_calib=weight_calib,alpha=alpha)
  ##   calib_term <- pmin(calib_term, 1)
  print("calib_term:")
  print(calib_term)

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
  print("median_test:")
  print(median_test)
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
  hist(lower_bnd)
  return(lower_bnd)
}

