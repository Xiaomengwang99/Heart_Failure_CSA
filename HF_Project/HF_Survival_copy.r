suppressPackageStartupMessages(library(conTree))
suppressPackageStartupMessages(library(grf))
source("cfsurv_paper-main/utils/source.R")

#seed = Sys.getenv('SGE_TASK_ID')
#seed = as.integer(seed)
seed = 123
set.seed(seed)


data = read.csv("final_df_May10.csv")
head(data)
colnames(data)


# Convert columns after censoring_time into true/false
data[, -(1:20)] <- lapply(data[, -(1:20)], as.factor)

other_variables = colnames(data)[c(21:33, 46:49)][-11]
all_variables = c("insurance", "language", "marital_status", "los", "gender", "anchor_age", "race", other_variables)
all_variables

data$censored_T = ifelse(is.na(data$dod_since_disch), data$censoring_time , data$dod_since_disch)
data = data[data$censored_T > 0,]
n = dim(data)[1]

X = data[, all_variables]
X[, "insurance"] = as.factor(X[, "insurance"])
X[, "language"] = as.factor(X[, "language"])
X[, "marital_status"] = as.factor(X[, "marital_status"])
X[, "gender"] = as.factor(X[, "gender"])
X[, "race"] = as.factor(X[, "race"])
censored_T = data$censored_T
C = data$censoring_time
event = censored_T < C
T = ifelse(is.na(data$dod_since_disch), Inf, data$dod_since_disch)


df = data.frame(X = X, C = C, censored_T = censored_T, event = event)

set.seed(123)


train_idx = sample(1:n, n/2)
predicted_lbd = rep(0, n)

for (i in 1:1){
    df_train = df[train_idx, ]
    #saveRDS(df_train[, 1:p], "df_train.rds")
    df_test = df[-train_idx, ]
    p = dim(X)[2]
    c = median(C)

    res <- cfsurv_train(p = p, c_list = c, Xtrain = df_train[, 1:p], C = df_train$C,
                        event = df_train$event, time = df_train$censored_T,
                        alpha = 0.1, model="distBoost_train", seed = seed)

    output <- cfsurv_predict(x = df_test[, 1:p], c = res$c,
                    alpha=0.1,
                    type="quantile",
                    seed = seed,
                    model = "distBoost_predict",
                    dist= "weibull",
                    I_fit = NULL,
                    ftol=.1,tol=.1,
                    n.tree=100, 
                    weight_calib = res$weight_calib, 
                    rf_model = res$rf_model,
                    DB_gbm_mdl = res$DB_gbm_mdl, 
                    DB_mdlrb = res$DB_mdlrb, 
                    score = res$score, 
                    res_T = res$res_T
                    )

    predicted_lbd[-train_idx] = output
}

for (i in 1:1){
    df_train = df[-train_idx, ]
    df_test = df[train_idx, ]
    p = dim(X)[2]
    c = median(C)

    res <- cfsurv_train(p = p, c_list = c, Xtrain = df_train[, 1:p], C = df_train$C,
                        event = df_train$event, time = df_train$censored_T,
                        alpha = 0.1, model="distBoost_train", seed = seed)

    output <- cfsurv_predict(x = df_test[, 1:p], c = res$c,
                    alpha=0.1,
                    type="quantile",
                    seed = seed,
                    model = "distBoost_predict",
                    dist= "weibull",
                    I_fit = NULL,
                    ftol=.1,tol=.1,
                    n.tree=100, 
                    weight_calib = res$weight_calib, 
                    rf_model = res$rf_model,
                    DB_gbm_mdl = res$DB_gbm_mdl, 
                    DB_mdlrb = res$DB_mdlrb, 
                    score = res$score, 
                    res_T = res$res_T
                    )

    predicted_lbd[train_idx] = output
}

df$predicted_lbd = predicted_lbd
head(df)

write.csv(df, "df_predicted_lbd.csv", row.names = FALSE)

threshold = 364
count = sum(df$predicted_lbd > threshold)
count

