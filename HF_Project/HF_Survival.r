suppressPackageStartupMessages(library(conTree))
suppressPackageStartupMessages(library(grf))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(quantreg))
source("HF_Project/cfsurv_paper-main/utils/source.R")

seed = Sys.getenv('SGE_TASK_ID')
seed = as.integer(seed)
#seed = 123
set.seed(seed)



data = read.csv("HF_Project/final_df_May10.csv")
head(data)
colnames(data)


# Convert columns after censoring_time into true/false
data[, -(1:20)] <- lapply(data[, -(1:20)], as.factor)

other_variables = colnames(data)[c(21:33, 46:49)][-11]
all_variables = c("insurance", "language", "marital_status", "los", "gender", "anchor_age", "race", other_variables)
all_variables

data$censored_T = ifelse((is.na(data$dod_since_disch)| (data$dod_since_disch + data$dischtime_since_admission) > 365), data$censoring_time , data$dod_since_disch + data$dischtime_since_admission)
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
#C = 365 - 2*X[, "los"] - X[, "anchor_age"] - runif(n, 0, 50)
#print(summary(C))
#C = runif(n, 100, 365)
censored_T = pmin(data$censored_T, C)
event = censored_T < C
T = ifelse(is.na(data$dod_since_disch), Inf, data$dod_since_disch + data$dischtime_since_admission)
mean(event)


df = data.frame(X = X, C = C, censored_T = censored_T, event = event)
df1 = data.frame(X = X, censored_T = censored_T, event = event)
p = dim(X)[2]
xnames = paste0("X", 1:p)
colnames(df) <- c(xnames, "C", "censored_T", "event")

cor(ifelse(df1$X.has_hypertension == "True", 1, 0), ifelse(df1$X.ACE.Inhibitor == "True", 1, 0))
# Check for missing data in df
#missing_data <- colSums(is.na(df))
#missing_data
#head(df1)
#mean(event)


num_trials = 1

predicted_lbd = matrix(0, nrow = n, ncol = num_trials)
predicted_lbd_cox = matrix(0, nrow = n, ncol = num_trials)

coverage = matrix(0, nrow = num_trials, ncol = 7)

n_train = 8000
n_test = n - n_train

alpha = 0.1

for (i in 1:num_trials){
    print(paste0("Iteration:", i))
    #set.seed(i+12345)

    train_idx = sample(1:n, n_train)
    df_train = df[train_idx, ]
    #saveRDS(df_train[, 1:p], "df_train.rds")
    df_test = df[-train_idx, ]
    
    column_types <- sapply(df_test, class)

    c = median(C)

    #hist(C)
    #abline(v = 365, col = "red")

    cat("Compute the result of Cox...")
    output_cox = rep(NA, n_test)
    fmla <- as.formula(paste("Surv(censored_T, event) ~ ",
                paste(xnames, collapse= "+")))
    mdl <- coxph(fmla, data = as.data.frame(df_train))
    for(j in 1:n_test){
        if (j %% 100 == 0){
            print(j)
        }
        res <- summary(survfit(mdl, newdata = df_test[j, ]))
        time_point <- res$time
        survcdf <- 1 - res$surv
        quant <- time_point[min(which(survcdf >= alpha))]
        output_cox[j] <- quant
    }
    output <- data.frame(cox.bnd = output_cox)
    cat("Done.\n")

    #output_cox2 = cfsurv(x = df_test[, 1:p], c_list = c, X = df_train[, 1:p], 
    #                    C = df_train$C, event = df_train$event, 
    #                    time = df_train$censored_T,
    #                    alpha= alpha,
    #                    model="cox", seed = seed)

    ## AFT + Weibull
    cat("Compute the result of AFT...")
    fmla <- as.formula(paste("Surv(censored_T, event) ~ ",
                paste(xnames, collapse= "+")))
    mdl <- survreg(fmla, data = df_train, dist = "weibull")
    res <- predict(mdl,
                newdata = df_test,
                type = "quantile",
                p = alpha)
    output$aft.bnd <- as.vector(res)
    cat("Done.\n")

    ## Powell
    cat("Compute the result of powell...")
    fmla <- as.formula(paste("Curv(censored_T, C, ctype = \"right\") ~ ",
                paste(xnames, collapse= "+")))
    mdl <- crq(fmla,data=df_train,taus = alpha,method = "Pow")
    res <- predict(mdl,
            newdata = df_test,
            type = "quantile")
    output$pow.bnd <- as.vector(res)
    cat("Done.\n")

    ## Portnoy 
    df_p = model.matrix(~ 0 + ., df[colnames(df) %in% xnames])
    df_train_p = df_p[train_idx, ]
    df_train_p = cbind(df_train_p, C = df_train$C, censored_T = df_train$censored_T, event = df_train$event)
    df_train_p = as.data.frame(df_train_p)
    df_test_p = df_p[-train_idx, ]
    xnames_p = colnames(df_p)

    cat("Compute the result of portnoy...")
    fmla <- as.formula(paste("Surv(censored_T, event) ~ ",
                paste(xnames_p, collapse= "+")))
    mdl <- crq(fmla, data = df_train_p, method = "Portnoy")
    mdl_coef <- coef(mdl,taus = alpha)
    res <- as.matrix(df_test_p) %*% 
    mdl_coef[-1] + mdl_coef[1] 
    output$por.bnd <- as.vector(res)
    cat("Done.\n")

    ## Peng and Huang
    cat("Compute the result of Peng and Huang...")
    fmla <- as.formula(paste("Surv(censored_T, event) ~ ",
                paste(xnames_p, collapse= "+")))
    mdl <- crq(fmla, data = df_train_p, method = "PengHuang")
    mdl_coef <- coef(mdl,taus = alpha)
    res <- as.matrix(df_test_p) %*% 
    mdl_coef[-1] + mdl_coef[1] 
    output$penghuang.bnd <- as.vector(res)
    cat("Done.\n")


    ## random forest
    cat("Compute the result of random forest...")
    df_rf = model.matrix(~ 0 + ., df[colnames(df) %in% xnames])
    df_train_rf = df_rf[train_idx, ]
    df_train_rf = cbind(df_train_rf, censored_T = df_train$censored_T, 
                        event = df_train$event)
    df_train_rf = as.data.frame(df_train_rf)
    df_test_rf = df_rf[-train_idx, ]
    df_test_rf = as.data.frame(df_test_rf)
    ntree <- 1000
    nodesize <- 80
    #colnames(df_train_rf)
    #colnames(df_test_rf)
    fmla <- as.formula(paste("censored_T ~ ",
                paste(xnames, collapse= "+")))
    mdl <- crf.km(fmla,
                ntree = ntree, 
                nodesize = nodesize,
                data_train = df_train_rf, 
                data_test = df_test_rf, 
                yname = 'censored_T', 
                iname = 'event',
                tau = alpha,
                method = "grf")
    output$rf.bnd  <- mdl$predicted
    #output$rf.bnd  <- rep(NA, n_test)
    cat("Done.\n")


    res <- cfsurv_train(p = p, c_list = c, Xtrain = df_train[, 1:p], C = df_train$C,
                        event = df_train$event, time = df_train$censored_T,
                        alpha = alpha, model="distBoost_train", seed = seed)

    #saveRDS(res, "trained_model2.rds")
    #res <- readRDS("trained_model1.rds")


    output_our <- cfsurv_predict(x = df_test[, 1:p], c = res$c,
                    alpha= alpha,
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
    
    output$our.bnd <- output_our
    #output$our.bnd <- rep(NA, n_test)
    
    #ret = cbind(output, output_cox)
    #write.csv(ret, paste0("output_HF/Mar5-seed-", seed, ".csv"), row.names=FALSE)

    #plot(output, output_cox)
    #hist((output - output_cox)/output_cox)
    #mean(output - output_cox/output_cox)
    #hist(output)

    #output_cox
    #output
    #predicted_lbd[-train_idx, i] = output
    #predicted_lbd_cox[-train_idx, i] = output_cox
    na_counts <- colSums(is.na(output))
    print(na_counts)
    output[is.na(output)] <- Inf
    #if (any(is.na(output_cox))) {
    #    print(paste0("Number of NA in output_cox:", sum(is.na(output_cox))))
    #    output_cox[is.na(output_cox)] <- Inf
    #}

    for (j in 1:7) {
        cov = sum(T[-train_idx] > output[, j])/n_test
        coverage[i, j] = cov
        print(paste0("Coverage for ", colnames(output)[j], ":", cov))
    }
    #cox_cov[i] = sum(T[-train_idx] > output_cox)/n_test
    #our_cov[i] = sum(T[-train_idx] > output)/n_test
    #print(paste0("Cox coverage:", cox_cov[i]))
    #print(paste0("Our coverage:", our_cov[i]))
}
print(colMeans(coverage))

output_table = as.matrix(output)
output_table = cbind(df_test, T[-train_idx], output_table)
colnames(output_table) = c(colnames(df_test), "T", colnames(output))

#write.csv(output_table, paste0("HF_Project/output_June4/output_June4_alpha_005_seed_", seed, ".csv"), row.names=FALSE)
#write.csv(predicted_lbd, "June4_predicted_lbd.csv", row.names=FALSE)
#write.csv(predicted_lbd_cox, "June4_predicted_lbd_cox.csv", row.names=FALSE)
#write.csv(cbind(cox_cov,our_cov), "June4_coverage.csv", row.names=FALSE)

#print("Cox coverage:")
#print(mean(cox_cov))
#print("Our coverage:")
#print(mean(our_cov))

write.csv(output_table, paste0("HF_Project/output_July18/July18_output_alpha_01_seed_", seed, ".csv"), row.names=FALSE)


#res_cox <- cfsurv_train(p = p, c_list = c, Xtrain = X[1:10000, ], C = C[1:10000],
#                    event = event[1:10000], time = censored_T[1:10000],
#                    alpha = 0.1, model= "cox_train")

#print(names(res_cox))

#output_cox <- cfsurv_predict(x = X[10001:n, ], c = res_cox$c,
#                   alpha=0.1,
#                   type="quantile",
#                   seed = 24601,
#                   model = "cox_predict",
#                   dist= "weibull",
#                   I_fit = NULL,
#                   ftol=.1,tol=.1,
#                   n.tree=100, 
#                   weight_calib = res_cox$weight_calib, 
#                   gbm_mdl = res_cox$gbm_mdl, 
#                   mdlrb = res_cox$mdlrb, 
#                   res_fit = res_cox$res_fit, 
#                   score = res_cox$score, 
#                   mdl = res_cox$mdl,
#                   bw = res_cox$bw
#                   )


# Fit Cox model
train_idx = sample(1:n, n_train)
df_train = df1[train_idx, ]
df_test = df1[-train_idx, ]
fit <- coxph(Surv(censored_T, event) ~ ., data = df1)
fit$
summary(fit)

head(df)

mean_event_hypertension_true <- mean(df$event[df$X10 == "True"])
mean_event_hypertension_false <- mean(df$event[df$X10 == "False"])



