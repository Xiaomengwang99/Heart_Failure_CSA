suppressPackageStartupMessages(library(conTree))
suppressPackageStartupMessages(library(grf))
suppressPackageStartupMessages(library(survival))
source("HF_Project/cfsurv_paper-main/utils/source.R")

#seed = Sys.getenv('SGE_TASK_ID')
#seed = as.integer(seed)
seed = 123
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
p = dim(X)[2]
xnames = paste0("X", 1:p)
colnames(df) <- c(xnames, "C", "censored_T", "event")


# Check for missing data in df
#missing_data <- colSums(is.na(df))
#missing_data
#head(df)
#mean(event)


num_trials = 20

predicted_lbd = matrix(0, nrow = n, ncol = num_trials)
predicted_lbd_cox = matrix(0, nrow = n, ncol = num_trials)

cox_cov = rep(0, num_trials)
our_cov = rep(0, num_trials)

n_train = 8000
n_test = n - n_train

alpha = 0.1

for (i in 1:num_trials){
    print(paste0("Iteration:", i))
    set.seed(i)

    train_idx = sample(1:n, n_train)
    df_train = df[train_idx, ]
    #saveRDS(df_train[, 1:p], "df_train.rds")
    df_test = df[-train_idx, ]
    
    column_types <- sapply(df_test, class)

    c = median(C)

    #hist(C)
    #abline(v = 365, col = "red")

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
    
    #output_cox2 = cfsurv(x = df_test[, 1:p], c_list = c, X = df_train[, 1:p], 
    #                    C = df_train$C, event = df_train$event, 
    #                    time = df_train$censored_T,
    #                    alpha= alpha,
    #                    model="cox", seed = seed)

    res <- cfsurv_train(p = p, c_list = c, Xtrain = df_train[, 1:p], C = df_train$C,
                        event = df_train$event, time = df_train$censored_T,
                        alpha = alpha, model="distBoost_train", seed = seed)

    #saveRDS(res, "trained_model2.rds")
    #res <- readRDS("trained_model1.rds")


    output <- cfsurv_predict(x = df_test[, 1:p], c = res$c,
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
    if (any(is.na(output_cox))) {
        print(paste0("Number of NA in output_cox:", sum(is.na(output_cox))))
        output_cox[is.na(output_cox)] <- Inf
    }
    cox_cov[i] = sum(T[-train_idx] > output_cox)/n_test
    our_cov[i] = sum(T[-train_idx] > output)/n_test
    print(paste0("Cox coverage:", cox_cov[i]))
    print(paste0("Our coverage:", our_cov[i]))
}

#write.csv(predicted_lbd, "June4_predicted_lbd.csv", row.names=FALSE)
#write.csv(predicted_lbd_cox, "June4_predicted_lbd_cox.csv", row.names=FALSE)
#write.csv(cbind(cox_cov,our_cov), "June4_coverage.csv", row.names=FALSE)

print("Cox coverage:")
print(mean(cox_cov))
print("Our coverage:")
print(mean(our_cov))

write.csv(cbind(cox_cov,our_cov), "July5_coverage2.csv", row.names=FALSE)


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

#sim_data <- data.frame(insurance = sample(levels(X[, "insurance"]), n, replace = TRUE),
#                       language = sample(levels(X[, "language"]), n, replace = TRUE),
#                       marital_status = sample(levels(X[, "marital_status"]), n, replace = TRUE),
#                       los = rnorm(n, mean = mean(X[, "los"]), sd = 2*sd(X[, "los"])),
#                       gender = sample(levels(X[, "gender"]), n, replace = TRUE),
#                       anchor_age = rnorm(n, mean = mean(X[, "anchor_age"]), sd = 2*sd(X[, "anchor_age"])))

#sim_data[, "insurance"] = as.factor(sim_data[, "insurance"])
#sim_data[, "language"] = as.factor(sim_data[, "language"])
#sim_data[, "marital_status"] = as.factor(sim_data[, "marital_status"])
#sim_data[, "gender"] = as.factor(sim_data[, "gender"])

#sapply(sim_data, class)

#ret <- cfsurv(x = X[10001:n, ],c_list = c, X = X[1:10000,], C = C[1:10000],event = event[1:10000], time = censored_T[1:10000] ,alpha= 0.1,model="distBoost_train")


#sim_data$predicted = res
#sim_data1 <- sim_data[sim_data$insurance == "Medicare", ]
#sim_data1 <- sim_data1[sim_data1$language == "ENGLISH", ]
#sim_data1 <- sim_data1[sim_data1$marital_status == "MARRIED", ]
#sim_data1 <- sim_data1[sim_data1$gender == "M", ]
#sim_data1 <- sim_data1[sim_data1$los < 10, ]
#sim_data1


#plot(sim_data1[, "anchor_age"], sim_data1$predicted)

