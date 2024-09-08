ret_mat_01 = matrix(NA, 250, 7)
ret_mat_005 = matrix(NA, 250, 7)
ret_mat_001 = matrix(NA, 250, 7)



for (i in 1:250){
  #if(i ==15 | i == 18|i ==24|i ==25|i ==27|i ==28|i ==32|i ==34|i ==38|i ==44|i ==45|i ==48|i ==49|i ==164|i ==176|i ==187){
  #  ret_mat[i,] = c(NA, NA, NA, NA, NA, NA, NA)
  #  next
  #}
    ret = read.csv(paste("HF_Project/output_July8/July8_coverage_seed_", i, ".csv", sep=""))
    ret = as.matrix(ret)
    ret_mat_01[i,] = ret
    
    ret = read.csv(paste("HF_Project/output_July18/July18_coverage_alpha_005_seed_", i, ".csv", sep=""))
    ret = as.matrix(ret)
    ret_mat_005[i,] = ret
    
    ret = read.csv(paste("HF_Project/output_July18/July18_coverage_alpha_001_seed_", i, ".csv", sep=""))
    ret = as.matrix(ret)
    ret_mat_001[i,] = ret
}
colnames(ret_mat_01) = c("Cox", "AFT", "Powell", 
            "Portnoy", "Peng+Huang", "RF", "CSA")
colnames(ret_mat_005) = c("Cox", "AFT", "Powell", 
                      "Portnoy", "Peng+Huang", "RF", "CSA")
colnames(ret_mat_001) = c("Cox", "AFT", "Powell", 
                      "Portnoy", "Peng+Huang", "RF", "CSA")
ret_mat_01 = ret_mat_01[, -5]
ret_mat_005 = ret_mat_005[, -5]
ret_mat_001 = ret_mat_001[, -5]
#if there is any NA, remove the row
#ret_mat = ret_mat[complete.cases(ret_mat),]
#ret_mat = ret_mat[, -c(1, 4)]

par(mfrow = c(1, 3))
boxplot(ret_mat_01, col = "lightblue", 
  xlab = "Method", ylab = "Empirical coverage", ylim = c(0.82, 0.92), main = "(a)")
abline(h = 0.9, col = "red")

boxplot(ret_mat_005, col = "lightblue", 
  xlab = "Method", ylim = c(0.93, 0.97), main = "(b)")
abline(h = 0.95, col = "red")

boxplot(ret_mat_001, col = "lightblue", 
  xlab = "Method", ylim = c(0.98, 1), main = "(c)")
abline(h = 0.99, col = "red")


##############################################################################
ret_mat_01 = matrix(NA, 250, 7)
ret_mat_005 = matrix(NA, 250, 7)
ret_mat_001 = matrix(NA, 250, 7)



for (i in 1:250){
  if(i ==15 | i == 18|i ==24|i ==25|i ==27|i ==28|i ==32|i ==34|i ==38|i ==44|i ==45|i ==48|i ==49|i ==164|i ==176|i ==187){
    ret_mat_01[i,] = c(NA, NA, NA, NA, NA, NA, NA)
  }else{
    ret = read.csv(paste("HF_Project/output_July18/July18_coverage_large_alpha01_seed_", i, ".csv", sep=""))
    ret = as.matrix(ret)
    ret_mat_01[i,] = ret
  }
  
  ret = read.csv(paste("HF_Project/output_July18/July18_coverage_large_alpha005_seed_", i, ".csv", sep=""))
  ret = as.matrix(ret)
  ret_mat_005[i,] = ret
  
  ret = read.csv(paste("HF_Project/output_July18/July18_coverage_large_alpha001_seed_", i, ".csv", sep=""))
  ret = as.matrix(ret)
  ret_mat_001[i,] = ret
}
colnames(ret_mat_01) = c("Cox", "AFT", "Powell", 
                         "Portnoy", "Peng+Huang", "RF", "CSA")
colnames(ret_mat_005) = c("Cox", "AFT", "Powell", 
                          "Portnoy", "Peng+Huang", "RF", "CSA")
colnames(ret_mat_001) = c("Cox", "AFT", "Powell", 
                          "Portnoy", "Peng+Huang", "RF", "CSA")
ret_mat_01 = ret_mat_01[, -5]
ret_mat_005 = ret_mat_005[, -5]
ret_mat_001 = ret_mat_001[, -5]
#if there is any NA, remove the row
ret_mat_01 = ret_mat_01[complete.cases(ret_mat_01),]


par(mfrow = c(1, 3))
boxplot(ret_mat_01, col = "lightblue", 
        xlab = "Method", ylab = "Empirical coverage", ylim = c(0, 0.92), main = "(a)")
abline(h = 0.9, col = "red")

boxplot(ret_mat_005, col = "lightblue", 
        xlab = "Method", ylim = c(0.6, 0.97), main = "(b)")
abline(h = 0.95, col = "red")

boxplot(ret_mat_001, col = "lightblue", 
        xlab = "Method", ylim = c(0.97, 1), main = "(c)")
abline(h = 0.99, col = "red")


ret_mat_01 = ret_mat_01[, -c(1, 4)]
ret_mat_005 = ret_mat_005[, -c(1, 4)]
ret_mat_001 = ret_mat_001[, -c(1, 4)]

par(mfrow = c(1, 3))
boxplot(ret_mat_01, col = "lightblue", main = "(a)",
        xlab = "Method", ylab = "Empirical coverage", ylim = c(0.875, 0.935))
abline(h = 0.9, col = "red")

boxplot(ret_mat_005, col = "lightblue", main = "(b)",
        xlab = "Method", ylim = c(0.93, 0.97))
abline(h = 0.95, col = "red")

boxplot(ret_mat_001, col = "lightblue", main = "(c)",
        xlab = "Method", ylim = c(0.97, 1))
abline(h = 0.99, col = "red")

################################################################################
temp = read.csv(paste("HF_Project/output_July18/July18_output_alpha_005_seed_", 1, ".csv", sep=""))
output = as.matrix(temp)

for (i in 2:250){
  #if(i ==15 | i == 18|i ==24|i ==25|i ==27|i ==28|i ==32|i ==34|i ==38|i ==44|i ==45|i ==48|i ==49|i ==164|i ==176|i ==187){
  #  ret_mat[i,] = c(NA, NA, NA, NA, NA, NA, NA)
  #  next
  #}
  ret = read.csv(paste("HF_Project/output_July18/July18_output_alpha_005_seed_", i, ".csv", sep=""))
  ret = as.matrix(ret)
  output = rbind(output, ret)
}

dim(output)

cox = as.numeric(output[, 28])
aft = as.numeric(output[, 29])
pow = as.numeric(output[, 30])
port = as.numeric(output[, 31])
#ph = as.numeric(output[, 32])
rf = as.numeric(output[, 33])
csa = as.numeric(output[, 34])

df = data.frame(cox, aft, pow, port, rf, csa)
colnames(df) = c("Cox", "AFT", "powell", 
                 "Portnoy", "RF", "CSA")
head(df)


num_inf = apply(df, 2, function(x) sum(is.infinite(x)))

replace(df, is.infinite(as.matrix(df)), NA)


quantiles <- apply(df, 2, quantile, probs = c(0.25, 0.5, 0.75))
quantiles_table <- t(quantiles)
colnames(quantiles_table) <- c("Q1", "Median", "Q3")
length_table = t(quantiles_table)
write.table(length_table, file = "/Users/xiaomengwang/Desktop/MT Survival/HF_Project/length_table.tex", sep = "&", col.names = TRUE, row.names = TRUE, quote = FALSE)

############################################################################

head(output)
True_T = as.numeric(output[, 27])
cox_cov = cox < True_T
aft_cov = aft < True_T
pow_cov = pow < True_T
port_cov = port < True_T
#ph_cov = ph < True_T
rf_cov = rf < True_T
csa_cov = csa < True_T
insurance = as.factor(output[, 1])
language = as.factor(output[, 2])
marital = as.factor(output[, 3])
los = as.numeric(output[, 4])
gender = as.factor(output[, 5])
age = as.integer(output[, 6])
age_group <- cut(age, breaks = c(0, 18, 25, 35, 45, 55,65,Inf), 
                    labels = c("0-18", "18-24", "25-34", "35-44", "45-54", "55-64", "65+"))
race = as.factor(output[, 7])

df_cov = data.frame(insurance, language, marital, los, gender, age,  age_group, race, cox_cov, aft_cov, pow_cov, port_cov, rf_cov, csa_cov)
head(df_cov)

race_mean <- matrix(NA, nrow = length(levels(df_cov$race)), ncol = 6)
for (i in 1:length(levels(df_cov$race))) {
  race <- levels(df_cov$race)[i]
  race_data <- df_cov[df_cov$race == race, ]
  race_mean[i, ] <- colMeans(race_data[, 9:14], na.rm = TRUE)
}
rownames(race_mean) <- levels(df_cov$race)
colnames(race_mean) <- c("Cox", "AFT", "Powell", "Portnoy", "RF", "CSA")
race_mean
library(ggplot2)
library(reshape2)
melted_data <- melt(race_mean)

p1 <- ggplot(data = melted_data, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  labs(x = "Race", y = "Method", title = "(a)")+
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.key.size = unit(0.5, "cm"),        # Adjusts the size of the legend keys
    legend.text = element_text(size = 8),     # Adjusts the size of the legend text
    legend.title = element_text(size = 10)    # Adjusts the size of the legend title
  )


gender_mean <- matrix(NA, nrow = length(levels(df_cov$gender)), ncol = 6)
for (i in 1:length(levels(df_cov$gender))) {
  gender <- levels(df_cov$gender)[i]
  gender_data <- df_cov[df_cov$gender == gender, ]
  gender_mean[i, ] <- colMeans(gender_data[, 9:14], na.rm = TRUE)
}
rownames(gender_mean) <- levels(df_cov$gender)
colnames(gender_mean) <- c("Cox", "AFT", "Powell", "Portnoy", "RF", "CSA")
gender_mean
melted_data_gender <- melt(gender_mean)
# Draw the heatmap
p2 <- ggplot(data = melted_data_gender, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  labs( x = "Sex", y = "", title = "(b)")+
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.key.size = unit(0.5, "cm"),        # Adjusts the size of the legend keys
    legend.text = element_text(size = 8),     # Adjusts the size of the legend text
    legend.title = element_text(size = 10)    # Adjusts the size of the legend title
  )


age_group_mean <- matrix(NA, nrow = length(levels(df_cov$age_group)), ncol = 6)
for (i in 1:length(levels(df_cov$age_group))) {
  age_group <- levels(df_cov$age_group)[i]
  age_group_data <- df_cov[df_cov$age_group == age_group, ]
  age_group_mean[i, ] <- colMeans(age_group_data[, 9:14], na.rm = TRUE)
}
rownames(age_group_mean) <- levels(df_cov$age_group)
colnames(age_group_mean) <- c("Cox", "AFT", "Powell", "Portnoy", "RF", "CSA")
age_group_mean
melted_data_age_group <- melt(age_group_mean)
# Draw the heatmap
p3 <- ggplot(data = melted_data_age_group, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  theme_minimal() +
  labs(x = "Age Group", y = "", title = "(c)")+
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.key.size = unit(0.5, "cm"),        # Adjusts the size of the legend keys
    legend.text = element_text(size = 8),     # Adjusts the size of the legend text
    legend.title = element_text(size = 10)    # Adjusts the size of the legend title
  )

library(gridExtra)
grid.arrange(p2, p1, p3, ncol = 3)


