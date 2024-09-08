library(ggplot2)
cov_lst = read.csv("July5_coverage.csv", header = TRUE)
cov_lst2 = read.csv("July5_coverage2.csv", header = TRUE)


cov_lst = rbind(cov_lst, cov_lst2)
cov_lst
colnames(cov_lst) = c("Cox", "Conformal")

df = read.csv("df_predicted_lbd.csv")
head(df)

p1 = boxplot(cov_lst, beside = TRUE, names.arg = c("Cox", "Conformal"), 
    xlab = "Method", ylab = "Coverage", main = "Coverage Comparison")
abline(h = 0.90, col = "red")

df$X.gender <- ifelse(df$X.gender == "F", "Female", "Male")

par(mfrow = c(1, 3), mar = c(4, 4, 2, 1) + 0.1, oma = c(0, 0, 2, 0))
boxplot(df$predicted_lbd ~ df$X.gender, xlab = "Sex", ylab = "Predicted Lower Bound (days)", col = "lightblue", main = "(a)")


boxplot(df$predicted_lbd ~ df$X.race, col = "lightblue", xlab = "Race", ylab = "", main = "(b)")

df$age_group <- cut(df$X.anchor_age, breaks = c(0, 18, 25, 35, 45, 55,65,Inf), 
                    labels = c("0-18", "18-24", "25-34", "35-44", "45-54", "55-64", "65+"))

boxplot(df$predicted_lbd ~ df$age_group, col = "lightblue", xlab = "Age group", ylab = "", main = "(c)")

summary(as.factor(df$X.race))




