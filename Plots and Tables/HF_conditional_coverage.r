
df = read.csv("June4_predicted_lbd.csv")
df_cox = read.csv("June4_predicted_lbd_cox.csv")

data = read.csv("final_df_May10.csv")
T = data$dod_since_disch
T = ifelse(is.na(data$dod_since_disch), Inf, data$dod_since_disch + data$dischtime_since_admission)

colnames(data)
n = nrow(df)
cov_lst = rep(0, n)
cov_lst_cox = rep(0, n)

lbd_lst = rep(0, n)
lbd_lst_cox = rep(0, n)

#mean(as.vector(df[1, ]-df_cox[1, ]))

for (i in 1:n){
    cov_lst[i] = mean(df[i, ][df[i, ] != 0] < T[i])
    cov_lst_cox[i] = mean(df_cox[i, ][df_cox[i, ] != 0] < T[i])
    lbd_lst[i] = mean(df[i, ][df[i, ] != 0])
    lbd_lst_cox[i] = mean(df_cox[i, ][df_cox[i, ] != 0])
}

mean(cov_lst)
mean(cov_lst_cox)

mean(lbd_lst - lbd_lst_cox)
hist(lbd_lst - lbd_lst_cox)

lbd_diff  = c()

#for (i in 1:n){
#  for (j in 1:100){
#    if (df[i, j] != 0){
#      diff = (df[i, j] - df_cox[i, j])
#      lbd_diff = c(lbd_diff, diff)
#    }
#  }
#}

#hist(lbd_diff)
#mean(lbd_diff)


data$cov = cov_lst
data$cox_cov = cov_lst_cox

boxplot(data$cov ~ data$race, main = "Predicted lower bound by sex", xlab = "Gender", ylab = "Conditional coverage")

boxplot(data$cox_cov ~ data$race, main = "Predicted lower bound by sex", xlab = "Gender", ylab = "Conditional coverage")

boxplot(T ~ data$race, main = "Predicted lower bound by sex", xlab = "Gender", ylab = "Survival time")


data$age_group <- cut(data$anchor_age, breaks = c(0, 18, 25, 35, 45, 55,65,Inf), 
                    labels = c("0-18", "18-24", "25-34", "35-44", "45-54", "55-64", "65+"))

boxplot(data$cov ~ data$age_group, main = "Predicted lower bound by age group", xlab = "Age group", ylab = "Predicted Lower Bound")

library(ggplot2)

heatmap_data <- aggregate(cov ~ age_group + race, data = data, FUN = mean)

ggplot(heatmap_data, aes(x = age_group, y = race, fill = cov)) +
  geom_tile() +
  labs(title = "Heatmap of Conditional Coverage by Age Group and Race",
       x = "Age Group",
       y = "Race",
       fill = "Conditional Coverage")

heatmap_data2 <- aggregate(cox_cov ~ age_group + race, data = data, FUN = mean)

ggplot(heatmap_data2, aes(x = age_group, y = race, fill = cox_cov)) +
  geom_tile() +
  labs(title = "Heatmap of Conditional Coverage (Cox) by Age Group and Race",
       x = "Age Group",
       y = "Race",
       fill = "Conditional Coverage")

heatmap_data3 <- aggregate(cov ~ gender + race, data = data, FUN = mean)

ggplot(heatmap_data3, aes(x = gender, y = race, fill = cov)) +
  geom_tile() +
  labs(title = "Heatmap of Conditional Coverage (Cox) by Age Group and Race",
       x = "Age Group",
       y = "Race",
       fill = "Conditional Coverage")

heatmap_data4 <- aggregate(cox_cov ~ gender + race, data = data, FUN = mean)

ggplot(heatmap_data4, aes(x = gender, y = race, fill = cox_cov)) +
  geom_tile() +
  labs(title = "Heatmap of Conditional Coverage (Cox) by Age Group and Race",
       x = "Age Group",
       y = "Race",
       fill = "Conditional Coverage")

