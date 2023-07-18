# 生成随机的三维坐标数据
data <- read.csv("data1.csv")
sub_data<-data[1:972,c("X","Y","Z","virusPercent")]
sub_data_1<-data[1:972,c("X","Y","Z")]
# 计算每个维度的基本统计信息
summary(sub_data)

# 绘制散点图
library(plot3D)
scatter3D(sub_data[, 1], sub_data[, 2], sub_data[, 3], col = "blue")

# 绘制直方图
hist(sub_data[, 1], breaks = 20)
hist(sub_data[, 2], breaks = 20)
hist(sub_data[, 3], breaks = 20)

# 绘制核密度估计图
plot(density(sub_data[, 1]))
plot(density(sub_data[, 2]))
plot(density(sub_data[, 3]))
# 绘制QQ图
qqnorm(sub_data[, 1])
qqline(sub_data[, 1])
qqnorm(sub_data[, 2])
qqline(sub_data[, 2])
qqnorm(sub_data[, 3])
qqline(sub_data[, 3])

#计算每个维度的均值和标准差
apply(sub_data,2,mean)
apply(sub_data,2,sd)

# 进行双样本t检验
t.test(sub_data[, 1], sub_data[, 2])
t.test(sub_data[, 1], sub_data[, 3])
t.test(sub_data[, 2], sub_data[, 3])
#如果p值很小，说明两个维度之间的均值和标准差之间存在显著差异，
#需要进行标准化处理等进一步操作

#方差分析空间差异性
data_df<-as.data.frame(sub_data_1)
#进行方差分析
#相关模型
model_1 <- aov(X ~Y + Z,data = data_df)
summary(model_1)

# 计算皮尔逊相关系数
cor_mat <- cor(data_df, method = "pearson")
# 查看相关系数矩阵
print(cor_mat)
#相关系数的范围在 -1 到 1 之间,
#绝对值越大则表示相关性越强，而正负号则表示相关性的方向
for (i in 1:(ncol(data_df) - 1)) {
  for (j in (i + 1):ncol(data_df)) {
    cor_test <- cor.test(data_df[, i], data_df[, j], method = "pearson")
    cat("correlation between", names(data_df)[i], "and", names(data_df)[j], "is", cor_test$estimate, "with p-value", cor_test$p.value, "\n")
  }
}

# 计算斯皮尔曼等级相关系数
cor_mat_1 <- cor(data_df, method = "spearman")
# 查看相关系数矩阵
print(cor_mat_1)
# 计算多个变量之间的斯皮尔曼等级相关系数及其显著性水平
for (i in 1:(ncol(data_df) - 1)) {
  for (j in (i + 1):ncol(data_df)) {
    cor_test <- cor.test(data_df[, i], data_df[, j], method = "spearman")
    cat("correlation between", names(data_df)[i], "and", names(data_df)[j], "is", cor_test$estimate, "with p-value", cor_test$p.value, "\n")
  }
}

# 计算肯德尔等级相关系数
cor_mat_2 <- cor(data_df, method = "kendall")
# 查看相关系数矩阵
print(cor_mat_2)
# 计算多个变量之间的肯德尔等级相关系数及其显著性水平
for (i in 1:(ncol(data_df) - 1)) {
  for (j in (i + 1):ncol(data_df)) {
    cor_test <- cor.test(data_df[, i], data_df[, j], method = "kendall")
    cat("correlation between", names(data_df)[i], "and", names(data_df)[j], "is", cor_test$estimate, "with p-value", cor_test$p.value, "\n")
  }
}

#如果显著性水平非常小（通常小于 0.05），
#则我们可以认为相关性是显著的，也就是说，在总体水平下，
#这两个变量之间存在显著的相关性
#根据皮尔逊和肯德尔相关性系数及显著性水平可看出，
#X和Y之间并未存在较大的相关性，需进一步探究

