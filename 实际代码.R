# Load necessary libraries
library(sp)
library(spgwr)
library(spdep)
library(readxl)
library(spatialreg)
library(gstat)
library(stats)
library(dbscan)
library(cluster)
library(ggplot2)
library(scatterplot3d)
library(rgl)
library(mclust)
library(plotly)
library(dplyr)
library(FNN)
library(matrixStats)
library(tidyr)
library(scales)
library(factoextra)
library(spatstat)
library(clValid)

# Load data with spatial coordinates and mutation rates
data<-as.data.frame(read.csv("data1.csv"))
data<-na.omit(data)
data_matrix<-as.matrix(data)
data2<-data[1:972,c("X","Y","Z","virusPercent")]
log_virus_percent <- log(data2$virusPercent)

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
#这两个变量之间存在显著的相关性。

#特征向量+Kmean聚类法模型

#抽提坐标数据
xyz<-data2[,c("X","Y","Z")]
#标准化xyz数据
xyz_norm<-scale(xyz)
# 将标准化处理过的坐标和突变频率组合为特征向量
feature_vector <- cbind(xyz_norm, data2$virusPercent)
# 进行 KMeans 聚类
set.seed(123)
kmeans_result <- kmeans(feature_vector, centers = 3, nstart = 20)
# 将聚类结果添加到数据框中
data2$cluster <- as.factor(kmeans_result$cluster)
# 可视化聚类结果
fviz_cluster(list(data = feature_vector, cluster = kmeans_result$cluster))
# 比较聚类结果和突变频率
cluster_freq <- data2 %>%
  group_by(cluster) %>%
  summarise(mean_freq = mean(virusPercent), count = n())

total_freq <- mean(data2$virusPercent)

print(cluster_freq)
print(total_freq)

#Permutation法(排名百分数假设检验)
# 计算原始数据的聚类结果
orig_clusters <- kmeans(cbind(xyz_norm,data2$virusPercent), 3)
orig_means <- tapply(data2[,ncol(data2)], orig_clusters$cluster, mean)
#这一步是为了计算原始数据中每个类的平均突变频率所处的百分位数。
#具体地，中包含原始数据中每个类突变频率的均值，那么 rank(orig_means) 返回原始数据中每个类的突变频率均值的排名。
#然后，通过除以 length(orig_means)，即原始数据被分成的类数，可以计算出每个类的均值的排名的比例。
#最后，乘以 100 并四舍五入到两位小数，即可得到排名。
orig_pcts <- round(rank(orig_means)/length(orig_means)*100, 2)

# 进行1000次突变概率重排并计算每次聚类结果
n_permutations <- 1000
perm_results <- matrix(NA, n_permutations, length(orig_means))
for (i in 1:n_permutations) {
  perm_data <- data2
  perm_data[,ncol(data2)] <- sample(data2[,ncol(data2)], replace=FALSE)
  perm_clusters <- kmeans(cbind(xyz_norm,perm_data[,ncol(data2)]), 3)
  perm_means <- tapply(perm_data[,ncol(data2)], perm_clusters$cluster, mean)
  perm_pcts <- round(rank(perm_means)/length(perm_means)*100, 2)
  #这行代码是用来计算每个聚类在重排结果中的百分位数。
  #具体来说，rank()函数计算了每个聚类均值在重排结果中的排名，
  #然后除以重排结果的总数，并乘以100，
  #得到了每个聚类在重排结果中的百分位数。
  #最后，round()函数将结果保留两位小数并四舍五入。
  #这个百分位数可以用来评估原始聚类结果中每个聚类的可能性是否显著高于随机聚类结果，
  #在第六步计算p值时会用到。
  perm_results[i,] <- perm_pcts
  #将每次重排得到的聚类结果的百分位数保存在perm_results矩阵中。
  #具体来说，perm_results矩阵的每一行都对应一次重排结果，
  #每一列对应一个聚类的百分位数。
  #在for循环中，第i次重排得到的三个聚类的百分位数保存在perm_pcts向量中，
  #然后将perm_pcts向量赋值给perm_results矩阵的第i行，
  #从而将这次重排的聚类结果存储在perm_results矩阵中。
}
# 计算每个聚类的p值
#检验原始数据中每个类的均值是否显著高原始数据的均值。
#具体来说，这个假设检验的零假设是：原始数据中每个类的均值与原始数据的均值相等，即不存在显著差异。
#备择假设是：原始数据中每个类的均值显著高于原始数据中的均值
p_values <- numeric(length(orig_means))
for (i in 1:length(orig_means)) {
  p_values[i] <- mean(perm_results[,i] >= orig_pcts[i])
}
# 输出结果
cat("原始聚类结果的均值：", orig_means, "\n")
cat("原始聚类结果的百分位数：", orig_pcts, "\n")
cat("p值：", p_values, "\n")










#Permutation法(均值假设检验)
# 计算原始聚类结果中每个聚类的均值
orig_means <- tapply(data2[, ncol(data2)], orig_clusters$cluster, mean)

# 确定突变高发区
high_freq_cluster <- which.max(orig_means)

# 初始化向量来保存每次重排后的最高平均值
max_means <- rep(NA, n_permutations)
p_v<-numeric(n_permutations)
p_v_1<-numeric(n_permutations)

# 进行一次1000次重排
for (i in 1:n_permutations) {
  # 对突变概率进行重排
  perm_data[,ncol(data2)] <- sample(data2[,ncol(data2)], replace=FALSE)
  
  # 对重排后的数据进行聚类
  perm_clusters <- kmeans(cbind(xyz_norm,perm_data[,ncol(data2)]), 3)
  
  # 计算每个聚类的平均突变频率
  perm_means <- tapply(perm_data[, ncol(perm_data)], perm_clusters$cluster, mean)
  
  # 保存最高平均值
  max_means[i] <- max(perm_means)
}

# 计算原始聚类结果中突变高发区的平均值在排序后的位置
#我们使用 rank 变量计算原始聚类结果中突变高发区的平均值在重排结果中的排名。
#具体来说，我们将所有重排结果中的最高平均值进行排序，
#并计算突变高发区的原始平均值在排序后的位置。
#如果突变高发区的原始平均值在排序后的位置较靠前，
#那么它在重排结果中的表现就比较好，
#表明突变高发区是有统计学意义的；
#那么我们就可以拒绝原始聚类结果中突变高发区的平均值在随机重排中产生的最大值的排序位置上的假设，
#即突变高发区的平均值是随机出现的。
#反之，如果它在排序后的位置较靠后，
#那么它在重排结果中的表现就比较差，
#表明突变高发区的结果可能是随机出现的。
rank <- sum(max_means >= orig_means[high_freq_cluster]) / n_permutations
#计算出超过原始数据中聚类结果中平均突变频率最高的聚类的平均值的随机重排数量
#将这个数量除以随机重排的次数num_permutations，
#得到原始聚类结果中平均突变频率最高的聚类的平均值在所有随机重排中出现的排序百分比rank。
print(rank)
#rank 的值为0.206，这意味着在 1000 次重排的结果中，
#有大约 20.6% 的结果中的突变高发区的平均值比原始聚类结果中的突变高发区的平均值更高，
#而有 79.4% 的结果中的突变高发区的平均值比原始聚类结果中的突变高发区的平均值更低。

# 计算p值
#在零假设下，我们假设原始突变高发区与其他区域之间没有显著的差异，
#也就是说原始突变高发区的平均值与其他区域的平均值相同。
#我们随机重排了数据，重新计算了每个区域的平均值，并记录每个区域在每个重排结果中的排名。
#然后，我们计算了原始突变高发区在所有重排结果中的排名比例，并将其与 1 减去排名比例中较小的那个值进行比较。
#如果这个比例小于显著性水平（例如0.05），那么我们可以拒绝零假设，认为突变高发区是有统计学意义的。
#在计算 rank 时，我们对原始聚类结果中的突变高发区进行了随机重排，
#重新计算了每个聚类的平均值，并记录每个聚类在每个重排结果中的平均值。
#然后，我们计算了原始聚类结果中突变高发区的平均值在所有重排结果中的排名比例，作为 rank 值。
#这种方法的科学性依据是基于随机化原理，即我们通过随机化来模拟零假设，
#即突变高发区的平均值与其他区域的平均值相同，并比较观察到的结果与随机化结果的分布，
#以评估观察结果的显著性。
p_value <- min(rank, 1-rank)
print(p_value)

#为了进一步验证原始聚类的非随机性，进行第二种p值计算
#我们使用 rank 值来计算双尾 p 值。
#具体而言，我们将 rank 值标准化为标准正态分布的分位数，并计算两个尾部的概率。
#这个 p 值告诉我们，如果原始数据中聚类结果中平均突变频率最高的聚类的平均值是完全随机的，那么它会在所有随机重排中出现的位置的概率是多少。
#这个 p 值可以用来评估原始聚类结果中的最大值是否显著高于随机重排的最大值，
#从而判断是否可以拒绝原始聚类结果中平均突变频率最高的聚类的平均值是随机出现的假设。
p_value_1 <- 2 * (1 - pnorm(abs(rank - 0.5) * sqrt(n_permutations / 12)))
print(p_value_1)

#上面只展示了进行一次1000次重排的代码，但由于Kmean聚类的随机性，
#我们每次会得到一个不一样或者一样的p值，所以单独一次p值大于或小于0.05不能绝对说明原始聚类的突变高发区是否和其他区域的突变频率直接存在显著差异。
#因此我打算进行100次每次1000次重排的实验，然后取所有得到的p值的平均值。
m <- 100
for (i in 1:m) {
  for (j in 1:n_permutations){
    # 对突变概率进行重排
    perm_data[,ncol(data2)] <- sample(data2[,ncol(data2)], replace=FALSE)
    
    # 对重排后的数据进行聚类
    perm_clusters <- kmeans(cbind(xyz_norm,perm_data[,ncol(data2)]), 3)
    
    # 计算每个聚类的平均突变频率
    perm_means <- tapply(perm_data[, ncol(perm_data)], perm_clusters$cluster, mean)
    
    # 保存最高平均值
    max_means[j] <- max(perm_means)
    
  }
  rank <- sum(max_means >= orig_means[high_freq_cluster]) / n_permutations
  p_v[i] <- min(rank, 1-rank)
  p_v_1[i]<-2 * (1 - pnorm(abs(rank - 0.5) * sqrt(n_permutations / 12)))
  
}

print(mean(p_v))
print(mean(p_v_1))

#设置显著性水平
alpha<-0.05

# 判断是否拒绝原假设
if (mean(p_v) < alpha) {
  cat("Reject null hypothesis with average p-value =", mean(p_v))
} else {
  cat("Fail to reject null hypothesis with average p-value =", mean(p_v))
}
if (mean(p_v_1) < alpha) {
  cat("Reject null hypothesis with average p-value =", mean(p_v_1))
} else {
  cat("Fail to reject null hypothesis with average p-value =", mean(p_v_1))
}