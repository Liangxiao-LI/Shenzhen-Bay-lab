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
t.test(sub_data[, 1], sub_data[, 4])
t.test(sub_data[, 2], sub_data[, 4])
t.test(sub_data[, 3], sub_data[, 4])
#如果p值很小，说明两个维度之间的均值和标准差之间存在显著差异，
#需要进行标准化处理等进一步操作

#方差分析空间差异性
data_df<-as.data.frame(sub_data_1)

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

#在进行标准化处理前，需要先观察数据在不同维度之间是否存在显著性差异，
#如果差异性较大，那么可以考虑进行标准化处理。
#但是，在进行标准化处理时，还需要观察标准化前后数据的分布形态和分布特征，
#以确定标准化处理是否适合该数据集.



#特征向量+Kmean聚类法模型

#抽提坐标数据
xyz<-data2[,c("X","Y","Z")]
#标准化xyz数据
xyz_norm<-scale(xyz)
# 将标准化处理过的坐标和突变频率组合为特征向量
feature_vector <- cbind(xyz_norm, data2$virusPercent)

#修改1 可视化标准化后的蛋白质结构
# 创建3D散点图
plot3d(feature_vector[,1], feature_vector[,2], feature_vector[,3], col = "blue", size = 2)
plot3d(data2[,1], data2[,2], data2[,3], col = "blue", size = 2)
#观察发现，原始数据是实际的972个点的三维坐标xyz和每个点的突变频率，可以画出来一个三维散点图，
#现在将xyz进行标准化处理后，得到了一个新的数据，又可以画出一个新的标准化后的三维散点图。
#前后两个散点图空间结构和形状上几乎相同。
#标准化处理会改变原始数据的数值范围和分布，但不会改变数据的结构和形状。
#因此，将原始数据进行标准化处理后，生成的新的标准化三维散点图与原始三维散点图相似，但数值范围和分布不同。

#原始的xyz和标准化后的xyz的数据范围和分布
# 查看原始数据的数据范围和分布情况
summary(data2[, 1:3])  # 显示前三列（x、y、z）的数据范围和分布情况
hist(data2[, 1])      # 绘制x列的直方图
hist(data2[, 2])      # 绘制y列的直方图
hist(data2[, 3])      # 绘制z列的直方图

# 查看标准化后的数据的数据范围和分布情况
summary(feature_vector[, 1:3])  # 显示前三列（x_norm、y_norm、z_norm）的数据范围和分布情况
hist(feature_vector[, 1])      # 绘制x_norm列的直方图
hist(feature_vector[, 2])      # 绘制y_norm列的直方图
hist(feature_vector[, 3])      # 绘制z_norm列的直方图

#验证标准化前后数据分布是否相似
#可以比较前后的直方图，观察它们的形状、中心位置、分散程度等特征是否相似
#如果两个直方图的形状和特征较为相似，那么它们的分布也可能相似
#另外，还可以使用一些统计方法来比较两个数据集的分布情况，例如Kolmogorov-Smirnov检验、Anderson-Darling检验，
#这些方法可以计算两个数据集之间的距离或差异程度，并给出相应的显著性水平，可以用来评估前后数据分布的相似性。

# 加载ks包，用于进行Kolmogorov-Smirnov检验
library(ks)

# 对比前后x列数据的分布
ks.test(data2[, "X"], feature_vector[, "X"])  # 进行Kolmogorov-Smirnov检验，比较x列的分布

# 对比前后y列数据的分布
ks.test(data2[, 2], feature_vector[, 2])  # 进行Kolmogorov-Smirnov检验，比较y列的分布

# 对比前后z列数据的分布
ks.test(data2[, 3], feature_vector[, 3])  # 进行Kolmogorov-Smirnov检验，比较z列的分布

#如果两个数据列的分布相似，则Kolmogorov-Smirnov检验的p值应该比较大（大于0.05），
#否则p值较小（小于0.05）。如果p值较小，则可以认为前后数据列的分布不相似

#如果标准化后的数据分布与原始数据分布不相似，
#那么将标准化后的数据进行聚类可能会导致结果不准确。
#因为聚类算法通常基于数据的距离或相似性进行聚类，如果数据分布不同，
#那么聚类结果可能会受到影响，导致聚类效果不佳。
#也就是说，如果标准化后数据的分布和标准化前的不相似的话，
#如果对标准化后的数据进行聚类，得到突变概率高发区所具有的点的话，
#不能直接说这些点在未标准化的数据中也对应突变高发区。

#另外，即使聚类结果在标准化后的数据上表现良好，也不一定能直接反映到原始数据上。
#因为标准化是一种线性变换，它可能会改变数据的分布、范围和形状等特征。
#如果将聚类结果直接应用到原始数据上，可能会导致聚类结果不准确。

#xyz三个值的数据类型是地理位置，而突变频率的数据类型是概率，
#将地理坐标和突变概率两种数据类型合在一起进行K均值聚类是不合适的，
#因为地理坐标和突变概率是不同类型的变量，它们具有不同的量纲和意义，
#合并在一起可能会导致聚类结果失真。
#为了解决这个问题，需要将地理坐标和突变概率转换为相同的尺度，
#例如将地理坐标转换为数值型变量，将突变概率转换为概率的对数或负数。

#如果数据已经是数值型的三维坐标，那么不需要再进行特殊的转换，可以直接使用这些数值型变量进行聚类分析，
#概率的对数或负数可以用来表示概率的大小，同时也可以用来比较不同概率之间的大小关系。
#例如，将突变概率p转换为概率的对数ln(p)或负数-np，可以将概率的范围从[0,1]转换为[-∞,0]或[0,∞]，
#使其与地理坐标在同一尺度下进行比较。
#这种转换方法可以将概率转换为数值型变量，使其与其他数值型变量一起参与聚类分析。

#坐标值类的数据可以被认为是数值型变量，是因为它们是由数值表达的。
#在空间坐标系中，每个点的位置都可以用一组数值来表示。
#经过标准化处理后的坐标值可以仍然被认为是数值型变量。
#标准化处理是一种线性变换，将原始数据映射到一个新的坐标系中，
#使得每个维度的数据具有相同的尺度和范围。
#在标准化处理后，坐标值的数值类型和数值特征都没有发生变化

#从密度图中观察密度曲线的形状和位置来判断标准化前后数据的分布差异性
library(ggplot2)
library(gridExtra)
# 原始数据的密度图
p1 <- ggplot(data2, aes(x = X)) +
  geom_density(fill = "steelblue", color = "white") +
  xlab("X") + ylab("Density") + ggtitle("Original X Density Distribution")

p2 <- ggplot(data2, aes(x = Y)) +
  geom_density(fill = "steelblue", color = "white") +
  xlab("Y") + ylab("Density") + ggtitle("Original Y Density Distribution")

p3 <- ggplot(data2, aes(x = Z)) +
  geom_density(fill = "steelblue", color = "white") +
  xlab("Z") + ylab("Density") + ggtitle("Original Z Density Distribution")
feature_vector<-data.frame(feature_vector)
# 标准化后的密度图
p4 <- ggplot(feature_vector, aes(x = X)) +
  geom_density(fill = "steelblue", color = "white") +
  xlab("X") + ylab("Density") + ggtitle("Standardized X Density Distribution")

p5 <- ggplot(feature_vector, aes(x = Y)) +
  geom_density(fill = "steelblue", color = "white") +
  xlab("Y") + ylab("Density") + ggtitle("Standardized Y Density Distribution")

p6 <- ggplot(feature_vector, aes(x = Z)) +
  geom_density(fill = "steelblue", color = "white") +
  xlab("Z") + ylab("Density") + ggtitle("Standardized Z Density Distribution")

# 将六张图合并为两行三列,形状，位置，峰度和偏度
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2, ncol = 3)

#如果在不同的变量上进行标准化，即使在每个变量上的密度图相似，
#也不能保证在所有变量上的分布相似。
#因此，如果使用标准化数据进行聚类分析，建议同时考虑多个变量的分布情况，而不仅仅是单个变量。
#如果在每个变量上的密度图都相似，但KS检验却拒绝原假设，可能是由于KS检验对样本量敏感，
#并且在样本量较大时可能会拒绝原假设。此外，KS检验可能会检测到一些细微的差异，这些差异可能在实际分析中并不重要。
#因此，还可以使用其他的检验方法，例如Anderson-Darling检验或Chi-squared检验等，以进一步探究数据的分布是否相似

# 循环绘制所有的二维多变量密度图
# 标准化前的二维密度图
library(ggplot2)
p1 <- ggplot(data2, aes(X, Y)) + 
  geom_density_2d(alpha = 0.6) +
  labs(x = "X", y = "Y") +
  theme_bw()

p2 <- ggplot(data2, aes(X, Z)) + 
  geom_density_2d(alpha = 0.6) +
  labs(x = "X", y = "Y") +
  theme_bw()

p3 <- ggplot(data2, aes(Y, Z)) + 
  geom_density_2d(alpha = 0.6) +
  labs(x = "Y", y = "Z") +
  theme_bw()
#标准化后的二维变量密度图
p4 <- ggplot(feature_vector, aes(X, Y)) + 
  geom_density_2d(alpha = 0.6) +
  labs(x = "X", y = "Y") +
  theme_bw()

p5 <- ggplot(feature_vector, aes(X, Z)) + 
  geom_density_2d(alpha = 0.6) +
  labs(x = "X", y = "Z") +
  theme_bw()

p6 <- ggplot(feature_vector, aes(Y, Z)) + 
  geom_density_2d(alpha = 0.6) +
  labs(x = "Y", y = "Z") +
  theme_bw()

# 将六张图合并为两行三列,形状，位置，峰度和偏度
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2, ncol = 3)

#将坐标和频率转化到同一尺度下
#将不同尺度的数据转化到同一尺度可以帮助我们消除不同类型的数据之间的差异，
#包括地理意义上的坐标和概率意义上的突变频率。
#这可以使得不同类型的数据在聚类分析中具有可比性和可解释性，从而提高聚类分析的准确性和有效性

#坐标值属于连续型数值尺度，通常表示物理空间中的位置或大小。
#在聚类分析中，坐标值通常被视为连续型数值特征，需要进行特征缩放或标准化处理，以便进行聚类分析。
#突变概率特征属于概率型数值尺度，通常表示某个事件发生的概率或可能性。
#在聚类分析中，突变概率特征通常被视为概率型数值特征，需要进行特征缩放或标准化处理，以便进行聚类分析

#特征缩放或标准化处理后，坐标值和突变概率特征都将被转化到相同的数值范围内，
#通常是[0,1]或[-1,1]的范围内。
#这样，它们就可以在聚类算法中进行比较和处理，而不会因为不同的度量单位和量纲而对聚类结果造成影响。

# 检查数据集中是否有NA和Infinity值
sum(is.na(data2$virusPercent))
sum(is.infinite(data2$virusPercent))
sum(is.na(scale_vector$V4))
sum(is.infinite(scale_vector$V4))

#转化方法1：同时标准化处理
#抽提突变频率
w<-data2[,c("virusPercent")]
#标准化xyz数据
w_norm<-scale(w)
# 将标准化处理过的坐标和突变频率组合为特征向量
scale_vector <- cbind(xyz_norm, w_norm)
#观察突变频率标准化前后的数据分布
p1 <- ggplot(data2, aes(x = virusPercent)) +
  geom_density(fill = "steelblue", color = "white") +
  xlab("virusPercent") + ylab("Density") + ggtitle("Original vP Density Distribution")
p1<-p1+xlim(0,0.02)

scale_vector<-data.frame(scale_vector)

p2 <- ggplot(scale_vector, aes(x = V4)) +
  geom_density(fill = "steelblue", color = "white") +
  xlab("virusPercent") + ylab("Density") + ggtitle("standardized vP Density Distribution")
p2<-p2+xlim(-0.2,-0.1)

# 将两张图合并为一行二列,形状，位置，峰度和偏度
grid.arrange(p1, p2, nrow = 1, ncol = 2)

#直方图
hist(data2[,4])
hist(scale_vector[,4])

# 对比前后突变频率的分布
ks.test(data2[, 4], scale_vector[, 4])  # 进行Kolmogorov-Smirnov检验，比较z列的分布
duplicated(data2[,4])
table(data2[,4])

#标准化处理可以将不同特征的数据转化到相同的尺度和分布范围内，
#消除不同特征之间的单位和量纲差异。
#因此，如果同时对坐标值和突变频率进行标准化处理，处理前后数据分布相似，
#那么可以认为已经将它们转化到了相同的尺度下。
#标准化处理将数据转化为均值为0，标准差为1的分布范围内，
#这意味着数据的取值范围已经被缩放到相同的尺度。
#如果处理前后数据分布相似，说明它们的方差和分布范围已经被调整到相同的尺度，因此可以认为它们已经被转化到了相同的尺度下。

#转化方法2：特征放缩
# 最小-最大缩放x，y，z坐标值
data2$x <- (data2$X - min(data2$X)) / (max(data2$X) - min(data2$X))
data2$y <- (data2$Y - min(data2$Y)) / (max(data2$Y) - min(data2$Y))
data2$z <- (data2$Z - min(data2$Z)) / (max(data2$Z) - min(data2$Z))
data2$mutation <- (data2$virusPercent - min(data2$virusPercent)) / (max(data2$virusPercent) - min(data2$virusPercent))

# 构建新的数据框只存放缩后的四列数据
new_data <- data.frame(x = data2$x, y = data2$y, z = data2$z, mutation = data2$mutation)

#验证特征放缩前后数据分布相似性

# 创建3D散点图
plot3d(new_data[,1], new_data[,2], new_data[,3], col = "blue", size = 2)

# 绘制散点图矩阵
pairs(data2[, c("X", "Y", "Z", "virusPercent")], main = "Scatterplot Matrix before")
pairs(new_data[,c("x","y","z","mutation")],main = "Scatterplot Matrix after")

# 比较缩放前后数据的直方图
par(mfrow = c(2, 4))   # 将图形区域划分为2行4列
hist(data2$virusPercent, main = "Mut Before Scale", xlab = "Mutation")
hist(data2$X, main = "X Before Scale", xlab = "X")
hist(data2$Y, main = "Y Before Scale", xlab = "Y")
hist(data2$Z, main = "Z Before Scale", xlab = "Z")
hist(new_data$mutation, main = "Mut After Scale", xlab = "Mutation")
hist(new_data$x, main = "X After Scale", xlab = "X")
hist(new_data$y, main = "Y After Scale", xlab = "Y")
hist(new_data$z, main = "Z After Scale", xlab = "Z")

# 比较缩放前后数据的密度图
par(mfrow = c(2, 4))   # 将图形区域划分为2行4列
plot(density(data2$virusPercent), main = "Mut Before Scale", xlab = "Mutation")
plot(density(data2$X), main = "X Before Scaling", xlab = "X")
plot(density(data2$Y), main = "Y Before Scaling", xlab = "Y")
plot(density(data2$Z), main = "Z Before Scaling", xlab = "Z")
plot(density(new_data$mutation), main = "Mut After Scale", xlab = "Mutation")
plot(density(new_data$x), main = "X After Scale", xlab = "X")
plot(density(new_data$y), main = "Y After Scale", xlab = "Y")
plot(density(new_data$z), main = "Z After Scale", xlab = "Z")

# 比较缩放前后数据的QQ图
par(mfrow = c(2, 4))   # 将图形区域划分为2行4列
qqnorm(data2$virusPercent, main = "Mut Before Scale")
qqline(data2$virusPercent)
qqnorm(data2$X, main = "X Before Scale")
qqline(data2$X)
qqnorm(data2$Y, main = "Y Before Scale")
qqline(data2$Y)
qqnorm(data2$Z, main = "Z Before Scale")
qqline(data2$Z)
qqnorm(new_data$mutation, main = "Mut After Scale")
qqline(new_data$mutation)
qqnorm(new_data$x, main = "X After Scale")
qqline(new_data$x)
qqnorm(new_data$y, main = "Y After Scale")
qqline(new_data$y)
qqnorm(new_data$z, main = "Z After Scale")
qqline(new_data$z)

#如果我现在如果发现标准化前后数据分布几乎相似的话，
#不能直接说明可以将数据标准化后去聚类得到的突变频率高的点可以直接对应到标准化前的点上。
#如果想要将标准化后的聚类结果对应到原始数据上，
#可以使用聚类算法的聚类中心或簇标记等信息，结合标准化前后的数据分布情况，来确定标准化后每个聚类簇所对应的原始数据点


# 进行 KMeans 聚类
set.seed(123)
kmeans_result <- kmeans(new_data, centers = 3, nstart = 20)
# 将聚类结果添加到数据框中
data2$cluster <- as.factor(kmeans_result$cluster)
# 可视化聚类结果&对聚类结果进一步分析
fviz_cluster(list(data = new_data, cluster = kmeans_result$cluster))
fviz_cluster(list(data = sub_data, cluster = kmeans_result$cluster))

library(dplyr)

# 统计每个聚类簇中数据点的数量
new_data$cluster <- as.factor(kmeans_result$cluster)
cluster_count <- new_data %>%
  group_by(cluster) %>%
  summarize(count=n())

# 计算每个聚类簇中数据点的比例
cluster_count$proportion <- cluster_count$count / nrow(new_data)
print(cluster_count)

library(ggplot2)

# 绘制聚类簇的频率分布直方图
ggplot(cluster_count, aes(x=cluster, y=count)) +
  geom_bar(stat="identity", fill="blue") +
  labs(x="Cluster", y="Count", title="Cluster Frequency Distribution")

# 绘制聚类簇的密度图
ggplot(new_data, aes(x=x, y=y, fill=cluster)) +
  geom_density_2d() +
  scale_fill_manual(values=c("red", "blue", "green", "yellow")) +
  labs(x="X", y="Y", fill="Cluster", title="Cluster Density Plot")

# 比较聚类结果和突变频率
cluster_freq <- data2 %>%
  group_by(cluster) %>%
  summarise(mean_freq = mean(virusPercent), count = n())

total_freq <- mean(data2$virusPercent)

print(cluster_freq)
print(total_freq)


#Permutation法(排名百分数假设检验+特征放缩处理)
# 计算原始数据的聚类结果
orig_clusters <- kmeans(new_data, 25)
data3<-data[1:972,c("X","Y","Z","virusPercent")]
orig_means <- tapply(new_data$mutation, orig_clusters$cluster, mean)
#这一步是为了计算原始数据中每个类的平均突变频率所处的百分位数。
#具体地，中包含原始数据中每个类突变频率的均值，那么 rank(orig_means) 返回原始数据中每个类的突变频率均值的排名。
#然后，通过除以 length(orig_means)，即原始数据被分成的类数，可以计算出每个类的均值的排名的比例。
#最后，乘以 100 并四舍五入到两位小数，即可得到排名。
orig_pcts <- round(rank(orig_means)/length(orig_means)*100, 2)

# 进行2000次突变概率重排并计算每次聚类结果
n_permutations <- 2000
perm_results <- matrix(NA, n_permutations, length(orig_means))
for (i in 1:n_permutations) {
  perm_data <- data3
  perm_data[,ncol(data3)] <- sample(data3[,ncol(data3)], replace=FALSE)
  mut_scale <- (perm_data[,ncol(data3)] - min(perm_data[,ncol(data3)])) / (max(perm_data[,ncol(data3)]) - min(perm_data[,ncol(data3)]))
  sta_scale <- data.frame(x = data2$x,y = data2$y,z = data2$z,mutation = mut_scale)
  perm_clusters <- kmeans(sta_scale, 25)
  perm_means <- tapply(sta_scale$mutation, perm_clusters$cluster, mean)
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
  #具体来说，对于原始数据中的每个聚类，该代码计算了重排结果中百分位数高于该聚类百分位数的比例（即p值）。
  #如果p值很小，说明原始聚类结果中该聚类的百分位数显著高于随机聚类结果，
  #从而支持该聚类是显著的。反之，如果p值很大，说明该聚类的百分位数可能与随机聚类结果中的百分位数相似，从而不支持该聚类是显著的
  p_values[i] <- mean(perm_results[,i] >= orig_pcts[i])
}
# 输出结果
cat("原始聚类结果的均值：", orig_means, "\n")
cat("原始聚类结果的百分位数：", orig_pcts, "\n")
cat("p值：", p_values, "\n")
cat("原始聚类结果中最高的平均突变频率：", max(orig_means), "\n")
cat("对应的p值：", p_values[which.max(orig_means)], "\n")

#每次重排都会对原始数据进行打乱，破坏原始数据之间的任何关联性，
#然后重新进行聚类，得到每个聚类的平均突变频率。
#重排得到的聚类结果是随机的，由于在每次重排中使用了相同的聚类算法和距离度量方法，
#因此重排结果可以用来估计原始聚类结果中每个聚类的可能性是否显著高于随机聚类结果。
#如果原始聚类结果中的每个聚类的百分位数显著高于随机聚类结果，那么说明该聚类在原始数据中可能具有显著的生物学意义

#我们想要看到的是1/3越多越好，2/3和3/3越少越好。所以当1000次中1/3的次数超过950时，即2/3和3/3的次数小于50时，可说明在该聚类相对排名的普适性。
#首先对原始数据进行聚类，计算每个聚类的平均突变频率，
#并将每个聚类的平均突变频率转换为对应的百分位数。
#然后，使用突变概率重排法重复执行聚类过程，得到1000个随机聚类结果，
#并计算每个聚类在随机聚类结果中的百分位数。
#最后，对于原始数据中的每个聚类，计算其在随机聚类结果中的排名比原始聚类结果高的比例，高意思就是平均突变率小。排名越高，对应类的平均突变频率越小。
#从而得到对应的p值。如果p值小于设定的显著性水平（通常是0.05），则拒绝原假设，原始聚类结果中该聚类不是随机出现的，否则接受原假设（即该聚类可能是由随机性产生的，没有显著性）

#如果原始数据聚类后的平均突变那一类的p值小于0.05的话，print出该类里的点
# 如果最高平均突变频率的类的p值小于0.05，则显示该类的点
if (p_values[which.max(orig_means)] < 0.05) {
  cat("原始聚类结果中平均突变频率最高的聚类中的点：\n")
  cat(paste("行号\tX\tY\tZ\tvirusPercent\n"))
  selected_rows <- which(orig_clusters$cluster == which.max(orig_means))
  for (i in selected_rows) {
    cat(paste(i, "\t", new_data[i, "x"], "\t", new_data[i, "y"], "\t", new_data[i, "z"], "\t", data[i, "virusPercent"], "\n"))
  }
  
  # 可视化所有点（用蓝色表示）和突出显示的点（用红色表示）
  significant_points <- which(orig_clusters$cluster == which.max(orig_means))  # 提取p值小于0.05的聚类中的所有数据点的索引
  print(significant_clusters)
  print(significant_points)
  # 提取p值小于0.05的聚类中的所有数据点的坐标
  significant_coordinates <- new_data[significant_points, c("x", "y", "z")]
  
  # 在三维点图中可视化所有数据点和p值小于0.05的聚类中的数据点
  fig <- plot_ly() %>%
    add_trace(data = new_data, x = ~x, y = ~y, z = ~z, type = "scatter3d", mode = "markers", name = "All Points") %>%
    add_trace(data = significant_coordinates, x = ~x, y = ~y, z = ~z, type = "scatter3d", mode = "markers", name = "Significant Clusters", marker = list(color = "red")) %>%
    layout(scene = list(xaxis = list(title = "X"), yaxis = list(title = "Y"), zaxis = list(title = "Z")))
  
  fig
}

#检验点的正确
cluster_data <- scale_vector[orig_clusters$cluster == which.max(orig_means), ]
print(cluster_data)




#Permutation法(排名百分数假设检验+标准化处理)
# 计算原始数据的聚类结果
orig_clusters <- kmeans(scale_vector, 20)
data3<-data[1:972,c("X","Y","Z","virusPercent")]
orig_means <- tapply(scale_vector$V4, orig_clusters$cluster, mean)
#这一步是为了计算原始数据中每个类的平均突变频率所处的百分位数。
#具体地，中包含原始数据中每个类突变频率的均值，那么 rank(orig_means) 返回原始数据中每个类的突变频率均值的排名。
#然后，通过除以 length(orig_means)，即原始数据被分成的类数，可以计算出每个类的均值的排名的比例。
#最后，乘以 100 并四舍五入到两位小数，即可得到排名。
orig_pcts <- round(rank(orig_means)/length(orig_means)*100, 2)

# 进行2000次突变概率重排并计算每次聚类结果
n_permutations <- 2000
perm_results <- matrix(NA, n_permutations, length(orig_means))
x_scale <- scale(data3[,1])
y_scale <- scale(data3[,2])
z_scale <- scale(data3[,3])
for (i in 1:n_permutations) {
  perm_data <- data3
  perm_data[,ncol(data3)] <- sample(data3[,ncol(data3)], replace=FALSE)
  mut_scale <- scale(perm_data[,ncol(data3)])
  sta_scale <- data.frame(x = x_scale,y = y_scale,z = z_scale,mutation = mut_scale)
  perm_clusters <- kmeans(sta_scale, 20)
  perm_means <- tapply(sta_scale$mutation, perm_clusters$cluster, mean)
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
  #具体来说，对于原始数据中的每个聚类，该代码计算了重排结果中百分位数高于该聚类百分位数的比例（即p值）。
  #如果p值很小，说明原始聚类结果中该聚类的百分位数显著高于随机聚类结果，
  #从而支持该聚类是显著的。反之，如果p值很大，说明该聚类的百分位数可能与随机聚类结果中的百分位数相似，从而不支持该聚类是显著的
  p_values[i] <- mean(perm_results[,i] >= orig_pcts[i])
}
# 输出结果
cat("原始聚类结果的均值：", orig_means, "\n")
cat("原始聚类结果的百分位数：", orig_pcts, "\n")
cat("p值：", p_values, "\n")
cat("原始聚类结果中最高的平均突变频率：", max(orig_means), "\n")
cat("对应的p值：", p_values[which.max(orig_means)], "\n")

#每次重排都会对原始数据进行打乱，破坏原始数据之间的任何关联性，
#然后重新进行聚类，得到每个聚类的平均突变频率。
#重排得到的聚类结果是随机的，由于在每次重排中使用了相同的聚类算法和距离度量方法，
#因此重排结果可以用来估计原始聚类结果中每个聚类的可能性是否显著高于随机聚类结果。
#如果原始聚类结果中的每个聚类的百分位数显著高于随机聚类结果，那么说明该聚类在原始数据中可能具有显著的生物学意义

#我们想要看到的是1/3越多越好，2/3和3/3越少越好。所以当1000次中1/3的次数超过950时，即2/3和3/3的次数小于50时，可说明在该聚类相对排名的普适性。
#首先对原始数据进行聚类，计算每个聚类的平均突变频率，
#并将每个聚类的平均突变频率转换为对应的百分位数。
#然后，使用突变概率重排法重复执行聚类过程，得到1000个随机聚类结果，
#并计算每个聚类在随机聚类结果中的百分位数。
#最后，对于原始数据中的每个聚类，计算其在随机聚类结果中的排名比原始聚类结果高的比例，高意思就是平均突变率小。排名越高，对应类的平均突变频率越小。
#从而得到对应的p值。如果p值小于设定的显著性水平（通常是0.05），则拒绝原假设，原始聚类结果中该聚类不是随机出现的，否则接受原假设（即该聚类可能是由随机性产生的，没有显著性）

#如果原始数据聚类后的平均突变那一类的p值小于0.05的话，print出该类里的点
# 如果最高平均突变频率的类的p值小于0.05，则显示该类的点
if (p_values[which.max(orig_means)] < 0.05) {
  cat("原始聚类结果中平均突变频率最高的聚类中的点：\n")
  cat(paste("行号\tX\tY\tZ\tvirusPercent\n"))
  selected_rows <- which(orig_clusters$cluster == which.max(orig_means))
  for (i in selected_rows) {
    cat(paste(i, "\t", scale_vector[i, "X"], "\t", scale_vector[i, "Y"], "\t", scale_vector[i, "Z"], "\t", scale_vector[i, "V4"], "\n"))
  }
  
  # 可视化所有点（用蓝色表示）和突出显示的点（用红色表示）
  significant_points <- which(orig_clusters$cluster == which.max(orig_means))  # 提取p值小于0.05的聚类中的所有数据点的索引
  print(significant_clusters)
  print(significant_points)
  # 提取p值小于0.05的聚类中的所有数据点的坐标
  significant_coordinates <- scale_vector[significant_points, c("X", "Y", "Z")]
  
  # 在三维点图中可视化所有数据点和p值小于0.05的聚类中的数据点
  fig <- plot_ly() %>%
    add_trace(data = scale_vector, x = ~X, y = ~Y, z = ~Z, type = "scatter3d", mode = "markers", name = "All Points") %>%
    add_trace(data = significant_coordinates, x = ~X, y = ~Y, z = ~Z, type = "scatter3d", mode = "markers", name = "Significant Clusters", marker = list(color = "red")) %>%
    layout(scene = list(xaxis = list(title = "X"), yaxis = list(title = "Y"), zaxis = list(title = "Z")))
  
  fig
}

#检验点的正确
cluster_data <- scale_vector[orig_clusters$cluster == which.max(orig_means), ]
print(cluster_data)










#Permutation法(均值假设检验+特征放缩处理)
# 计算原始聚类结果中每个聚类的均值
orig_clusters <- kmeans(new_data, 20)
orig_means <- tapply(new_data$mutation, orig_clusters$cluster, mean)
n_permutations <- 1000
data3<-data[1:972,c("X","Y","Z","virusPercent")]
# 确定突变高发区
high_freq_cluster <- which.max(orig_means)

# 初始化向量来保存每次重排后的最高平均值
max_means <- rep(NA, n_permutations)
p_v<-numeric(n_permutations)
p_v_1<-numeric(n_permutations)

# 进行一次1000次重排
for (i in 1:n_permutations) {
  perm_data <- data3
  perm_data[,ncol(data3)] <- sample(data3[,ncol(data3)], replace=FALSE)
  mut_scale <- (perm_data[,ncol(data3)] - min(perm_data[,ncol(data3)])) / (max(perm_data[,ncol(data3)]) - min(perm_data[,ncol(data3)]))
  sta_scale <- data.frame(x = data2$x,y = data2$y,z = data2$z,mutation = mut_scale)
  # 对重排后的数据进行聚类
  perm_clusters <- kmeans(sta_scale, 20)
  # 计算每个聚类的平均突变频率
  perm_means <- tapply(new_data$mutation, perm_clusters$cluster, mean)
  # 保存最高平均值
  max_means[i] <- max(perm_means)
}

# 计算原始聚类结果中突变高发区的平均值在排序后的位置
#我们使用 rank 变量计算原始聚类结果中突变高发区的平均值在重排结果中的排名。
#具体来说，我们将所有重排结果中的最高平均值进行排序，并计算突变高发区的原始平均值在排序后的位置。
#如果突变高发区的原始平均值在排序后的位置较靠前，
#那么它在重排结果中的表现就比较好，表明突变高发区是有统计学意义的；
#那么我们就可以拒绝原始聚类结果中突变高发区的平均值在随机重排中产生的最大值的排序位置上的假设，
#即突变高发区的平均值是随机出现的。
#反之，如果它在排序后的位置较靠后，那么它在重排结果中的表现就比较差，
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
    perm_data <- data3
    perm_data[,ncol(data3)] <- sample(data3[,ncol(data3)], replace=FALSE)
    mut_scale <- (perm_data[,ncol(data3)] - min(perm_data[,ncol(data3)])) / (max(perm_data[,ncol(data3)]) - min(perm_data[,ncol(data3)]))
    sta_scale <- data.frame(x = data2$x,y = data2$y,z = data2$z,mutation = mut_scale)
    
    # 对重排后的数据进行聚类
    perm_clusters <- kmeans(sta_scale, 20)
    
    # 计算每个聚类的平均突变频率
    perm_means <- tapply(new_data$mutation, perm_clusters$cluster, mean)
    
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
if(mean(p_v) < alpha){
  # 可视化所有点（用蓝色表示）和突出显示的点（用红色表示）
  significant_points <- which(orig_clusters$cluster == which.max(orig_means))  # 提取p值小于0.05的聚类中的所有数据点的索引
  print(significant_points)
  # 提取p值小于0.05的聚类中的所有数据点的坐标
  significant_coordinates <- new_data[significant_points, c("x", "y", "z")]
  
  # 在三维点图中可视化所有数据点和p值小于0.05的聚类中的数据点
  fig <- plot_ly() %>%
    add_trace(data = new_data, x = ~x, y = ~y, z = ~z, type = "scatter3d", mode = "markers", name = "All Points") %>%
    add_trace(data = significant_coordinates, x = ~x, y = ~y, z = ~z, type = "scatter3d", mode = "markers", name = "Significant Clusters", marker = list(color = "red")) %>%
    layout(scene = list(xaxis = list(title = "X"), yaxis = list(title = "Y"), zaxis = list(title = "Z")))
  
  fig
}




#Permutation法(均值假设检验+标准化处理)
# 计算原始聚类结果中每个聚类的均值
orig_clusters <- kmeans(scale_vector, 20)
orig_means <- tapply(scale_vector$V4, orig_clusters$cluster, mean)
n_permutations <- 1000
data3<-data[1:972,c("X","Y","Z","virusPercent")]
# 确定突变高发区
high_freq_cluster <- which.max(orig_means)

# 初始化向量来保存每次重排后的最高平均值
max_means <- rep(NA, n_permutations)
p_v<-numeric(n_permutations)
p_v_1<-numeric(n_permutations)
m <- 100
x_scale <- scale(data3[,1])
y_scale <- scale(data3[,2])
z_scale <- scale(data3[,3])
for (i in 1:m) {
  for (j in 1:n_permutations){
    perm_data <- data3
    perm_data[,ncol(data3)] <- sample(data3[,ncol(data3)], replace=FALSE)
    mut_scale <- scale(perm_data[,ncol(data3)])
    sta_scale <- data.frame(x = x_scale,y = y_scale,z = z_scale,mutation = mut_scale)
    
    # 对重排后的数据进行聚类
    perm_clusters <- kmeans(sta_scale, 20)
    
    # 计算每个聚类的平均突变频率
    perm_means <- tapply(scale_vector$V4, perm_clusters$cluster, mean)
    
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
if(mean(p_v) < alpha){
  # 可视化所有点（用蓝色表示）和突出显示的点（用红色表示）
  significant_points <- which(orig_clusters$cluster == which.max(orig_means))  # 提取p值小于0.05的聚类中的所有数据点的索引
  print(significant_points)
  # 提取p值小于0.05的聚类中的所有数据点的坐标
  significant_coordinates <- scale_vector[significant_points, c("X", "Y", "Z")]
  
  # 在三维点图中可视化所有数据点和p值小于0.05的聚类中的数据点
  fig <- plot_ly() %>%
    add_trace(data = scale_vector, x = ~X, y = ~Y, z = ~Z, type = "scatter3d", mode = "markers", name = "All Points") %>%
    add_trace(data = significant_coordinates, x = ~X, y = ~Y, z = ~Z, type = "scatter3d", mode = "markers", name = "Significant Clusters", marker = list(color = "red")) %>%
    layout(scene = list(xaxis = list(title = "X"), yaxis = list(title = "Y"), zaxis = list(title = "Z")))
  
  fig
}


# 计算原始聚类结果中每个聚类的均值
orig_clusters <- kmeans(new_data, 20)
n_permutations <- 1000
data3<-data[1:972,c("X","Y","Z","virusPercent")]

# 初始化向量来保存每次重排后的最高平均值
m <- 100
# 存储每个聚类的p-value
n_cluster<-length(unique(orig_clusters$cluster)) 
p_values <- matrix(NA, nrow=n_cluster, ncol=m)
p_values_1 <- matrix(NA, nrow = n_cluster, ncol = m)
# 存储p-value小于0.05的聚类的索引
significant_clusters <- c()

for (i in 1:n_cluster) {
  # 确定当前聚类的高频区域
  high_freq_cluster <- i
  
  # 计算当前聚类的平均突变频率
  orig_means <- tapply(new_data$mutation, orig_clusters$cluster, mean)[high_freq_cluster]
  
  # 初始化向量来保存每次重排后的最高平均值
  max_means <- rep(NA, n_permutations)
  
  for (u in 1:m){
    for (j in 1:n_permutations) {
      perm_data <- data3
      perm_data[,ncol(data3)] <- sample(data3[,ncol(data3)], replace=FALSE)
      mut_scale <- (perm_data[,ncol(data3)] - min(perm_data[,ncol(data3)])) / (max(perm_data[,ncol(data3)]) - min(perm_data[,ncol(data3)]))
      sta_scale <- data.frame(x = data2$x,y = data2$y,z = data2$z,mutation = mut_scale)
      # 对重排后的数据进行聚类
      perm_clusters <- kmeans(sta_scale, n_cluster)
      
      # 计算每个聚类的平均突变频率
      perm_means <- tapply(new_data$mutation, perm_clusters$cluster, mean)
      
      # 保存最高平均值
      max_means[j] <- max(perm_means)
    }
    
    # 计算p-value
    rank <- sum(max_means >= orig_means) / n_permutations
    p_value <- min(rank, 1-rank)
    p_value2 <- 2 * (1 - pnorm(abs(rank - 0.5) * sqrt(n_permutations / 12)))
    p_values[i,u] <- p_value
    p_values_1[i,u]<- p_value2
  }
}

for (i in 1:n_cluster){
  # 如果p-value小于0.05，则将当前聚类的索引添加到significant_clusters向量中
  if (mean(p_value[i,]) < 0.05 && mean(p_values_1[i,] < 0.05)) {
    significant_clusters <- c(significant_clusters, i)
  }
}

# 初始化三维散点图
fig <- plot_ly() %>%
  add_trace(data = new_data, x = ~x, y = ~y, z = ~z, type = "scatter3d", mode = "markers", name = "All Points")

# 循环遍历每个聚类
colors <- heat.colors(n_cluster)
for (i in 1:n_cluster){
  # 如果p-value小于0.05，则将当前聚类的索引添加到significant_clusters向量中
  if (mean(p_value[i,]) < 0.05 && mean(p_values_1[i,] < 0.05)) {
    # 提取属于当前聚类的数据点的索引
    cluster_points <- which(orig_clusters$cluster == i)
    # 提取属于当前聚类的数据点的坐标
    cluster_coordinates <- new_data[cluster_points, c("x", "y", "z")]
    # 为聚类指定颜色
    color <- color[i]
    # 在三维散点图中添加聚类
    fig <- fig %>% add_trace(data = cluster_coordinates, x = ~x, y = ~y, z = ~z, type = "scatter3d", mode = "markers", name = paste("Cluster ", i), marker = list(color = color))
  }
}

# 设定散点图的布局
fig <- fig %>% layout(scene = list(xaxis = list(title = "X"), yaxis = list(title = "Y"), zaxis = list(title = "Z")))

# 显示散点图
fig