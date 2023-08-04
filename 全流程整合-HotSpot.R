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
library(factoextra)
library(fpc)
library(NbClust)
library(readxl)
library(plot3D)


# Read Data ---------------------------------------------------------------

# Load data with spatial coordinates and mutation rates
#Pipe
data<-as.data.frame(read.csv("data1.csv"))
#Xiao
setwd("/Users/ryan/Documents/GitHub/Shenzhen Bay lab/spike_protein_3structures_markerd_with_variants_from_cncb")
data <- read_excel("6vxx_variants.xls")


# Cluster Pipe -----------------------------------------------------------------

data<-na.omit(data) #去除na值
data_matrix<-as.matrix(data)
data2<-data[1:972,c("X","Y","Z","virusPercent")]
log_virus_percent <- log(data2$virusPercent)


# 生成随机的三维坐标数据
#Pipe
data <- read.csv("data1.csv")
#Xiao
data <- read_excel("6vxx_variants.xls")


sub_data<-data[1:972,c("X","Y","Z","virusPercent")]
sub_data_1<-data[1:972,c("X","Y","Z")]
# 计算每个维度的基本统计信息
summary(sub_data)

##############
# 绘制散点图
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

###############

#计算每个维度的均值和标准差
apply(sub_data,2,mean)
apply(sub_data,2,sd)


# 进行双样本t检验
t.test(sub_data[, 1], sub_data[, 2])
t.test(sub_data[, 1], sub_data[, 3])
t.test(sub_data[, 2], sub_data[, 3])
#H0：两个独立样本的均值相等。
#H1：两个独立样本的均值不相等。
#如果p值很小，说明两个维度之间的均值和标准差之间存在显著差异，
#需要进行标准化处理等进一步操作

###########

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
