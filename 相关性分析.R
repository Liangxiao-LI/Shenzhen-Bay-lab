#探究1000个点的空间相关性，可以使用空间统计分析方法，
#例如空间自相关分析、空间克里金插值等等。
#这些方法可帮助识别地理空间数据中的空间相关性和模式，从而更好地理解和分析数据。
#其中，空间自相关分析是一种常用的空间统计方法，它可以用来评估空间数据中的局部和全局空间相关性。
#该方法通过计算每个点与其邻近点之间的空间距离和属性值之间的关系，来确定空间相关性。
#其中，局部空间自相关分析可以帮助您确定空间数据中的局部空间相关性和空间聚集，
#而全局空间自相关分析可以帮助您确定整个空间数据中的全局空间相关性和空间分布模式。

#此外，空间克里金插值是一种用于空间数据插值的方法，
#它可以帮助您估计空间数据在未观测位置的值。
#该方法通过使用邻近点之间的空间距离和属性值之间的关系来插值未观测点的值，
#从而估计整个空间数据的空间分布模式。

#使用突变频率作为评判点于点之间空间相关性的依据。
#由于突变频率是每个点周围的突变数量与总细胞数量之比，
#因此可以反映每个点周围的突变情况。通过比较不同点的突变频率，
#可以评估它们之间的空间相关性和空间分布模式。

#空间相关性可以体现在以下几个方面：

#1.空间聚集：如果这1000个点中的一部分点在空间上彼此靠近，形成了一个空间聚集的模式，那么这些点之间就存在着空间相关性。
#2.空间分散：相反，如果这1000个点的分布相对均匀、没有明显的空间聚集模式，那么这些点之间就不存在明显的空间相关性。
#3.空间自相关：如果这1000个点的属性值在空间上彼此相关，即相邻的点具有相似的属性值，那么这些点之间就存在着空间自相关性。
#可以通过计算每个点周围其他点的属性值和距离之间的关系来评估空间自相关性。
#4.空间异质性：如果这1000个点的属性值在不同的空间位置上具有不同的变化模式，
#即空间异质性，那么这些点之间也存在着空间相关性。
#可以通过计算每个点周围其他点的属性值和距离之间的关系，并根据属性值的变化模式来评估空间异质性。

#探究思路：手上已有的数据是1000个点的xyz坐标和突变频率，
#先只探究所有点的空间位置之间是否存在空间相关性，
#然后再单独探究所有点的突变频率之间的分布，
#最后探究空间位置之间存在的空间相关性是否会对突变频率的分布产生影响，
#即位置时间的空间相关性是否和突变频率的分布之间存在某种联系


#探究所有点的空间位置之间是否存在空间相关性。
#可以使用全局空间自相关分析方法(例如Moran's I指数)来评估所有点的空间相关性。
#可帮助了解整个区域的空间分布模式和空间聚集情况。

#探究所有点的突变频率之间的分布。
#可以使用直方图、密度图或箱线图等方法来描述突变频率的分布情况。
#可帮助了解突变频率的中心趋势、分散程度和异常值情况。

#探究空间位置之间存在的空间相关性是否会对突变频率的分布产生影响。
#可以使用空间回归分析方法（例如空间误差模型或空间Durbin模型）来评估空间位置和突变频率之间的关系。
#具体来说，可以将突变频率作为响应变量，距离和其他可能影响突变频率的因素（例如环境因素）作为解释变量，
#然后使用空间回归模型进行建模和分析，可帮助确定空间依赖性和空间自相关性，并评估各个解释变量对突变频率的影响程度。

library(spdep)
library(sp)
# 导入数据
data <- read.csv("data1.csv")
data2<-data[1:972, c("X", "Y", "Z","virusPercent")]

set.seed(123)
# 提取 x、y、z 列
coords <- data2[, c("X", "Y", "Z")]

# 计算欧几里得距离
#欧式距离：欧式距离是最常见的距离度量方法之一，
#它计算两点之间的直线距离。欧式距离越小，表示两点之间越接近，
#反之则越远离。欧式距离可以告诉我们点之间的相对距离和空间分布范围，
#以及在某些情况下可能存在的空间集聚或分散现象。
dist_euclidean <- dist(coords)  # 默认使用欧几里得距离
dist_euclidean_m<-as.matrix(dist_euclidean)
dist_euclidean_m[1:10, 1:10]  # 查看前 10 行 10 列的距离矩阵

# 进行空间聚类分析
hc_euclidean <- hclust(dist_euclidean)  # 使用层次聚类算法
plot(hc_euclidean)  # 绘制聚类树图

# 绘制三维散点图
library(plotly)
plot_ly(data2, x = ~X, y = ~Y, z = ~Z, type = "scatter3d", mode = "markers")

# 计算曼哈顿距离
#曼哈顿距离计算两点之间的曼哈顿距离，即两点在 x、y、z 等方向上的距离之和。
#曼哈顿距离可以告诉我们点之间的空间分布规律
dist_manhattan <- as.matrix(dist(coords, method = "manhattan"))
dist_manhattan[1:10, 1:10]  # 查看前 10 行 10 列的距离矩阵

# 绘制散点图
#曼哈顿距离的散点图可以用来分析点之间的相似性和空间分布规律。
#点的颜色表示点之间的曼哈顿距离，颜色越深表示距离越远，颜色越浅表示距离越近。
#通过观察散点图，我们可以得到以下信息：

#点的颜色趋于相似，表示这些点在空间上更加接近。
#点的颜色趋于不相似，表示这些点在空间上更加分散。
#颜色分布越均匀，表示空间分布越平均。
#颜色分布越集中，表示空间分布越集聚。
plot(data2, col = dist_manhattan, pch = 20, main = "Manhattan Distance Scatterplot")

# 绘制热力图
#曼哈顿距离的热力图可以用来更直观地分析点之间的相似性和空间分布规律。
#颜色越深表示相应的距离越远，颜色越浅表示距离越近。
#通过观察热力图，我们可以得到以下信息：

#色块越小，表示相应的点之间距离越近。
#色块越大，表示相应的点之间距离越远。
#颜色分布越均匀，表示空间分布越平均。
#颜色分布越集中，表示空间分布越集聚。
library(ggplot2)
ggplot(as.data.frame(as.table(dist_manhattan)), aes(Var1, Var2, fill = Freq)) +
  geom_tile() + scale_fill_gradient(low = "white", high = "red") +
  labs(title = "Manhattan Distance Heatmap", x = "", y = "")



# 构建空间点对象
#pts <- SpatialPoints(coords)

# 定义插值参数
#idw <- gstat::idw(value ~ 1, pts, nmax = 10, idp = 2)


#探究每个点周围存在的其他突变点的个数和该点突变频率之间的关系

#构建空间点对象
library(sp)
coords <- cbind(data2$X, data2$Y, data2$Z)
pts <- SpatialPoints(coords)

#计算点之间的距离
library(spdep)
dist_mat <- nbdists(pts)

# 定义半径范围
radius_range <- seq(0.1, 1, by = 0.1) # 从0.1到1，每次增加0.1

# 计算每个半径下的相邻点个数
neighbours_list <- lapply(radius_range, function(radius) {
  neighbours <- apply(dist_mat, 1, function(x) sum(x <= radius))
  return(neighbours)
})

# 绘制相邻点个数随半径变化的折线图
library(ggplot2)
# 整理数据
df <- data.frame(radius = rep(radius_range, each = nrow(data2)),
                 neighbours = unlist(neighbours_list),
                 freq = rep(data2$virusPercent, length(radius_range)))

# 绘制图形
library(ggplot2)
ggplot(df, aes(x = radius, y = neighbours, color = freq)) +
  geom_line() +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  labs(x = "Radius", y = "Neighbours", color = "Mutation Frequency")


# 绘制相邻点个数随半径变化的折线图
library(ggplot2)
df <- data.frame(radius = rep(radius_range, each = nrow(data2)),
                 neighbours = unlist(neighbours_list),
                 freq = rep(data2$virusPercent, times = length(radius_range)))
ggplot(df, aes(x = radius, y = neighbours, group = freq)) +
  geom_line(aes(color = freq)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  labs(x = "Radius", y = "Neighbours", color = "Mutation Frequency")

#确定每个点周围存在的其他突变点的个数
radius<-10
neighbours_1 <- apply(dist_mat, 1, function(x) sum(x <= radius))

#探究点周围突变个数和突变频率之间的关系
plot(neighbours_1, data2$virusPercent)
