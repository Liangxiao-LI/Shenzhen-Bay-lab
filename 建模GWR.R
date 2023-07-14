
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
data3<-as.matrix(data2)
data_subset1<-data[1:972,6]
data_subset2<-data[1:972,7]
data_subset3<-data[1:972,8]
data_subset4<-data[1:972,10]
data_subset5<-as.matrix(data_subset4)
log_virus_percent <- log(data2$virusPercent)

set.seed(123)
data4<-data.frame(x = data_subset1,y = data_subset2,z = data_subset3,w = data_subset4)
q1 <- quantile(data_subset4, 0.25)
q2 <- quantile(data_subset4, 0.5)
q3 <- quantile(data_subset4, 0.75)
data_smaller_Q1 <- data4[data4$w < q1,]
data_between_Q1_Q2 <- data4[data4$w > q1 & data$w < q2, ]
data_between_Q2_Q3 <- data4[data4$w > q2 & data$w < q3, ]
data_greater_Q4 <- data4[data4$w > q3, ]

data4$category <- ifelse(data4$w < q1, 1,
                        ifelse(data4$w < q2, 2,
                               ifelse(data4$w < q3, 3, 4)))
colors <- c("red", "green", "blue", "purple")
col_data <- colors[data4$category]

plot3d(data4$x, data4$y, data4$z, col = col_data)
legend3d("topright", legend = c("Category 1", "Category 2", "Category 3", "Category 4"), 
         col = colors, pch = 16, cex = 1.5)

hist(data_subset4)

scatterplot3d(data_subset1, data_subset2, data_subset3, color = "blue", pch = 16, main = "3D Scatter Plot")

kmeans_result <- kmeans(cbind(data_subset4), centers = 3)
scatterplot3d(data_subset1, data_subset2, data_subset3, color = kmeans_result$cluster, pch = 16, 
              main = "3D Scatter Plot with Cluster Colors")
plot3d(data_subset1, data_subset2, data_subset3, type = "s", size = 2, col = kmeans_result$cluster, xlab = "X", ylab = "Y", zlab = "Z")
rglwidget(elementId = "rglplot", spin = TRUE)

#Kmean聚类和散点图可视化hotspot
coords <- data2 %>%
  select(X, Y, Z)
mutation_freq <- log_virus_percent
num_clusters <- 5
kmeans_model <- kmeans(mutation_freq, centers = num_clusters)
cluster_labels <- kmeans_model$cluster
clustered_data <- data.frame(coords, cluster_labels)
plot_ly(clustered_data, x = ~X, y = ~Y, z = ~Z, color = ~factor(cluster_labels),
        type = "scatter3d", mode = "markers", marker = list(size = 3)) %>%
  layout(scene = list(xaxis = list(title = "X"),
                      yaxis = list(title = "Y"),
                      zaxis = list(title = "Z")))
cluster1_data <- data2[kmeans_model$cluster == 1, ]
cluster2_data <- data2[kmeans_model$cluster == 2, ]
cluster3_data <- data2[kmeans_model$cluster == 3, ]
cluster4_data <- data2[kmeans_model$cluster == 4, ]
cluster5_data <- data2[kmeans_model$cluster == 5, ]
d1<-mean(cluster1_data$virusPercent)
d2<-mean(cluster2_data$virusPercent)
d3<-mean(cluster3_data$virusPercent)
d4<-mean(cluster4_data$virusPercent)
d5<-mean(cluster5_data$virusPercent)

#DBSCAN模型聚类
#代码1
set.seed(123)
iris<-data2$virusPercent
iris<-as.matrix(iris)
dbscan_model <- dbscan(iris, eps = 0.5, minPts = 5)
data_clustered <- data.frame(data2, cluster = dbscan_model$cluster)
scatter_data <- data.frame(
  x = data2$X,
  y = data2$Y,
  z = data2$Z,
  cluster = factor(data_clustered$cluster)
)
plot_ly(scatter_data, x = ~x, y = ~y, z = ~z, color = ~cluster, type = "scatter3d", mode = "markers")%>%
  layout(scene = list(xaxis = list(title = "X"), yaxis = list(title = "Y"), zaxis = list(title = "Z")))

#代码2
coords <- data2 %>%
  select(X, Y, Z)
mutation_freq <- log_virus_percent
mutation_freq_1<-as.matrix(mutation_freq)
dbscan_model <- dbscan(mutation_freq_1, eps = 0.05, minPts = 5)
cluster_labels <- dbscan_model$cluster
clustered_data <- data.frame(coords, cluster_labels)
plot_ly(clustered_data, x = ~X, y = ~Y, z = ~Z, color = ~factor(cluster_labels),
        type = "scatter3d", mode = "markers", marker = list(size = 3)) %>%
  layout(scene = list(xaxis = list(title = "X"),
                      yaxis = list(title = "Y"),
                      zaxis = list(title = "Z")))
cluster6_data <- data2[dbscan_model$cluster == 1, ]
cluster7_data <- data2[dbscan_model$cluster == 2, ]
cluster8_data <- data2[dbscan_model$cluster == 3, ]
cluster9_data <- data2[dbscan_model$cluster == 4, ]
cluster10_data <- data2[dbscan_model$cluster == 5, ]
cluster11_data <- data2[dbscan_model$cluster == 6, ]
cluster12_data <- data2[dbscan_model$cluster == 7, ]
cluster13_data <- data2[dbscan_model$cluster == 8, ]
cluster14_data <- data2[dbscan_model$cluster == 9, ]
cluster15_data <- data2[dbscan_model$cluster == 10, ]
d6<-mean(cluster6_data$virusPercent)
d7<-mean(cluster7_data$virusPercent)
d8<-mean(cluster8_data$virusPercent)
d9<-mean(cluster9_data$virusPercent)
d10<-mean(cluster10_data$virusPercent)
d11<-mean(cluster11_data$virusPercent)
d12<-mean(cluster12_data$virusPercent)
d13<-mean(cluster13_data$virusPercent)
d14<-mean(cluster14_data$virusPercent)
d15<-mean(cluster15_data$virusPercent)

#高斯混合模型
xyz<-data2[,c("X","Y","Z")]
model<-Mclust(mutation_freq)
clusters <- model$classification
s3d <- scatterplot3d(xyz, type = "p", color = clusters, 
                     pch = 16, main = "Mutation frequency clusters", 
                     xlab = "X", ylab = "Y", zlab = "Z")
plot3d(xyz, type = "s", size = 1, col = clusters)

#多元聚类模型
#标准化xyz数据
xyz_norm<-scale(xyz)
# 合并 xyz 坐标和突变频率数据
data2_cbind <- cbind(xyz_norm, mutation_freq)
# 计算距离矩阵
dist_matrix <- dist(data2_cbind)
# 进行多元聚类分析
model_mul <- pam(dist_matrix, k = 3)
# 绘制三维散点图
plot3d(xyz, col = heat.colors(nlevels(factor(model_mul$clustering)))[model_mul$clustering], type = "s", size = 1)
legend3d("topright", legend = levels(factor(model_mul$clustering)), 
         col = heat.colors(nlevels(factor(model_mul$clustering))), pch = 16, cex = 1.2)
data_mul_1<-data2[model_mul$clustering ==1,]
data_mul_2<-data2[model_mul$clustering ==2,]
data_mul_3<-data2[model_mul$clustering ==3,]
mul_1<-mean(data_mul_1$virusPercent)
mul_2<-mean(data_mul_2$virusPercent)
mul_3<-mean(data_mul_3$virusPercent)

#PCA模型
# 进行主成分分析
pca <- prcomp(xyz_norm)
# 取出前两个主成分
pca_data <- pca$x[, 1:2]
# 合并主成分和突变频率数据
data2_cbind_1 <- cbind(pca_data, mutation_freq)
# 定义颜色向量
colors <- rainbow(100)
# 计算距离矩阵
dist_matrix_1 <- dist(data2_cbind_1)
model_mul_1 <- pam(dist_matrix_1, k = 3)
# 绘制三维散点图
plot3d(xyz, col = heat.colors(nlevels(factor(model_mul_1$clustering)))[model_mul_1$clustering], type = "s", size = 1)
legend3d("topright", legend = levels(factor(model_mul_1$clustering)), 
         col = heat.colors(nlevels(factor(model_mul_1$clustering))), pch = 16, cex = 1.2)
data_pca_1<-data2[model_mul_1$clustering ==1,]
data_pca_2<-data2[model_mul_1$clustering ==2,]
data_pca_3<-data2[model_mul_1$clustering ==3,]
pca_1<-mean(data_pca_1$virusPercent)
pca_2<-mean(data_pca_2$virusPercent)
pca_3<-mean(data_pca_3$virusPercent)

#分区域假设检验
#方法一
result <- pam(dist(xyz_norm), k = 5)
# 取出聚类结果中心点
centers <- result$medoids
# 取出所有点中心点
global_center <- colMeans(xyz_norm)
# 计算每个区域的中心点
region_centers <- aggregate(xyz_norm, by = list(result$clustering), FUN = "mean")[,-1]
# 计算每个区域的中心点到整体中心点的距离
dist_to_global_center <- colSds(t(region_centers)) * sqrt(nrow(region_centers))
# 进行假设检验
t.test(dist_to_global_center, mu = 0)
#我们拒绝掉了H0，说明这种方法对点进行区域划分并不均匀
cluster_labels_test_1 <- as.factor(result$cluster)
plot_ly(data2, x = ~X, y = ~Y, z = ~Z, color = ~factor(result$cluster), size = 1, type = "scatter3d") %>%
  layout(title = "Clustered Data", scene = list(camera = list(eye = list(x = -1.8, y = -1.8, z = 0.5))))
mut_rate_1 <- aggregate(log_virus_percent, by = list(result$cluster), FUN = mean)
colnames(mut_rate_1) <- c("cluster", "avg_mut_rate")
combined_data_11 <- cbind(mut_rate_1, dist_to_global_center)
ggplot(combined_data_11, aes(x = dist_to_global_center, y = avg_mut_rate)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Distance to Center") +
  ylab("Average Mutation Rate") +
  ggtitle("Relationship between Distance and Mutation Rate")

#方法2 不标准化坐标
# 计算所有点的中心点
center_point <- colMeans(data2[, 1:3])
# 将所有点划分成几个区域
result_1 <- kmeans(data2[, 1:3], 200, nstart = 20)
# 计算每个区域的中心点
region_centers_1 <- result_1$centers
# 计算每个区域的中心点到整体中心点的距离
dist_to_center <- apply(region_centers_1, 1, function(x) dist(rbind(x, center_point)))
# 使用假设检验检验每个区域的中心点到整体中心点的距离是否相等
t.test(dist_to_center, mu = 0)
cluster_labels_test <- as.factor(result_1$cluster)
plot_ly(data2, x = ~X, y = ~Y, z = ~Z, color = ~factor(result_1$cluster), size = 1, type = "scatter3d") %>%
  layout(title = "Clustered Data", scene = list(camera = list(eye = list(x = -1.8, y = -1.8, z = 0.5))))
mut_rate <- aggregate(log_virus_percent, by = list(result_1$cluster), FUN = mean)
colnames(mut_rate) <- c("cluster", "avg_mut_rate")
combined_data_10 <- cbind(mut_rate, dist_to_center)
ggplot(combined_data_10, aes(x = dist_to_center, y = avg_mut_rate)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Distance to Center") +
  ylab("Average Mutation Rate") +
  ggtitle("Relationship between Distance and Mutation Rate")

#现在有一个问题，比如说我现在得到了距离和突变频率的关系，但在三维空间里，
#到一个点的距离一样的点可以形成一个球面，
#我不能说这球面上所有点的突变概率都是一样因为大概率他们是不一样的
#所以接下来进行等距度分析，更细致地探究距离和突变频率之间的关系
# 创建距离和突变率的数据框
# 计算所有点的中心点
#center_point <- colMeans(data2[, 1:3])
# 计算每个点到中心点的距离
#dist_to_center <- apply(data2[, 1:3], 1, function(x) sqrt(sum((x - center_point)^2)))
# 计算每个点的突变率
#mut_rate <- log_virus_percent
# 创建距离和突变率的数据框
#dismut_data <- data.frame(dist_to_center, mut_rate)
# 使用等距度分析估计距离上的突变率分布
#iso_fit <- gstat::idw(mut_rate ~ dist_to_center, dismut_data)
#dis_mut_grid <- data.frame(dist_to_center = seq(min(dist_to_center), max(dist_to_center), length.out = 100))
#mut_rate_grid <- predict(iso_fit, newdata = dis_mut_grid)
#iso_data <- data.frame(dist_to_center = dis_mut_grid$dist_to_center, mut_rate_grid)

#三维空间点转化为特征向量聚类
X1 <- as.matrix(data2[, 1:4])
# 将坐标和突变频率组合为特征向量
feature_vectors <- apply(X1, 1, function(x) c(x[1:3], x[4]))
# 对特征向量进行聚类分析
cluster_labels <- kmeans(feature_vectors, centers = 3)
# 将聚类结果可视化为三维散点图
colors <- c("red", "blue", "green") # 每个聚类的颜色
points3d(data2$X, data2$Y, data2$Z, col = colors[cluster_labels$cluster], size = 5)

#特征向量法2
# 将坐标进行标准化处理
scaled_coord <- data2 %>%
  select(X, Y, Z) %>%
  as.matrix() %>%
  rescale()
# 将坐标和突变频率组合为特征向量
feature_vector <- cbind(scaled_coord, data2$virusPercent)
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
# 比较聚类结果和已知的突变高发区
known_hotspot <- data2 %>%
  filter(X > 0.5 & Y > 0.5 & Z > 0.5)

print(known_hotspot)
# 比较不同聚类算法和参数的效果
fviz_nbclust(feature_vector, kmeans, method = "silhouette")
# 输出聚类效果评估指标
#silhouette <- silhouette(kmeans_result$cluster, feature_vector)
#calinski_harabasz <- calinski_harabasz(feature_vector, kmeans_result$cluster)
#davies_bouldin <- davies_bouldin(feature_vector, kmeans_result$cluster)

#print(silhouette)
#print(calinski_harabasz)
#print(davies_bouldin)


#空间插值聚类
# 创建空间数据框
coordinates(data2) <- c("X", "Y", "Z")

# 测试不同的半方差函数类型
vgm_sph <- vgm(psill = 1, model = "Sph", range = 100, nugget = 0.1)
vgm_exp <- vgm(psill = 1, model = "Exp", range = 100, nugget = 0.1)
vgm_gau <- vgm(psill = 1, model = "Gau", range = 100, nugget = 0.1)

# 创建空间插值模型(有误)
#fit_sph <- fit.variogram(variogram(data2$virusPercent ~ 1, data = data2), model = vgm_sph)
#fit_exp <- fit.variogram(variogram(data2$virusPercent ~ 1, data = data2), model = vgm_exp)
fit_gau <- fit.variogram(variogram(data2$virusPercent ~ 1, data = data2), model = vgm_gau)

# 查看每个半方差函数类型的 psill 值
print(c("Spherical" = vgm_sph$psill, "Exponential" = vgm_exp$psill, "Gaussian" = vgm_gau$psill))

# 进行空间插值
grd <- expand.grid(x = seq(min(data2$X), max(data2$X), by = 0.5),
                   y = seq(min(data2$Y), max(data2$Y), by = 0.5),
                   z = seq(min(data2$Z), max(data2$Z), by = 0.5))
coordinates(grd) <- c("x", "y", "z")
pred_sph <- predict(fit_sph, grd)
pred_exp <- predict(fit_exp, asgrd)
pred_gau <- predict(fit_gau, grd)

# 可视化插值结果
par(mfrow = c(1, 3))
slice_sph <- pred_sph[z == median(data2$Z), ]
image(slice_sph[, c("x", "y")], col = terrain.colors(100), xlab = "x", ylab = "y", main = "Spherical")
slice_exp <- pred_exp[z == median(data2$Z), ]
image(slice_exp[, c("x", "y")], col = terrain.colors(100), xlab = "x", ylab = "y", main = "Exponential")
slice_gau <- pred_gau[z == median(data2$Z), ]
image(slice_gau[, c("x", "y")], col = terrain.colors(100), xlab = "x", ylab = "y", main = "Gaussian")





#活性中心聚类
#以突变频率最小的5个点作为活性中心
smallest_freq <- data2[order(data2$virusPercent)[1:5], ]
smallest_coords <- smallest_freq[, c("X", "Y", "Z")]
print(smallest_freq)
print(smallest_coords)
scatter_data <- data.frame(
  x = data2$X,
  y = data2$Y,
  z = data2$Z,
  marker_color = ifelse(
    rownames(data2) %in% rownames(smallest_freq),
    "red",
    "blue"
  )
)
plot_ly(scatter_data, x = ~x, y = ~y, z = ~z, color = ~marker_color, type = "scatter3d", mode = "markers") %>%
  layout(scene = list(xaxis = list(title = "X"), yaxis = list(title = "Y"), zaxis = list(title = "Z")))
#从十个点中选点使得其能均匀分布在所有点中
# Calculate the range of each coordinate
range_x <- max(smallest_coords$X) - min(smallest_coords$X)
range_y <- max(smallest_coords$Y) - min(smallest_coords$Y)
range_z <- max(smallest_coords$Z) - min(smallest_coords$Z)

# Calculate the average distance between points in each dimension
avg_dist_x <- range_x / (length(smallest_coords$X) - 1)
avg_dist_y <- range_y / (length(smallest_coords$Y) - 1)
avg_dist_z <- range_z / (length(smallest_coords$Z) - 1)

# Calculate the maximum number of points that can be selected
max_points <- min(length(smallest_coords$X), floor(range_x / avg_dist_x),
                  floor(range_y / avg_dist_y), floor(range_z / avg_dist_z))

# Select points that are evenly spaced in each dimension
selected_coords <- smallest_coords%>%
  slice(as.integer(seq(1, nrow(smallest_coords), length.out = max_points)))

# Print the selected coordinates
print(selected_coords)


# Create SpatialPointsDataFrame object with existing coordinates
coordinates(data2) <- c("X", "Y", "Z")

# Create plot of the data
plot(data2$virusPercent ~ data2$X + data2$Y + data2$Z, col = "blue", pch = 19, main = "Mutation Frequency by 3D Spatial Coordinates")

#Spatial autoregression analysis
# Create spatial weights matrix
nn <- knn2nb(knearneigh(coordinates(data2), k = 4), sym = TRUE)
w <- nb2listw(nn)

# Calculate Moran's I index for mutation rates
moran_index <- moran.test(data2$virusPercent, w)

# Print results
cat("Moran's I index:", moran_index$estimate, "\n")
cat("p-value:", moran_index$p.value, "\n")



#CLUSTER THE DATA

#K-mean cluster
# Perform k-means clustering with k = 3
kmeans_result <- kmeans(data_subset4, centers = 3)
# Print cluster assignments
print(kmeans_result$cluster)
# Create cluster plot
clusplot(data_subset4, kmeans_result$cluster, color = TRUE, shade = TRUE)
cluster1_data <- data2[kmeans_result$cluster == 1, ]
cluster2_data <- data2[kmeans_result$cluster == 2, ]
cluster3_data <- data2[kmeans_result$cluster == 3, ]

#Hierarchical cluster
# Perform hierarchical clustering with complete linkage
dist_matrix <- dist(data_subset4)
hclust_result <- hclust(dist_matrix, method = "complete")
# Print cluster assignments
print(cutree(hclust_result, k = 3))
# Create dendrogram plot
plot(hclust_result, hang = -1, main = "Dendrogram of Mutation Frequencies", 
     sub = "Complete Linkage", xlab = "Samples", ylab = "Distance")


