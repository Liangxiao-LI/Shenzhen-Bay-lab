
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


