library(readxl)

setwd("/Users/ryan/Documents/GitHub/Shenzhen Bay lab/spike_protein_3structures_markerd_with_variants_from_cncb")

data <- read_excel("6vxx_variants.xls")

library(cluster)

library(dbscan)

data_mutation<-data[1:972,10]

k <- 3  # 聚类数目
kmeans_result <- kmeans(data_mutation, centers = k)

# 获取聚类结果
cluster_labels <- kmeans_result$cluster

#data_mutation$cluster <- as.factor(cluster_labels)

# 可视化聚类结果

clusplot(data_mutation, cluster_labels, color = TRUE)

#合成数据

merged_data <- cbind(data_mutation, cluster_labels)

high_mutation_data <- merged_data[grep("1|2", merged_data$cluster_labels), ]

dbscan_result <- dbscan(high_mutation_data, eps = 0.5, minPts = 2)

plot(dbscan_result, high_mutation_data)

plot(merged_data$virusPercent ~ merged_data$X + data2$Y + data2$Z, col = "blue", pch = 19, main = "Mutation Frequency by 3D Spatial Coordinates")