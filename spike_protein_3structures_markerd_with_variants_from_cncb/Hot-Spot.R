library(readxl)

setwd("/Users/ryan/Documents/GitHub/Shenzhen Bay lab/spike_protein_3structures_markerd_with_variants_from_cncb")

data <- read_excel("6vxx_variants.xls")

library(cluster)

data_mutation<-data[1:972,10]

k <- 3  # 聚类数目
kmeans_result <- kmeans(data_mutation, centers = k)

# 获取聚类结果
cluster_labels <- kmeans_result$cluster

#data_mutation$cluster <- as.factor(cluster_labels)

# 可视化聚类结果

clusplot(data_mutation, cluster_labels, color = TRUE, labels = 1,shape = TRUE)

