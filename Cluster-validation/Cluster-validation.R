library(factoextra)
library(fpc)
library(NbClust)
library(readxl)

# Excluding the column "Species" at position 5
df_0 <- iris[, -5]
# Standardize
df <- scale(df_0)

# Retrieve the scaling factors and column means
scale_factor <- attr(df, "scaled:scale")
col_means <- attr(df, "scaled:center")

# Convert df_0 back to the original data frame
df_new <- as.data.frame(t( t(df)* scale_factor + col_means))

# K-means clustering
km.res <- eclust(df, "kmeans", k = 3, nstart = 25, graph = FALSE)
# Visualize k-means clusters(2d)
fviz_cluster(km.res, geom = "point", ellipse.type = "norm",
             palette = "jco", ggtheme = theme_minimal())

# Hierarchical clustering
hc.res <- eclust(df, "hclust", k = 3, hc_metric = "euclidean", 
                 hc_method = "ward.D2", graph = FALSE)

# Visualize dendrograms
fviz_dend(hc.res, show_labels = FALSE,
          palette = "jco", as.ggplot = TRUE)

#
fviz_silhouette(km.res)


# Silhouette information
silinfo <- km.res$silinfo
names(silinfo)
# Silhouette widths of each observation
head(silinfo$widths[, 1:3], 10)
# Average silhouette width of each cluster
silinfo$clus.avg.widths
# The total average (mean of all individual silhouette widths)
silinfo$avg.width
# The size of each clusters
km.res$size

# Silhouette width of observation
sil <- km.res$silinfo$widths[, 1:3]
# Objects with negative silhouette
neg_sil_index <- which(sil[, 'sil_width'] < 0)
sil[neg_sil_index, , drop = FALSE]

#Validation statistics
library(fpc)
# Statistics for k-means clustering
km_stats <- cluster.stats(dist(df),  km.res$cluster)
# Dun index
km_stats$dunn
km_stats



########
library(scatterplot3d)

setwd("/Users/ryan/Documents/GitHub/Shenzhen Bay lab/spike_protein_3structures_markerd_with_variants_from_cncb")

data <- read_excel("6vxx_variants.xls")
df <- data[1:972,6:8]

df <- scale(df)

# Retrieve the scaling factors and column means
scale_factor <- attr(df, "scaled:scale")
col_means <- attr(df, "scaled:center")

# Method1: K-means clustering
km.res <- eclust(df, "kmeans", k = 100, nstart = 25, graph = FALSE)

# Visualize k-means clusters(3d)

scatterplot3d(df[, 1], df[, 2], df[, 3], color = km.res$cluster)

# Method2: Hierarchical clustering

hc.res <- eclust(df, "hclust", k = 3, hc_metric = "euclidean", 
                 hc_method = "ward.D2", graph = FALSE)

scatterplot3d(df[, 1], df[, 2], df[, 3], color = hc.res$cluster)

#It can be seen that several samples, in cluster 2, have a negative silhouette coefficient. This means that they are not in the right cluster.
#Si 为负的点说明聚类错误，在后面有提及
#横坐标是每一个数据点，纵坐标是Si，总平均S以及每一个Si越接近1越好

fviz_silhouette(km.res)
fviz_silhouette(hc.res)

#km

# Silhouette information
silinfo <- km.res$silinfo
names(silinfo)
# Silhouette widths of each observation
head(silinfo$widths[, 1:3], 100)
# Average silhouette width of each cluster
silinfo$clus.avg.widths
# The total average (mean of all individual silhouette widths)
silinfo$avg.width
# The size of each clusters
km.res$size

#Find negative observations
# Silhouette width of observation
sil <- km.res$silinfo$widths[, 1:3]
# Objects with negative silhouette
neg_sil_index <- which(sil[, 'sil_width'] < 0)
sil[neg_sil_index, , drop = FALSE]

#Validation statistics ( Super important )

# Statistics for k-means clustering
km_stats <- cluster.stats(dist(df),  km.res$cluster)
# Dunn index
dunn <- km_stats$dunn
# Calinski-Harabasz index
ch <- km_stats$ch
# Average silhouette width
sil <- km_stats$avg.silwidth
#Entropy index
ent <- km_stats$entropy

# number of noise points :Noise points are data points that do not belong to any specific cluster or are considered outliers
km_stats$noisen
# vector of clusterwise within cluster average distances. (这个是scaled version的average distance， non-scaled version到时候还得再算一下)
km_stats$average.distance
# Calculate the separation matrix. provides a summary of the between-cluster distances and can help identify clusters that are well-separated or overlapping.
#The separation matrix is a symmetric matrix that quantifies the degree of separation between clusters based on some distance measure. 
separation_matrix <- km_stats$separation.matrix

index <- data.frame()
index <- rbind(index, c(dunn, ch, sil, ent))
index <- rbind(index, c(dunn*2, ch*3, sil*4, ent*5))
index
colnames(index) <- c("Dunn Index", "Calinski-Harabasz Index", "Average Silhouette Width", "Entropy")
index

par(mfrow = c(2, 2))

plot(index$`Dunn Index`, type = "l", xlab = "Observation", ylab = "Dunn Index", main = "Dunn Index")

# Line chart for Calinski-Harabasz Index
plot(index$`Calinski-Harabasz Index`, type = "l", xlab = "Observation", ylab = "Calinski-Harabasz Index", main = "Calinski-Harabasz Index")

# Line chart for Average Silhouette Width
plot(index$`Average Silhouette Width`, type = "l", xlab = "Observation", ylab = "Average Silhouette Width", main = "Average Silhouette Width")

# Line chart for Entropy
plot(index$Entropy, type = "l", xlab = "Observation", ylab = "Entropy", main = "Entropy")
