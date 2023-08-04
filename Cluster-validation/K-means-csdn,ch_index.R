#K-means function 创建
K_means <- function(data,k,max.iter = 10){
  rows <- nrow(data)                 #获取行数
  cols <- ncol(data)                 #获取列数
  within <- numeric(k)               #用于存储组类平方和
  lable_matrix <- matrix(0,rows,2)   #用于存储类标签及到类中心的距离
  centers <- matrix(0,cols,k)        #用于存储类中心
  centers_matrix <- matrix(0,rows,k) #用于存储初始确定初始类中心时到类中心的距离
  iter <- 0                          #迭代次数
  random <- sample(1:rows,1)
  centers[,1] <- as.matrix(data[random,])
  for(j in 2:k){
    for(i in 1:rows){
      centers_matrix[i,j] <- sum((data[i,] - centers[,j-1])^2)+centers_matrix[i,j-1]
    }
    centers[,j] <- as.matrix(data[which(centers_matrix[,j] == max(centers_matrix[,j])),])
  }                                  #计算初始类中心
  changed <- TRUE                    #用于判断数据的类标签是否发生改变
  while(changed){
    if(iter >= max.iter){
      changed <- FALSE
      break
    }
    for(j in 1:rows){
      updata <- 1000000000
      for( i in 1:k){
        updatas <- sum((data[j,]-centers[,i])^2)
        if(updatas < updata){
          updata <- updatas
          lable_matrix[j,1] <- updatas
          lable_matrix[j,2] <- i
        }
      }
    }                                 #更新到类中心的距离以及类标签
    center <- centers
    for(i in 1:k){
      centers[,i] <- colMeans(data[lable_matrix[,2]==i,])
    }                                 #更新类中心
    changed <- !all(center == centers)#判断类中心是否变化
    iter <- iter + 1
  }
  ###计算函数返回值：
  totss <- sum((t(data[,])-colMeans(data))^2)
  withinss <- numeric()
  for(i in 1:k){
    withinss[i] <- sum((t(data[lable_matrix[,2]==i,])-(centers[,i]))^2)
  }
  tot.withinss <- sum(withinss)
  betweenss<-0
  for(i in 1:3){
    betweenss <- betweenss + sum(nrow(data[lable_matrix[,2]==i,])*(rowMeans(t(data[lable_matrix[,2]==i,]))-colMeans(data))^2)
  }
  size <- aggregate(lable_matrix[,2], by=list(lable=lable_matrix[,2]),length)[,2]
  
  centers_matrix <- t(centers)
  colnames(centers_matrix) <- colnames(data)
  result <- list(cluster = lable_matrix[,2],centers = centers_matrix,totss = totss,withinss = withinss,tot.withinss = tot.withinss
                 ,betweenss = betweenss,size = size,iter = iter)
  return(result)
}

df <- kmeans(iris[,1:4],3)
#参数：data：用于聚类的数据；k：用于聚类的数目；max.iter：最大迭代次数

#输出结果（均与R语言中自带kmeans函数输出结果命名一致）：
df$cluster
#cluster：聚类结果，即类标签；

df$centers
#centers：聚类中心；

df$totss
#totss：总平方和；

df$withinss
#withinss：各组内的平方和；

df$tot.withinss
#tot.withinss:组内平方和；

df$betweenss
#betweenss：组间平方和；

df$size
#size：各类的数量；

df$iter
#iter ：迭代次数；

K_means(iris[,1:4],3)

#可视化模版
library(ggplot2)
library(rgl)

mycolors <- c('royalblue1', 'darkcyan', 'oldlace')
iris$color <- mycolors[ as.numeric(iris$Species) ]
#绘制
plot3d( 
  x=iris$`Sepal.Length`, y=iris$`Sepal.Width`, z=iris$`Petal.Length`, 
  col = iris$color, 
  type = 's', 
  radius = .1,
  xlab="Sepal Length", ylab="Sepal Width", zlab="Petal Length")

#旋转图片至合适角度截图
rgl.postscript("plot1.pdf", "pdf", drawText = TRUE)
rgl.snapshot( "snapshot.png", fmt = "png")




#获取数据
library(readxl)
library(cluster)
library(fpc)


setwd("/Users/ryan/Documents/GitHub/Shenzhen Bay lab/spike_protein_3structures_markerd_with_variants_from_cncb")

data <- read_excel("6vxx_variants.xls")
data <- data[1:972,]

df <- kmeans(data[,6:8],50)
df$cluster
df$centers
df$totss
df$withinss
df$tot.withinss
df$betweenss
df$size
df$iter

ch_index <- cluster.stats(data[,6:8], as.matrix(df$cluster))$ch


