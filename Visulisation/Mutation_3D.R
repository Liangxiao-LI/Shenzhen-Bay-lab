library(readxl)
library(dplyr)
library(ggplot2)
library(rgl)
library(grDevices)

setwd("/Users/ryan/Documents/GitHub/Shenzhen Bay lab/spike_protein_3structures_markerd_with_variants_from_cncb")

data <- read_excel("6vxx_variants.xls")

x <- data[1:972,6]
y <- data[1:972,7]
z <- data[1:972,8]
values <- data[1:972,9]

x <- unlist(x)
x <- as.numeric(x)
y <- unlist(y)
y <- as.numeric(y)
z <- unlist(z)
z <- as.numeric(z)
values <- unlist(values)
values <- as.numeric(values)
values <- log(values)

colors <- colorRampPalette(c("lightblue", "darkblue"))(max(values))
col_vec <- numeric()
for (i in 1:972){ 
  
  col_vec <- c(col_vec, colors[values[i]])
  
  }

# 绘制散点图
plot3d(x, y, z, col = colors, size = 2)
