#空间相关性分析
library(spatstat)
library(ggplot2)
library(leaps)

# 读取数据
data <- read.csv("data1.csv")
data2<-data[1:972,c("X","Y","Z","virusPercent")]

# 分割数据集
library(caTools)
set.seed(123) # 设置随机数种子，以便结果可重复
#在机器学习中，我们通常需要将数据分为训练集和测试集，
#以便训练模型和评估模型的性能。
#如果我们将整个数据集分割为训练集和测试集，
#那么训练集和测试集中的数据将可能来自不同的数据分布，
#导致模型在测试集上的表现不够准确。
#因此，我们需要根据某一列数据（比如标签数据）进行分割，
#以保证训练集和测试集中的数据分布相同。
#具体来说，我们通常会按照某一列数据的取值来划分数据集，
#比如将数据集按照标签值进行划分。这样，训练集和测试集中的数据分布就相同了，从而可以更准确地评估模型的性能。
#在回归问题中，我们通常会将数据集按照目标变量y的取值进行分割，
#以保证训练集和测试集中的目标变量分布相同。
split <- sample.split(data2$virusPercent, SplitRatio = 0.7) # 将数据分为70%的训练集和30%的测试集
train_data <- subset(data2, split == TRUE) # 选择训练集
test_data <- subset(data2, split == FALSE) # 选择测试集

#可视化单独每个坐标和突变频率之间的散点图
plot(train_data)

a<-cor(train_data$virusPercent,train_data$X)
b<-cor(train_data$virusPercent,train_data$Y)
c<-cor(train_data$virusPercent,train_data$Z)
#做一个可视化表格

#find best subsets and original mdoel
reg_sub<-regsubsets(virusPercent~.,data = train_data)
summary.out<-summary(reg_sub)
summary.out
summary.out$cp
plot(reg_sub,scale = 'Cp')

plot(reg_sub,scale = 'bic')

fit1<-lm(virusPercent~.,data = train_data)
fit_step1<-step(fit1)
fit2<-lm(virusPercent~1,data = train_data)
fit_step2<-step(fit2,scope = virusPercent~X+Y+Z)
predictions<-predict(fit_step1,newdata = select(test_data,-virusPercent))
mse_step1<-mean((predictions-BrozekResponses)^2)
predictions<-predict(fit_step2,newdata = select(TrainData,-brozek))
mse_step2<-mean((predictions-BrozekResponses)^2)
mse_step1
mse_step2