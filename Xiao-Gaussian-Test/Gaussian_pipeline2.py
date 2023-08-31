#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 17:12:00 2023

@author: ryan
"""

#%%

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.colors as colors
import pandas as pd
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
from sklearn.model_selection import GridSearchCV
from sklearn.gaussian_process.kernels import RBF, Matern
import random
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error
import math
from sklearn.gaussian_process.kernels import RBF, Matern,  ExpSineSquared, DotProduct
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler
from sklearn.gaussian_process.kernels import RationalQuadratic,WhiteKernel
#from imblearn.under_sampling import RandomUnderSampler
#from imblearn.over_sampling import SMOTERegression
from sklearn.ensemble import RandomForestRegressor
import seaborn as sns


#%%
#region Parameter change
vp_threshold = 0.2 #! 用来调节必须归入test set的vp值域
k = 5  #! 后面生成test集的预测值比较分成几副小图，在'Test set检验模型好坏'的region里面
K_num_list = [10] #! 高斯模型遍历parameter K-fold cross validation 的parameter（不过现在采用的是log marginial likelihood，所以不需要设置这个，随便取一个数字就行）
length_scale_list = [1] #! 高斯模型遍历parameter （这个parameter取决于kernel，如果更换了kernel那这个也得变）
alpha_num_list = [1] #! 高斯模型遍历parameter （这个parameter取决于kernel，如果更换了kernel那这个也得变）
noise_level = 0.01 #! noise level for the gaussian process model(0.1 means standard deviation)，如果noise过高会导致residue过大，mse过大
sigma_filter = 0.15 #! 这个是最后一个region中筛选的参数，大于sigma_filter的数值都会被筛选出来
#endregion
#%% Oversampling part 1 
#region Oversampling sample code(*)
import random
from sklearn.neighbors import NearestNeighbors
import numpy as np

#! 实话实说我也不知道这块在干什么，网上找到的SMOTE样板代码，而且这一块和下面的代码也没有任何关系
#! 这一块自定义了一个Smote,看个乐呵就好哈哈哈哈哈哈哈哈哈哈哈哈
#! 下一个section用imblearn import的SMOTE才比较重要
class Smote:
    def __init__(self,samples,N=10,k=5):
        self.n_samples,self.n_attrs=samples.shape
        self.N=N
        self.k=k
        self.samples=samples
        self.newindex=0
       # self.synthetic=np.zeros((self.n_samples*N,self.n_attrs))

    def over_sampling(self):
        N=int(self.N/100)
        self.synthetic = np.zeros((self.n_samples * N, self.n_attrs))
        neighbors=NearestNeighbors(n_neighbors=self.k).fit(self.samples)
        print ('neighbors'),neighbors
        for i in range(len(self.samples)):
            nnarray=neighbors.kneighbors(self.samples[i].reshape(1,-1),return_distance=False)[0]
            #print nnarray
            self._populate(N,i,nnarray)
        return self.synthetic


    # for each minority class samples,choose N of the k nearest neighbors and generate N synthetic samples.
    def _populate(self,N,i,nnarray):
        for j in range(N):
            nn=random.randint(0,self.k-1)
            dif=self.samples[nnarray[nn]]-self.samples[i]
            gap=random.random()
            self.synthetic[self.newindex]=self.samples[i]+gap*dif
            self.newindex+=1
a=np.array([[1,2,3],[4,5,6],[2,3,1],[2,1,2],[2,3,4],[2,3,4]])
s=Smote(a,N=100)
print (s.over_sampling())

#%% Oversampling part 2
#! SMOTE样板案例应用
from imblearn.over_sampling import SMOTE
from sklearn.datasets import make_classification
import matplotlib.pyplot as plt

# 生成一个不平衡的样本数据集
X, y = make_classification(n_classes=2, class_sep=2, weights=[0.1, 0.9], n_informative=3,
                           n_redundant=1, flip_y=0, n_features=20, n_clusters_per_class=1,
                           n_samples=1000, random_state=42)

# 绘制原始数据和过采样后的数据
plt.figure(figsize=(10, 6))
plt.scatter(X[:, 0], X[:, 1], c=y, marker='o', label='Original Data')
#plt.scatter(X_resampled[:, 0], X_resampled[:, 1], c=y_resampled, marker='x', label='Resampled Data')
plt.legend()
plt.title('SMOTE Over-sampling')
plt.show()

# 创建 SMOTE 对象
smote = SMOTE(sampling_strategy='auto', random_state=42)

# 使用 SMOTE 进行过采样
X_resampled, y_resampled = smote.fit_resample(X, y)

# 绘制原始数据和过采样后的数据
plt.figure(figsize=(10, 6))
plt.scatter(X[:, 0], X[:, 1], c=y, marker='o', label='Original Data')
plt.scatter(X_resampled[:, 0], X_resampled[:, 1], c=y_resampled, marker='x', label='Resampled Data')
plt.legend()
plt.title('SMOTE Over-sampling')
plt.show()
#endregion
#%% Load data
#region Load data
#Parameter change
#纯手工K-Fold section K：

df = pd.read_excel('/Users/ryan/Documents/GitHub/Shenzhen Bay lab/spike_protein_3structures_markerd_with_variants_from_cncb/6vxx_variants.xls')
df_CYS = pd.read_excel('/Users/ryan/Documents/GitHub/Shenzhen Bay lab/spike_protein_3structures_markerd_with_variants_from_cncb/6vxx_variants_CYS.xls')

sub_df = df.iloc[:972, 5:10]
sub_df_CYS = df.iloc[:972, 5:10]

# Extract the x, y, z, and m columns from the subset DataFrame
x = sub_df.iloc[:972,0]
y = sub_df.iloc[:972,1]
z = sub_df.iloc[:972,2]
vp = sub_df.iloc[:972,4]

#endregion

#region SMOTE过采样
#!（oversampling是针对classification problem时候遇到imbalance set设计的一种方法）
#Create the label array based on the conditions
#labels = np.where(vp > 0.5, 2, np.where((vp > 0.1) & (vp <= 0.5), 1, 0))

labels = np.where(vp > 0.5, 3, 
                  np.where((vp > 0.3) & (vp <= 0.5), 2,
                  np.where((vp > 0.1) & (vp <= 0.3), 1, 0)))


X = np.column_stack((x, y, z))

# Initialize SMOTE with the desired sampling strategy
smote = SMOTE(sampling_strategy='auto', random_state=42)

# Perform SMOTE oversampling
X_resampled, labels_resampled = smote.fit_resample(X, labels)

#we discover that the added oversampling values are all behind the original values
matching_indices = np.where((X[:, None] == X_resampled).all(-1))[1]

vp_resampled = []
for i in X_resampled:
    
    if i in X:
        
        indices = (np.where(np.all(X == i, axis=1))[0])[0]
        
        vp_resampled.append(vp[indices])
    
    else:
        
        indices = (np.where(np.all(X_resampled == i, axis=1))[0])[0]
        
        # if labels_resampled[indices] == 2:
            
        #     vp_resampled.append(np.random.uniform(0.5, 1))
        
        # if labels_resampled[indices] == 1:
            
        #     vp_resampled.append(np.random.uniform(0.1, 0.5))
        
        # if labels_resampled[indices] == 0:
            
        #     vp_resampled.append(np.random.uniform(0, 0.1))
            
        if labels_resampled[indices] == 3:
            
            vp_resampled.append(np.random.uniform(0.5, 1))
        
        if labels_resampled[indices] == 2:
            
            vp_resampled.append(np.random.uniform(0.3, 0.5))
        
        if labels_resampled[indices] == 1:
            
            vp_resampled.append(np.random.uniform(0.1, 0.3))
        
        if labels_resampled[indices] == 0:
            
            vp_resampled.append(np.random.uniform(0, 0.1))
            
     
vp_resampled = np.array(vp_resampled)
        
#endregion

#%% 划分train set与test set
#region 划分train set与test set
#! 注意这个部分里有一个额外操作，就是把大于vp_threshold的train set中的数值给转移到了test set中，如果不需要的话可以comment掉

score = -1000

for size in [0.05]:
    
    # Splitting the dataset into training and test sets
    # Random state is currently fixed
    #seed = random.randint(1, 100)
    seed = 20
    
    X_train, X_test, vp_train, vp_test = train_test_split(X_resampled, vp_resampled, test_size=size, random_state= seed)

#we discover that the added oversampling values are all behind the original values
vp_selected = vp[vp > vp_threshold]
X_selected = X[vp_selected.index]
matching_indices = np.where((X_selected[:, None] == X_train).all(-1))[1]

# 将这些行添加到X_test和vp_test
X_test = np.concatenate((X_test, X_train[matching_indices]), axis=0)
vp_test = np.concatenate((vp_test, vp_train[matching_indices]), axis=0)

# 从X_train和vp_train中删除matching_indices指定的行
X_train = np.delete(X_train, matching_indices, axis=0)
vp_train = np.delete(vp_train, matching_indices, axis=0)

#endregion
#%%  GP model train
#region GP model train

#! 这里是用这个循环而不是使用GridSearchCV或者RandomizedSearchCV的原因是：Grid或者Randomized中不包含log marginal likelihood作为scoring function
#! 这边我用循环写了一下这个hyperparameter tuning，以logmarginal likelihood作为K-fold cross validation的score
for K_num in K_num_list:
    for l in length_scale_list:
        for alpha_num in alpha_num_list:
            
            #! Kernel是非常重要的一个参数，这个也是hyperparameter调整的重要依据
            #! 其他可选的kernel比如有
            ker = RationalQuadratic(length_scale= l ,alpha=alpha_num,length_scale_bounds=(1e-100,1e100))  # You can choose other kernels as well
            
            #! alpha是GP model的噪点
            gp = GaussianProcessRegressor(kernel=ker,
                                            alpha= noise_level,
                                            optimizer='fmin_l_bfgs_b',  # Use L-BFGS-B optimizer
                                            n_restarts_optimizer=3) 
            
            kf = KFold(n_splits=K_num)
            temp_scores = []
            
            for train_index, test_index in kf.split(X_train):
                X_train_kf, X_test_kf = X_train[train_index], X_train[test_index]
                vp_train_kf, vp_test_kf = vp_train[train_index], vp_train[test_index]
                
                gp.fit(X_train_kf, vp_train_kf)
                
                
                temp_scores.append(gp.log_marginal_likelihood())
            
            temp_score = np.mean(temp_scores)
            
            if score < temp_score:
                    
                best_gp = gp
                score = temp_score
                K = K_num

#endregion

#%% Take into test set
#region Test set 检验模型好坏

vp_pred,sigma = best_gp.predict(X_test,return_std=True)
#vp_pred = best_rf.predict(X_test)


mse = mean_squared_error(vp_test, vp_pred)
#mse_noi = mean_squared_error(vp_test_noi, vp_pred)
r2 = r2_score(vp_test, vp_pred)

print(f"The mse is {mse}, this is the mse between vp_pred and vp_test")
#print(f"The mse_noi is {mse_noi},this is the mse between vp_pred and vp_test_noi")
print(f"The R-squared is {r2}")
print(f"The noise_level is {noise_level}")
print(f"The seed is {seed}")
print(f"The best_model is {best_gp}")
print(f"The log_marginal_likelihood is {score}")
#print(f"The best_model is {best_rf}")
print(f"The best_cross_val K-Fold is {K}")
print(f"The best_size is {size}")
#print(f"The noise is mu={mu}, std={std}")

plt.hist(vp, bins=20, edgecolor='black')
plt.xlabel('Log 10 Virus Number')
plt.ylabel('Frequency')
plt.title('Distribution of Virus Number')
plt.show()
count = np.sum(vp_test > 0.2)
print("Number of observations in vp where vp_test > 0.5:", count)


#%% Plot

# Create a figure with two subplots
fig, axs = plt.subplots(1, 2, figsize=(12, 6))

# Scatter plot for vp
axs[0].scatter(range(len(vp_test)), vp_test, color='blue', marker='o')
axs[0].set_xlabel('Index')
axs[0].set_ylabel('vp')
axs[0].set_title('Scatter Plot of vp_test')

# Scatter plot for vp_pred

axs[1].scatter(range(len(vp_pred)), vp_pred, color='blue', marker='o')
axs[1].set_xlabel('Index')
axs[1].set_ylabel('vp_pred')
axs[1].set_title('Scatter Plot of vp_pred')

# Add a main title
plt.suptitle(f"Test pred on test data", fontsize=16)

plt.tight_layout()
plt.show()



#%% 生成predict和test set之间的差距比较

fig, ax = plt.subplots() # 创建图实例
ax.errorbar(range(len(vp_pred)), vp_pred.ravel(), yerr=sigma,ecolor='r', color='b', fmt='o', label='Uncertainty')
ax.scatter(range(len(vp_test)), vp_test, color='g', marker='o', label='vp_test')

ax.set_xlabel('Index')
ax.set_ylabel('vp_pred')
ax.set_title('Predicted vp with Uncertainty')
ax.legend()
plt.show()

#%% 把上面这幅图分成k等份（几份小图）

#k = 5 # 分成k幅图
n = len(vp_pred) // k

matching_indices_test = np.where((X_selected[:, None] == X_test).all(-1))[1]

fig, axs = plt.subplots(k, 1, figsize=(10, k*5))  # 创建k个子图

for i in range(k):
    start = i * n
    end = start + n if i != k - 1 else len(vp_pred)
    
    for j in range(start, end):
        marker = '*' if j in matching_indices_test else 'o'
        color = 'c' if j in matching_indices_test else 'g'
        axs[i].scatter(j, vp_test[j], color=color, marker=marker, label='vp_test')
    
    axs[i].errorbar(range(start, end), vp_pred[start:end].ravel(), yerr=sigma[start:end], ecolor='r', color='b', fmt='o', label='Uncertainty')
    
    axs[i].set_xlabel('Index')
    axs[i].set_ylabel('vp_pred')
    axs[i].set_title(f'Predicted vp with Uncertainty (Part {i+1})')
    #axs[i].legend()

plt.tight_layout()
plt.show()
#%% visulise sigma

plt.scatter(sigma, vp_pred, alpha=0.5)
plt.xlabel('sigma')
plt.ylabel('vp_pred')
plt.title('Scatter plot of sigma vs vp_pred')
plt.show()

plt.hist(sigma, bins=20, edgecolor='black')
plt.xlabel('sigma')
plt.ylabel('Frequency')
plt.title('Distribution of sigma')
plt.show()

plt.boxplot(sigma)
plt.ylabel('sigma')
plt.title('Boxplot of sigma')
plt.show()

# Create a new figure
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Scatter plot with sigma controlling the color
scatter = ax.scatter(X_test[:, 0], X_test[:, 1], X_test[:, 2], c=sigma, cmap='coolwarm')

# Add a colorbar
cbar = plt.colorbar(scatter)
cbar.set_label('sigma')

# Set labels and title
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.title('3D Scatter Plot with sigma')

# Show the plot
plt.show()

#%% show the correlation between each covariates

# Calculate the correlation matrix
correlation_matrix = np.corrcoef(sigma, X_test.T)

# Create a heatmap
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt=".2f")

# Set the axis labels and title
plt.xlabel('Features of X_test')
plt.ylabel('sigma')
plt.title('Correlation Heatmap between sigma and X_test')

# Show the plot
plt.show()

#endregion
#%% Take into original non-resampled set
#region 将整个Train set放入模型中，查看uncertainty与residue
vp_pred,sigma = gp.predict(X,return_std=True)
#vp_pred = best_rf.predict(X_test)

mse = mean_squared_error(vp, vp_pred)
#mse_noi = mean_squared_error(vp_test_noi, vp_pred)
r2 = r2_score(vp, vp_pred)

print(f"The noise_level is {noise_level}")
print(f"The mse is {mse}, this is the mse between vp_pred and vp")
#print(f"The mse_noi is {mse_noi},this is the mse between vp_pred and vp_test_noi")
print(f"The R-squared is {r2}")

print(f"The seed is {seed}")
print(f"The best_model is {gp}")
#print(f"The best_model is {best_rf}")
print(f"The best_cross_val K-Fold is {K}")
print(f"The best_size is {size}")
#print(f"The noise is mu={mu}, std={std}")

plt.hist(vp, bins=20, edgecolor='black')
plt.xlabel('Log 10 Virus Number')
plt.ylabel('Frequency')
plt.title('Distribution of Virus Number')
plt.show()
count = np.sum(vp_test > 0.2)
print("Number of observations in vp where vp_test > 0.5:", count)
#endregion
#%% Plot
#region sigma探索分析
# Create a figure with two subplots
fig, axs = plt.subplots(1, 2, figsize=(12, 6))

# Scatter plot for vp
axs[0].scatter(range(len(vp)), vp, color='blue', marker='o')
axs[0].set_xlabel('Index')
axs[0].set_ylabel('vp')
axs[0].set_title('Scatter Plot of vp_test')

# Scatter plot for vp_pred

axs[1].scatter(range(len(vp_pred)), vp_pred, color='blue', marker='o')
axs[1].set_xlabel('Index')
axs[1].set_ylabel('vp_pred')
axs[1].set_title('Scatter Plot of vp_pred')

# Add a main title
plt.suptitle(f"Test pred on test data", fontsize=16)

plt.tight_layout()
plt.show()

#%%

fig, ax = plt.subplots() # 创建图实例
ax.errorbar(range(len(vp_pred)), vp_pred.ravel(), yerr=sigma,ecolor='r', color='b', fmt='o', label='Uncertainty')

ax.scatter(range(len(vp)), vp, color='g', marker='o', label='vp')

ax.set_xlabel('Index')
ax.set_ylabel('vp_pred')
ax.set_title('Predicted vp with Uncertainty')
ax.legend()
plt.show()

#%% visulise sigma

plt.scatter(sigma, vp, alpha=0.5)
plt.xlabel('sigma')
plt.ylabel('vp')
plt.title('Scatter plot of sigma vs vp')
plt.show()

plt.hist(sigma, bins=20, edgecolor='black')
plt.xlabel('sigma')
plt.ylabel('Frequency')
plt.title('Distribution of sigma')
plt.show()

plt.boxplot(sigma)
plt.ylabel('sigma')
plt.title('Boxplot of sigma')
plt.show()

# Create a new figure
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Scatter plot with sigma controlling the color
scatter = ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=sigma, cmap='coolwarm')

# Add a colorbar
cbar = plt.colorbar(scatter)
cbar.set_label('sigma')

# Set labels and title
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.title('3D Scatter Plot with sigma')

# Show the plot
plt.show()

#%% show the correlation between each covariates

# Calculate the correlation matrix
correlation_matrix = np.corrcoef(sigma, X.T)

# Create a heatmap
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt=".2f")

# Set the axis labels and title
plt.xlabel('Features of X')
plt.ylabel('sigma')
plt.title('Correlation Heatmap between sigma and X_test')

# Show the plot
plt.show()

#%% show the correlation between each covariates

# Calculate the correlation matrix
correlation_matrix = np.corrcoef(sigma, vp)

# Create a heatmap
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt=".2f")

# Set the axis labels and title
plt.xlabel('Features of X')
plt.ylabel('sigma')
plt.title('Correlation Heatmap between sigma and vp')

# Show the plot
plt.show()

#endregion

# %% 计算出residue
#region residue探索分析
residue = vp_pred - vp

# 创建散点图
plt.figure(figsize=(10, 6))
plt.scatter(residue, vp, alpha=0.5)
plt.xlabel('residue')
plt.ylabel('vp')
plt.title('Scatter plot of residue vs vp')
plt.show()

# 创建sigma的直方图
plt.figure(figsize=(10, 6))
plt.hist(residue, bins=30, edgecolor='black')
plt.xlabel('residue')
plt.ylabel('Frequency')
plt.title('Histogram of residue')
plt.show()

plt.boxplot(residue)
plt.ylabel('residue')
plt.title('Boxplot of residue')
plt.show()
#%% show the correlation between each covariates

# Calculate the correlation matrix
correlation_matrix = np.corrcoef(residue, vp)

# Create a heatmap
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt=".2f")

# Set the axis labels and title
plt.xlabel('vp')
plt.ylabel('residue')
plt.title('Correlation Heatmap between residue and vp')

# Show the plot
plt.show()

#%% show the correlation between each covariates

# Calculate the correlation matrix
correlation_matrix = np.corrcoef(residue, X.T)

# Create a heatmap
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt=".2f")

# Set the axis labels and title
plt.xlabel('Features of X')
plt.ylabel('residue')
plt.title('Correlation Heatmap between residue and Features of X')

# Show the plot
plt.show()
# %%
# Create a new figure
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Scatter plot with sigma controlling the color
scatter = ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=residue, cmap='coolwarm')

# Add a colorbar
cbar = plt.colorbar(scatter)
cbar.set_label('residue')

# Set labels and title
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.title('3D Scatter Plot with residue')

# Show the plot
plt.show()

#endregion

# %% save output data
#region 筛选出sigma大于0.15的数据

#!数据整合， 把residue，sigma和感兴趣的数据整合在一块，输出成predict_on_X
temp_vp = sub_df.iloc[:972,4]

data = {
    'X':x,
    'Y':y,
    'Z':z,
    'vp':temp_vp,
    'residue': residue,
    'predict_vp': vp_pred,
    'sigma': sigma
}

df_combined = pd.DataFrame(data)
df_combined.to_csv('predict_on_X.csv', index=False)
# %% Greater_sigma 大于0.15

#!筛选出 value > sigma_filter

sigma_list = data['sigma'].tolist()

#! sigma filter 在这里
filtered_data = [value for value in sigma_list if value > sigma_filter]

indices = [index for index, value in enumerate(sigma_list) if value in filtered_data]

selected_vp = [vp[index] for index in indices]

#! 探索性分析: 看看这些筛选出来的sigma较大的数据是不是都来自test set或者是train set，结果发现都有

intersection_test = set(selected_vp) & set(vp_test)
intersection_train = set(selected_vp) & set(vp_train)
count_test = len(intersection_test)
count_train = len(intersection_train)

print(f"The number of elements in 'sigma > 0.15' also present in 'test set' is: {count_test}")
print(f"The number of elements in 'sigma > 0.15' also present in 'train set' is: {count_train}")
print(f"The size of the 'sigma > 0.15' is {len(filtered_data)}")

#%%

sigma_list = data['sigma'].tolist()

filtered_data = [value for value in sigma_list if  0.07<value <0.11]

indices = [index for index, value in enumerate(sigma_list) if value in filtered_data]

selected_vp = [vp[index] for index in indices]

intersection_test = set(vp) & set(vp_test)
intersection_train = set(vp) & set(vp_train)
count_test = len(intersection_test)
count_train = len(intersection_train)

print(f"The number of elements in 'vp_true' also present in 'test set' is: {count_test}")
print(f"The number of elements in 'vp_true' also present in 'train set' is: {count_train}")
print(f"The size of the 'vp_true' is {len(vp)}")
#%%

df_greater_sigma = pd.read_csv('/Users/ryan/Documents/GitHub/Shenzhen Bay lab/Xiao-Gaussian-Test/greater_sigma.csv')

x_g = df_greater_sigma.iloc[:,0]
y_g = df_greater_sigma.iloc[:972,1]
z_g = df_greater_sigma.iloc[:972,2]
sigma_g = df_greater_sigma.iloc[:972,6]

X_g = np.column_stack((x_g, y_g, z_g))

# Create a new figure
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Scatter plot with sigma controlling the color
scatter = ax.scatter(X_g[:, 0], X_g[:, 1], X_g[:, 2], c=sigma_g, cmap='coolwarm')

# Add a colorbar
cbar = plt.colorbar(scatter)
cbar.set_label('sigma')

# Set labels and title
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.title('3D Scatter Plot with greater_sigma')

# Show the plot
plt.show()
# %%

df_greater_sigma = df_greater_sigma.dropna()


# %%
df = pd.read_excel('/Users/ryan/Documents/GitHub/Shenzhen Bay lab/spike_protein_3structures_markerd_with_variants_from_cncb/6vxx_variants.xls')

A_chain = df.iloc[:972,]

merged_df = df_greater_sigma.merge(A_chain, how='left')

merged_df.to_csv('greater_sigma_full_data.csv', index=False)


# %% Greater_sigma 小于0

sigma_list = data['sigma'].tolist()

filtered_data = [value for value in sigma_list if value > 0.15]

indices = [index for index, value in enumerate(sigma_list) if value in filtered_data]

selected_vp = [vp[index] for index in indices]

intersection_test = set(selected_vp) & set(vp_test)
intersection_train = set(selected_vp) & set(vp_train)
count_test = len(intersection_test)
count_train = len(intersection_train)

print(f"The number of elements in 'sigma > 0.15' also present in 'test set' is: {count_test}")
print(f"The number of elements in 'sigma > 0.15' also present in 'train set' is: {count_train}")
print(f"The size of the 'sigma > 0.15' is {len(filtered_data)}")

#endregion