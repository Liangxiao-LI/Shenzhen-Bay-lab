#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 16:32:10 2023

@author: ryan
"""


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

#%%define add noise function

# x is my training data
# mu is the mean
# std is the standard deviation
mu=0.0
std = 0.1
def gaussian_noise(x,mu,std):
    noise = np.random.normal(mu, std, size = x.shape)
    x_noisy = x + noise
    return x_noisy 

#%%Input data

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
vp = sub_df.iloc[:972,3]



x = np.array(x)
y = np.array(y)
z = np.array(z)

#Scale
scaler = StandardScaler()

vp = np.log10(vp)
#col = np.arange(972)

vp = np.array(vp)
vp= scaler.fit_transform(vp.reshape(-1, 1))
vp_noisy = gaussian_noise(vp,mu,std)

# Create the input feature matrix
X = np.column_stack((x, y, z))
X = scaler.fit_transform(X)

# Transform X_scaled back to the original scale
#X_original = scaler.inverse_transform(X_scaled)

#%%


# 计算 'vp' 的百分位数，找到那个较大的阈值
percentile_threshold = np.percentile(vp, 90)

# 创建样本权重，将大于阈值的 'vp' 值设置为较大的权重，其余为默认权重
sample_weights = np.where(vp > percentile_threshold, 1, 0.01)


score = -1000

for size in [0.05]:
    
    # Splitting the dataset into training and test sets
    # Random state is currently fixed
    #seed = random.randint(1, 100)
    seed = 20

    # 拆分数据集为训练集和测试集
    X_train, X_test, vp_train, vp_test, sample_weights_train, sample_weights_test = train_test_split(X, vp, sample_weights, test_size=size, random_state= seed)



#%% Random Forest Prediction

for K_num in [10]:
    
    for num in [30,50,100]:

        # Initialize the Random Forest Regressor model
        rf_model = RandomForestRegressor(random_state=seed,n_estimators = num,criterion='absolute_error')
        
        vp_train = vp_train.ravel()
        
        # Train the Random Forest model on the training data
        rf_model.fit(X_train, vp_train,sample_weight=sample_weights_train.ravel())
    
        temp_score = np.mean(cross_val_score(rf_model, X_train, vp_train, cv= K_num, scoring='r2'))
    
        if score < temp_score:
                    
            best_rf = rf_model
            score = temp_score
            K = K_num

#%% Take into test set

#vp_pred,sigma = gp.predict(X_test,return_std=True)
vp_pred = best_rf.predict(X_test)

mse = mean_squared_error(vp_test, vp_pred)
#mse_noi = mean_squared_error(vp_test_noi, vp_pred)
r2 = r2_score(vp_test, vp_pred)

print(f"The mse is {mse}, this is the mse between vp_pred and vp_test")
#print(f"The mse_noi is {mse_noi},this is the mse between vp_pred and vp_test_noi")
print(f"The R-squared is {r2}")

print(f"The seed is {seed}")
#print(f"The best_model is {gp}")
print(f"The best_model is {best_rf}")
print(f"The best_cross_val K-Fold is {K}")
print(f"The best_size is {size}")
print(f"The noise is mu={mu}, std={std}")

plt.hist(vp, bins=20, edgecolor='black')
plt.xlabel('Log 10 Virus Number')
plt.ylabel('Frequency')
plt.title('Distribution of Virus Number')
plt.show()
count = np.sum(vp_test > 0.2)
print("Number of observations in vp where vp_test > 0.5:", count)

#%% Plot vp_test,vp_pred

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


#%% Plot vp_train, vp_pred

best_rf.fit(X_train, vp_train)
vp_pred= best_rf.predict(X_train)

# Create a figure with two subplots
fig, axs = plt.subplots(1, 2, figsize=(12, 6))

# Scatter plot for vp
axs[0].scatter(range(len(vp_train)), vp_train, color='blue', marker='o')
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

#%%By this section we see that gp is perfectly matched if adding a noise
