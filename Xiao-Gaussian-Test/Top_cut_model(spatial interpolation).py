#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 11:49:35 2023

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
vp_top = vp

col = np.arange(972)

x = np.array(x)
y = np.array(y)
z = np.array(z)

vp = np.array(vp)

# 计算均值和标准差
mean_vp = np.mean(vp)
std_vp = np.std(vp)

# 计算均值加减 3 倍标准差的值
upper_limit = mean_vp + 0.5 * std_vp
lower_limit = mean_vp - 0.5 * std_vp

plt.hist(vp,10)
plt.show()

# 将大于均值加 3 倍标准差的值变为均值加 3 倍标准差
vp[vp > upper_limit] = upper_limit

# 将小于均值减 3 倍标准差的值变为均值减 3 倍标准差
vp[vp < lower_limit] = lower_limit

plt.hist(vp,10)
plt.show()

vp = np.log10(vp)

plt.hist(vp,10)
plt.show()

#Scale
scaler = StandardScaler()

vp= scaler.fit_transform(vp.reshape(-1, 1))

plt.hist(vp,10)
plt.show()

# Create the input feature matrix
X = np.column_stack((x, y, z))


#%% Cross validation find best gp

score = -1000

for size in [0.1]:
    
    # Splitting the dataset into training and test sets
    # Random state is currently fixed
    #seed = random.randint(1, 100)
    seed = 200
    
    #X_train, X_test, vp_train_noi, vp_test_noi = train_test_split(X, vp_noisy, test_size=size, random_state= seed)
    
    X_train, X_test, vp_train, vp_test = train_test_split(X, vp, test_size=size, random_state= seed)






#%%
for K_num in [10]:
    for l in [1]:
        for alpha_num in [1]:
            
            ker = RationalQuadratic(length_scale= l ,alpha=alpha_num,length_scale_bounds=(1e-100,1e100))  # You can choose other kernels as well
            
            gp = GaussianProcessRegressor(kernel=ker,
                                          optimizer='fmin_l_bfgs_b',  # Use L-BFGS-B optimizer
                                          n_restarts_optimizer=3) 
            
            #gp.fit(X_train, vp_train_noi)
            gp.fit(X_train,vp_train)
            
            print(f"The model is {gp}")
            print(f"The K-Fold is {K_num}")
            print(f"The size is {size}")
           
            
            #neg_mean_squared_error --> super smart way, since most optimization target at maximizing the score, 
            #but we want to minimize the mean_square error, which is the same as maximizing the negative mean squared error
            
            temp_score = np.mean(cross_val_score(gp, X_train, vp_train, cv= K_num, scoring='neg_mean_squared_error'))
            
            if score < temp_score:
                    
                best_gp = gp
                score = temp_score
                K = K_num



#%% Take into test set

vp_pred,sigma = gp.predict(X_test,return_std=True)
#vp_pred = best_rf.predict(X_test)


mse = mean_squared_error(vp_test, vp_pred)
#mse_noi = mean_squared_error(vp_test_noi, vp_pred)
r2 = r2_score(vp_test, vp_pred)
log= gp.log_marginal_likelihood()

print(f"The mse is {mse}, this is the mse between vp_pred and vp_test")
#print(f"The mse_noi is {mse_noi},this is the mse between vp_pred and vp_test_noi")
print(f"The R-squared is {r2}")

print(f"The seed is {seed}")
print(f"The best_model is {gp}")
#print(f"The best_model is {best_rf}")
print(f"The log_marginal_likelihood is {log}")
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


#%%

fig, ax = plt.subplots() # 创建图实例
ax.errorbar(range(len(vp_pred)), vp_pred.ravel(), yerr=sigma,ecolor='r', color='b', fmt='o', label='Uncertainty')
ax.scatter(range(len(vp_test)), vp_test, color='g', marker='o', label='vp_test')

ax.set_xlabel('Index')
ax.set_ylabel('vp_pred')
ax.set_title('Predicted vp with Uncertainty')
ax.legend()
plt.show()







#%% Top_cut model
vp_top = np.array(vp_train)

# 对vp_top数组进行从大到小排序
sorted_vp_top = np.sort(vp_top)[::-1]  # [::-1] 反转数组

# 取从大到小排序后的第10%个数
percent = 0.1

#计算ze
top_cut = sorted_vp_top[ int( np.ceil(len(vp_train)* percent))  ]  # 注意索引是从0开始的

#计算m+

# 使用布尔索引获取所有大于top_cut的值
values_above_cut = sorted_vp_top[sorted_vp_top > top_cut]


# 找到vp_train中所有对应values_above_cut的值的序号
indices = np.where(vp_train > top_cut)[0]

# 计算大于top_cut的值的平均值
m = np.mean(values_above_cut)
