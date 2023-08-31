#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 09:10:18 2023

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
#region 输入数据

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
#endregion

#%% Just a visulisation
# Create a scatter plot
plt.scatter(range(len(vp)), vp, color='blue', marker='o', label='vp')

# Add labels and title
plt.xlabel('Index')
plt.ylabel('vp')
plt.title('Scatter Plot of vp')
plt.legend()

# Show the plot
plt.show()


#%%

#vp = np.log10(vp)
vp = np.log(vp) / np.log(10)

#%% Just a visulisation
# Create a scatter plot
plt.scatter(range(len(vp)), vp, color='blue', marker='o', label='vp')

# Add labels and title
plt.xlabel('Index')
plt.ylabel('vp')
plt.title('Scatter Plot of vp')
plt.legend()

# Show the plot
plt.show()

#%%

col = np.arange(972)

x = np.array(x)
y = np.array(y)
z = np.array(z)

#Scale
scaler = StandardScaler()

vp = np.array(vp)
vp= scaler.fit_transform(vp.reshape(-1, 1))

#%% just a visulisation
# Create a scatter plot
plt.scatter(range(len(vp)), vp, color='blue', marker='o', label='vp')

# Add labels and title
plt.xlabel('Index')
plt.ylabel('vp')
plt.title('Scatter Plot of vp')
plt.legend()

# Show the plot
plt.show()

#%%

vp_noisy = gaussian_noise(vp,mu,std)

# Create the input feature matrix
X = np.column_stack((x, y, z))
X = scaler.fit_transform(X)

# Transform X_scaled back to the original scale
#X_original = scaler.inverse_transform(X_scaled)

#%% Cross validation find best gp

score = -1000

for size in [0.05]:
    
    # Splitting the dataset into training and test sets
    # Random state is currently fixed
    #seed = random.randint(1, 100)
    seed = 20
    
    X_train, X_test, vp_train_noi, vp_test_noi = train_test_split(X, vp_noisy, test_size=size, random_state= seed)
    
    X_train, X_test, vp_train, vp_test = train_test_split(X, vp, test_size=size, random_state= seed)

    

#%% Gaussian_model Best is Rational quadratic

#Rational Quadratic best: size = 0.1,seed = 20, alpha = 1, length_scale = 1
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

# Matern R-squared always negative

    # for K_num in [10]:
    #     for l in [1]:
    #         for nu_num in [2.5]:
                
    #             ker = Matern(length_scale= l , nu= nu_num,length_scale_bounds=(1e-100,1e100))  # You can choose other kernels as well
                
    #             gp = GaussianProcessRegressor(kernel=ker,
    #                                           optimizer='fmin_l_bfgs_b',  # Use L-BFGS-B optimizer
    #                                           n_restarts_optimizer=3) 
                
    #             gp.fit(X_train, vp_train_noi)
                
    #             #neg_mean_squared_error --> super smart way, since most optimization target at maximizing the score, 
    #             #but we want to minimize the mean_square error, which is the same as maximizing the negative mean squared error
    #             temp_score = np.mean(cross_val_score(gp, X_train, vp_train_noi, cv= K_num, scoring='neg_mean_squared_error'))
                
    #             if score < temp_score:
                    
    #                 best_gp = gp
    #                 score = temp_score
    #                 K = K_num

# Dot Product not good
# for K_num in [10]:
#     for sig in [5]:
        
#         ker = DotProduct(sigma_0= sig, sigma_0_bounds=(1e-05, 1e5))  # You can choose other kernels as well
        
#         gp = GaussianProcessRegressor(kernel=ker,
#                                       optimizer='fmin_l_bfgs_b',  # Use L-BFGS-B optimizer
#                                       n_restarts_optimizer=3) 
        
#         #gp.fit(X_train, vp_train_noi)
#         gp.fit(X_train,vp_train_noi)
        
#         print(f"The model is {gp}")
#         print(f"The K-Fold is {K_num}")
#         print(f"The size is {size}")
#         print(f"The sig is {sig}")
       
        
#         #neg_mean_squared_error --> super smart way, since most optimization target at maximizing the score, 
#         #but we want to minimize the mean_square error, which is the same as maximizing the negative mean squared error
        
#         temp_score = np.mean(cross_val_score(gp, X_train, vp_train_noi, cv= K_num, scoring='r2'))
        
#         if score < temp_score:
            
#             best_gp = gp
#             score = temp_score
#             K = K_num

# DotProduct kernel not good

# White Kernel not good

# Compound Kernel

#%%#By this section we see that gp is perfectly matched if adding a noise

#vp_pred,sigma = gp.predict(X_train,return_std=True)

#mse = mean_squared_error(vp_train, vp_pred)
#mse_noi = mean_squared_error(vp_train_noi, vp_pred)

#print(f"The mse is {mse}")
#print(f"The mse_noi is {mse_noi}")

#%% Random Forest Prediction

# for K_num in [10]:

#     # Initialize the Random Forest Regressor model
#     rf_model = RandomForestRegressor(random_state=seed,n_estimators = 100,criterion='absolute_error')
    
#     vp_train = vp_train.ravel()
    
#     # Train the Random Forest model on the training data
#     rf_model.fit(X_train, vp_train)

#     temp_score = np.mean(cross_val_score(rf_model, X_train, vp_train, cv= K_num, scoring='r2'))

#     if score < temp_score:
                
#         best_rf = rf_model
#         score = temp_score
#         K = K_num


#%% Take into test set

vp_pred,sigma = gp.predict(X_test,return_std=True)
#vp_pred = best_rf.predict(X_test)


mse = mean_squared_error(vp_test, vp_pred)
mse_noi = mean_squared_error(vp_test_noi, vp_pred)
r2 = r2_score(vp_test, vp_pred)
log= gp.log_marginal_likelihood()

print(f"The mse is {mse}, this is the mse between vp_pred and vp_test")
print(f"The mse_noi is {mse_noi},this is the mse between vp_pred and vp_test_noi")
print(f"The R-squared is {r2}")

print(f"The seed is {seed}")
print(f"The best_model is {gp}")
print(f"The log_marginal_likelihood is {log}")
#print(f"The best_model is {best_rf}")
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


#%% Test pred on train data for rf

# Train final model with best kernel and hyperparameters
# best_rf.fit(X_train, vp_train)

# vp_pred= best_rf.predict(X_train)
# mse = mean_squared_error(vp_train, vp_pred)
# print("mse between vp_test and vp_predict(on X_train) is" ,mse)

# fig, axs = plt.subplots(1, 2, figsize=(12, 6))
# # Scatter plot for vp
# axs[0].scatter(range(len(vp_train)), vp_train, color='blue', marker='o')
# axs[0].set_xlabel('Index')
# axs[0].set_ylabel('vp')
# axs[0].set_title('Scatter Plot of vp_train')

# # Scatter plot for vp_pred
# axs[1].scatter(range(len(vp_pred)), vp_pred, color='blue', marker='o')
# axs[1].set_xlabel('Index')
# axs[1].set_ylabel('vp_pred')
# axs[1].set_title('Scatter Plot of vp_pred')

# # Add a main title
# plt.suptitle(f"Test pred on train data", fontsize=16)

# plt.tight_layout()
# plt.show()


