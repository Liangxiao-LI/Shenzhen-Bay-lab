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
from sklearn.decomposition import PCA
#%%
noise_level = 0.01 #! noise level for gaussian process (Since our dataset is small, it would be better to choose a rather small noise level)
K_num_list = [10] #! list of K-fold that you want to iterate with
length_scale_list = [1] #! 高斯模型Rational quadatic kernel遍历parameter （取决于kernel）
alpha_num_list = [1] #! 高斯模型Rational quadatic kernel遍历parameter （取决于kernel）
cv_size = [0.05] #! cross validation test_size
capping_percent = 90 #! This means the top 10% value will be capped to the same value as the 90th percentage

#%%Input data
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
vn = sub_df.iloc[:972,3]
#endregion

#%% Just a visulisation
# Create a scatter plot
#region Data preprocessing

plt.scatter(range(len(vn)), vn, color='blue', marker='o', label='vn')

# Add labels and title
plt.xlabel('Index')
plt.ylabel('vn')
plt.title('Scatter Plot of vn')
plt.legend()

# Show the plot
plt.show()


#%%
#region log transform
#vn = np.log10(vn)
vn = np.log10(vn)

#%% Just a visulisation
# Create a scatter plot
plt.scatter(range(len(vn)), vn, color='blue', marker='o', label='vn')

# Add labels and title
plt.xlabel('Index')
plt.ylabel('vn')
plt.title('Scatter Plot of log10(vn)')
plt.legend()

# Show the plot
plt.show()
#endregion
#%% z-score normalization
#region z-score
col = np.arange(972)

x = np.array(x)
y = np.array(y)
z = np.array(z)

#Scale
scaler = StandardScaler()

vn = np.array(vn)
vn= scaler.fit_transform(vn.reshape(-1, 1))
#%% just a visulisation for the standarization
# Create a scatter plot
plt.scatter(range(len(vn)), vn, color='blue', marker='o', label='vn')

# Add labels and title
plt.xlabel('Index')
plt.ylabel('vn')
plt.title('Scatter Plot of standardized log10(vn)')
plt.legend()

# Show the plot
plt.show()

#%%
plt.hist(vn, bins=10, color='blue', edgecolor='black')

# Add labels and title
plt.xlabel('vn')
plt.ylabel('Frequency')
plt.title('Distribution of vn')

# Show the plot
plt.show()

#endregion
#%% capping
#region capping
# Calculate the upper 90% value
upper = np.percentile(vn, capping_percent)

# Set the top 10% values to be equal to the upper 90% value with a random fluctuation(in order to prevent later cancelation)
for i in range(len(vn)):
    if vn[i][0] > upper:
        vn[i] = upper + 0.0001*np.random.random()

#%% just a visulisation for the capping
# Create a scatter plot
plt.scatter(range(len(vn)), vn, color='blue', marker='o', label='vn')

# Add labels and title
plt.xlabel('Index')
plt.ylabel('vn')
plt.title('Scatter Plot of capped standardized log10(vn)')
plt.legend()

# Show the plot
plt.show()
#%%
plt.hist(vn, bins=10, color='blue', edgecolor='black')

# Add labels and title
plt.xlabel('vn')
plt.ylabel('Frequency')
plt.title('Distribution of vn')

# Show the plot
plt.show()

#%%
# Create the input feature matrix
X = np.column_stack((x, y, z))
X = scaler.fit_transform(X)

#Transform X_scaled back to the original scale
X_original = scaler.inverse_transform(X)
#endregion
#endregion

#%% MGP-PCA (M)

# Combine X and vn into a single matrix
data_PCA = np.column_stack((X, vn))

# Create an instance of PCA
pca = PCA()

# Fit the PCA model to the data
pca.fit(data_PCA)

# Transform the data using the PCA model
transformed_data = pca.transform(data_PCA)

# Extract the transformed X and vn
transformed_X = transformed_data[:, :3]  # Assuming X has 3 columns
transformed_vn = transformed_data[:, 3]  # Assuming vn is the last column

X = transformed_X
vn = transformed_vn

# Print the explained variance ratio
print("Explained Variance Ratio:", pca.explained_variance_ratio_)


#%% Cross validation find best gp
#region GP model
score = -1000
for size in cv_size: 
    # Splitting the dataset into training and test sets
    # Random state is currently fixed
    #seed = random.randint(1, 100)
    seed = 20 
    X_train, X_test, vn_train, vn_test = train_test_split(X, vn, test_size=size, random_state= seed)

#%% Gaussian_model Best is Rational quadratic

#Rational Quadratic best: size = 0.1,seed = 20, alpha = 1, length_scale = 1
# for K_num in [10]:
#     for l in [1]:
#         for alpha_num in [1]:
            
#             ker = RationalQuadratic(length_scale= l ,alpha=alpha_num,length_scale_bounds=(1e-100,1e100))  # You can choose other kernels as well
            
#             gp = GaussianProcessRegressor(kernel=ker,
#                                             optimizer='fmin_l_bfgs_b',  # Use L-BFGS-B optimizer
#                                             alpha = noise_level,
#                                             n_restarts_optimizer=3) 
            
#             #gp.fit(X_train, vn_train_noi)
#             gp.fit(X_train,vn_train)
            
#             print(f"The model is {gp}")
#             print(f"The K-Fold is {K_num}")
#             print(f"The size is {size}")
            
#             #neg_mean_squared_error --> super smart way, since most optimization target at maximizing the score, 
#             #but we want to minimize the mean_square error, which is the same as maximizing the negative mean squared error
            
#             temp_score = np.mean(cross_val_score(gp, X_train, vn_train, cv= K_num, scoring='neg_mean_squared_error'))
            
#             if score < temp_score:
                    
#                 best_gp = gp
#                 score = temp_score
#                 K = K_num
for K_num in K_num_list:
    for l in length_scale_list:
        for alpha_num in alpha_num_list:
            
            ker = RationalQuadratic(length_scale= l ,alpha=alpha_num,length_scale_bounds=(1e-100,1e100))  # You can choose other kernels as well
            
            gp = GaussianProcessRegressor(kernel=ker,
                                            alpha= noise_level,
                                            optimizer='fmin_l_bfgs_b',  # Use L-BFGS-B optimizer
                                            n_restarts_optimizer=3) 
            
            kf = KFold(n_splits=K_num)
            temp_scores = []
            
            for train_index, test_index in kf.split(X_train):
                X_train_kf, X_test_kf = X_train[train_index], X_train[test_index]
                vn_train_kf, vn_test_kf = vn_train[train_index], vn_train[test_index]
                
                gp.fit(X_train_kf, vn_train_kf)
                
                temp_scores.append(gp.log_marginal_likelihood())
            
            temp_score = np.mean(temp_scores)
            
            if score < temp_score:
                    
                best_gp = gp
                score = temp_score
                K = K_num
                print(temp_scores)

# Matern R-squared always negative

    # for K_num in [10]:
    #     for l in [1]:
    #         for nu_num in [2.5]:
                
    #             ker = Matern(length_scale= l , nu= nu_num,length_scale_bounds=(1e-100,1e100))  # You can choose other kernels as well
                
    #             gp = GaussianProcessRegressor(kernel=ker,
    #                                           optimizer='fmin_l_bfgs_b',  # Use L-BFGS-B optimizer
    #                                           n_restarts_optimizer=3) 
                
    #             gp.fit(X_train, vn_train_noi)
                
    #             #neg_mean_squared_error --> super smart way, since most optimization target at maximizing the score, 
    #             #but we want to minimize the mean_square error, which is the same as maximizing the negative mean squared error
    #             temp_score = np.mean(cross_val_score(gp, X_train, vn_train_noi, cv= K_num, scoring='neg_mean_squared_error'))
                
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
        
#         #gp.fit(X_train, vn_train_noi)
#         gp.fit(X_train,vn_train_noi)
        
#         print(f"The model is {gp}")
#         print(f"The K-Fold is {K_num}")
#         print(f"The size is {size}")
#         print(f"The sig is {sig}")
       
        
#         #neg_mean_squared_error --> super smart way, since most optimization target at maximizing the score, 
#         #but we want to minimize the mean_square error, which is the same as maximizing the negative mean squared error
        
#         temp_score = np.mean(cross_val_score(gp, X_train, vn_train_noi, cv= K_num, scoring='r2'))
        
#         if score < temp_score:
            
#             best_gp = gp
#             score = temp_score
#             K = K_num

# DotProduct kernel not good

# White Kernel not good

# Compound Kernel

#%%#By this section we see that gp is perfectly matched if adding a noise

#vn_pred,sigma = gp.predict(X_train,return_std=True)

#mse = mean_squared_error(vn_train, vn_pred)
#mse_noi = mean_squared_error(vn_train_noi, vn_pred)

#print(f"The mse is {mse}")
#print(f"The mse_noi is {mse_noi}")

#%% Random Forest Prediction

# for K_num in [10]:

#     # Initialize the Random Forest Regressor model
#     rf_model = RandomForestRegressor(random_state=seed,n_estimators = 100,criterion='absolute_error')
    
#     vn_train = vn_train.ravel()
    
#     # Train the Random Forest model on the training data
#     rf_model.fit(X_train, vn_train)

#     temp_score = np.mean(cross_val_score(rf_model, X_train, vn_train, cv= K_num, scoring='r2'))

#     if score < temp_score:
                
#         best_rf = rf_model
#         score = temp_score
#         K = K_num
#endregion
#%% Take into test set
#region GP model predict
best_gp.fit(X_train, vn_train)
vn_pred,sigma = best_gp.predict(X_test,return_std=True)
#vn_pred = best_rf.predict(X_test)

mse = mean_squared_error(vn_test, vn_pred)
r2 = r2_score(vn_test, vn_pred)
log= gp.log_marginal_likelihood()

print(f"The mse is {mse}, this is the mse between vn_pred and vn_test")
print(f"The R-squared is {r2}")

print(f"The seed is {seed}")
print(f"The best_model is {gp}")
print(f"The log_marginal_likelihood is {log}")
#print(f"The best_model is {best_rf}")
print(f"The best_cross_val K-Fold is {K}")
print(f"The best_size is {size}")
print(f"The noise is noise_level = {noise_level}")

plt.hist(vn, bins=20, edgecolor='black')
plt.xlabel('Log 10 Virus Number')
plt.ylabel('Frequency')
plt.title('Distribution of Virus Number')
plt.show()
count = np.sum(vn_test > 0.2)
print("Number of observations in vn where vn_test > 0.5:", count)
#endregion 
#%% Plot

# # Create a figure with two subplots
# fig, axs = plt.subplots(1, 2, figsize=(12, 6))

# # Scatter plot for vn
# axs[0].scatter(range(len(vn_test)), vn_test, color='blue', marker='o')
# axs[0].set_xlabel('Index')
# axs[0].set_ylabel('vn')
# axs[0].set_title('Scatter Plot of vn_test')

# # Scatter plot for vn_pred

# axs[1].scatter(range(len(vn_pred)), vn_pred, color='blue', marker='o')
# axs[1].set_xlabel('Index')
# axs[1].set_ylabel('vn_pred')
# axs[1].set_title('Scatter Plot of vn_pred')

# # Add a main title
# plt.suptitle(f"Test pred on test data", fontsize=16)

# plt.tight_layout()
# plt.show()

#%%

fig, ax = plt.subplots() # 创建图实例
ax.errorbar(range(len(vn_pred)), vn_pred.ravel(), yerr=sigma,ecolor='r', color='b', fmt='o', label='Uncertainty')
ax.scatter(range(len(vn_test)), vn_test, color='g', marker='o', label='vn_test')

ax.set_xlabel('Index')
ax.set_ylabel('vn_pred')
ax.set_title('Predicted vn with Uncertainty')
ax.legend()
plt.show()

#%% Test pred on train data for rf

# Train final model with best kernel and hyperparameters
# best_rf.fit(X_train, vn_train)

# vn_pred= best_rf.predict(X_train)
# mse = mean_squared_error(vn_train, vn_pred)
# print("mse between vn_test and vn_predict(on X_train) is" ,mse)

# fig, axs = plt.subplots(1, 2, figsize=(12, 6))
# # Scatter plot for vn
# axs[0].scatter(range(len(vn_train)), vn_train, color='blue', marker='o')
# axs[0].set_xlabel('Index')
# axs[0].set_ylabel('vn')
# axs[0].set_title('Scatter Plot of vn_train')

# # Scatter plot for vn_pred
# axs[1].scatter(range(len(vn_pred)), vn_pred, color='blue', marker='o')
# axs[1].set_xlabel('Index')
# axs[1].set_ylabel('vn_pred')
# axs[1].set_title('Scatter Plot of vn_pred')

# # Add a main title
# plt.suptitle(f"Test pred on train data", fontsize=16)

# plt.tight_layout()
# plt.show()

#%% First time 
#region Delete non-predictable data points
vn_pred = vn_pred.reshape(sigma.shape)
vn_test = vn_test.reshape(sigma.shape)

#Find indices of points in vn_test
indices_to_remove = np.where(np.logical_or(vn_pred + sigma < vn_test, vn_pred - sigma > vn_test))[0]
#find indices of points in vn
matching_indices = np.where(np.isin(vn, vn_test[indices_to_remove]))[0]

#delete those unpredictable points
vn_filtered = np.delete(vn, matching_indices)
X_filtered = np.delete(X, matching_indices, axis=0)
print(len(X_filtered))

# %% iterate canceling points

for i in range(1,10):

    #! 这个900的话是删除点的阈值，如果说删除点删太多了就自动停止
    #! 不得不吐槽一下我其实觉得这种删点也挺不合理的。。。毕竟本来data就这么点
    if len(X_filtered) < 900:
        break

    for size in cv_size:
        # Splitting the dataset into training and test sets
        # Random state is currently fixed
        #seed = random.randint(1, 100)
        seed = 20
        X_train, X_test, vn_train, vn_test = train_test_split(X_filtered, vn_filtered, test_size=size, random_state= seed)

    #refit model
    best_gp.fit(X_train, vn_train)
    vn_pred,sigma = best_gp.predict(X_test,return_std=True)

    #evaluate model
    mse = mean_squared_error(vn_test, vn_pred)
    r2 = r2_score(vn_test, vn_pred)
    log= best_gp.log_marginal_likelihood()

    print(f"The mse is {mse}, this is the mse between vn_pred and vn_test")
    print(f"The R-squared is {r2}")
    print(f"The log_marginal_likelihood is {log}")

    #draw figure
    fig, ax = plt.subplots() # 创建图实例
    ax.errorbar(range(len(vn_pred)), vn_pred.ravel(), yerr=sigma,ecolor='r', color='b', fmt='o', label='Uncertainty')
    ax.scatter(range(len(vn_test)), vn_test, color='g', marker='o', label='vn_test')

    ax.set_xlabel('Index')
    ax.set_ylabel('vn_pred')
    ax.set_title('Predicted vn with Uncertainty')
    ax.legend()
    plt.show()

    #remove points
    vn_pred = vn_pred.reshape(sigma.shape)
    vn_test = vn_test.reshape(sigma.shape)

    #Find indices of points in vn_test
    indices_to_remove = np.where(np.logical_or(vn_pred + sigma < vn_test, vn_pred - sigma > vn_test))[0]

    #find indices of points in vn
    matching_indices = np.where(np.isin(vn_filtered, vn_test[indices_to_remove]))[0]

    #delete those unpredictable points
    vn_filtered = np.delete(vn_filtered, matching_indices)
    X_filtered = np.delete(X_filtered, matching_indices, axis=0)
    print(len(X_filtered))
  
#endregion
# %%
