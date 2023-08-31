#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  6 15:18:31 2023

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
import random
from sklearn.gaussian_process.kernels import RBF, Matern, RationalQuadratic, ExpSineSquared, DotProduct
from sklearn.model_selection import train_test_split



#%%define add noise function

# x is my training data
# mu is the mean
# std is the standard deviation
mu=0.0
std = 1
def gaussian_noise(x,mu,std):
    noise = np.random.normal(mu, std, size = x.shape)
    x_noisy = x + noise
    return x_noisy 

#%%
#Implementation plot

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
vp = np.log10(vp)
col = np.arange(972)

x = np.array(x)
y = np.array(y)
z = np.array(z)

vp = np.array(vp)
vp_noisy = gaussian_noise(vp,mu,std)

# Create the input feature matrix
X = np.column_stack((x, y, z))

#%%
# Splitting the dataset into training and test sets
X_train, X_test, vp_train_noi, vp_test_noi = train_test_split(X, vp_noisy, test_size=0.2, random_state=42)

X_train, X_test, vp_train, vp_test = train_test_split(X, vp, test_size=0.2, random_state=42)



#%% 纯手工K-Fold

def divide_integers_into_sets(x, K):
    num_integers_per_set = math.ceil(x / K)
    all_integers = list(range(x))
    random.shuffle(all_integers)  # Shuffle the integers to ensure randomness
    integer_sets = [all_integers[i:i+num_integers_per_set] for i in range(0, x, num_integers_per_set)]
    return integer_sets

k = 8
lowest_mse = 100
#Divide into K sets
K_Sets = divide_integers_into_sets(len(X_train),k)

test_index = K_Sets[Knum]
train_index = [num for num in range(len(x) ) if num not in test_index]


for Knum in range(0,k):
    #This part tries to find the best model
    for i in [1,10]:
        
        #Fit a gp
        best_kernel = RBF(length_scale=i)  # You can choose other kernels as well
    
        gp = GaussianProcessRegressor(kernel=best_kernel,
                                      optimizer='fmin_l_bfgs_b',  # Use L-BFGS-B optimizer
                                      n_restarts_optimizer=3,   # Number of optimizer restarts
                                      alpha=0.1,
                                      random_state=0) 
        
        gp.fit(X_train, vp_train_noi)
    
        #K-Fold cross validation
        #predict it on the 
    
        vp_pred,sigma = gp.predict(X_test,return_std=True)
        mse = mean_squared_error(vp_test, vp_pred)
        
        print("mse is",mse)
        
        if lowest_mse > mse:
            final_gp = gp
            lowest_mse = mse
            
    
    for i in [1,10]:
        for nu_num in [0.5,1.5]:
        # Train final model with best kernel and hyperparameters
            kernel = Matern(length_scale=i, nu= nu_num)  # You can adjust nu as needed
        
            gp = GaussianProcessRegressor(kernel=kernel,
                                          optimizer='fmin_l_bfgs_b',  # Use L-BFGS-B optimizer
                                          n_restarts_optimizer=3,   # Number of optimizer restarts
                                          random_state=0)            # For reproducibility
            
            gp.fit(X_train, vp_train_noi)
        
            vp_pred, sigma = gp.predict(X_test, return_std=True)
            mse = mean_squared_error(vp_test, vp_pred)
            
            print(f"mse for length_scale {i},nu {nu_num} is {mse}")
            
            if mse < lowest_mse:
                final_gp = gp
                lowest_mse = mse
                best_kernel = kernel
    
    #Evaluate the model by K-Fold

    
    #This part assigns values
    x_train, x_test = x[train_index], x[test_index]
    y_train, y_test = y[train_index], y[test_index]
    z_train, z_test = z[train_index], z[test_index]
    vp_train,vp_train_noi, vp_test = vp[train_index],vp_noisy[train_index], vp[test_index]

    X_train = np.column_stack((x_train, y_train, z_train))
    X_test = np.column_stack((x_test, y_test, z_test))


            
    
print("mse between vp_test and vp_predict(X_test) is" ,lowest_mse)

#%% 
#create a scatter plot of vp_pred and vp_test

# Assuming vp_pred and vp_test are your predicted and test values
x_row_numbers = range(len(X_test))

plt.scatter(x_row_numbers, vp_test, color='blue', marker='o', label='vp_test')
plt.scatter(x_row_numbers, vp_pred, color='red', marker='x', label='vp_pred')

# Add labels and title
plt.xlabel('Row Number of X_test')
plt.ylabel('Values')
plt.title('Scatter Plot: vp_test and vp_pred')
plt.legend()

plt.show()

#%%
import matplotlib.pyplot as plt

# Create a figure with two subplots
fig, axs = plt.subplots(1, 2, figsize=(12, 6))

# Scatter plot for vp
axs[0].scatter(range(len(vp_test)), vp_test, color='blue', marker='o')
axs[0].set_xlabel('Index')
axs[0].set_ylabel('vp')
axs[0].set_title('Scatter Plot of vp_test')

# Scatter plot for vp_pred


vp_pred_up = [ai + bi for ai, bi in zip(vp_pred, sigma)]

vp_pred_down = [ai - bi for ai, bi in zip(vp_pred, sigma)]


axs[1].scatter(range(len(vp_pred)), vp_pred, color='blue', marker='o')
axs[1].plot(range(len(vp_pred)), vp_pred_up, color='red', label='vp_pred_up')
axs[1].plot(range(len(vp_pred)), vp_pred_down, color='green', label='vp_pred_down')
axs[1].set_xlabel('Index')
axs[1].set_ylabel('vp_pred')
axs[1].set_title('Scatter Plot of vp_pred')

# Add a main title
plt.suptitle(f"Test pred on test data", fontsize=16)

plt.tight_layout()
plt.show()

#%%

import matplotlib.pyplot as plt

# Split vp_pred into segments
segment_size = 10  # Adjust this based on your preference (should be multiple of 2)
num_figs = segment_size // 2

for fig_idx in range(num_figs):
    # Calculate the indices for the current figure's segments
    start_idx = fig_idx * 2 * segment_size
    end_idx = (fig_idx * 2 + 2) * segment_size
    segment_vp_pred = vp_pred[start_idx:end_idx]
    segment_sigma = sigma[start_idx:end_idx]
    segment_vp_test = vp_test[start_idx:end_idx]

    # Create a figure with two subplots
    fig, axs = plt.subplots(2, 1, figsize=(8, 12))

    for segment_idx in range(2):
        segment_start = segment_idx * segment_size
        segment_end = (segment_idx + 1) * segment_size
        segment_vp_pred = vp_pred[start_idx + segment_start : start_idx + segment_end]
        segment_sigma = sigma[start_idx + segment_start : start_idx + segment_end]
        segment_vp_test = vp_test[start_idx + segment_start : start_idx + segment_end]  

        # Scatter plot for vp_pred
        vp_pred_up = [ai + bi for ai, bi in zip(segment_vp_pred, segment_sigma)]
        vp_pred_down = [ai - bi for ai, bi in zip(segment_vp_pred, segment_sigma)]

        axs[segment_idx].scatter(range(len(segment_vp_pred)), segment_vp_pred, color='blue', marker='o', label='vp_pred')
        axs[segment_idx].scatter(range(len(segment_vp_test)), segment_vp_test, color='red', marker='x', label='vp_test')
        axs[segment_idx].plot(range(len(segment_vp_pred)), vp_pred_up, color='red', label='vp_pred_up')
        axs[segment_idx].plot(range(len(segment_vp_pred)), vp_pred_down, color='green', label='vp_pred_down')
        axs[segment_idx].set_xlabel('Index')
        axs[segment_idx].set_ylabel('vp_pred')
        axs[segment_idx].set_title(f'Scatter Plot of vp_pred (Segment {segment_idx+1})')
        axs[segment_idx].legend()

    # Add a main title to the figure
    fig.suptitle(f"Test pred on test data (Figure {fig_idx+1})", fontsize=16)

    plt.tight_layout()
    plt.show()









#%%
# Create a scatter plot for the actual data X and vp
fig = plt.figure(figsize=(12, 6))

# Plot for actual data
ax1 = fig.add_subplot(121, projection='3d')
sc1 = ax1.scatter(X_test[:, 0], X_test[:, 1], X_test[:, 2], c=vp_test, cmap='viridis', s=50)
cbar1 = plt.colorbar(sc1)
cbar1.set_label('Actual vp')

ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.set_zlabel('Z')
ax1.set_title('Actual Data - X_test,vp_test')

# Create a scatter plot for the predicted data X_pred and vp_pred
ax2 = fig.add_subplot(122, projection='3d')
sc2 = ax2.scatter(X_test[:, 0], X_test[:, 1], X_test[:, 2], c=vp_pred, cmap='viridis', s=50)
cbar2 = plt.colorbar(sc2)
cbar2.set_label('Predicted vp_pred')

ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_zlabel('Z')
ax2.set_title('Predicted Data - X_test,vp_predict')

# Add a main title
plt.suptitle(f"3D scatter plot comparing predict and test", fontsize=16)


plt.tight_layout()
plt.show()

#%%

# Create a scatter plot for sigma
plt.scatter(range(len(sigma)), sigma, color='blue', marker='o')

# Add labels and title
plt.xlabel('Index')
plt.ylabel('Sigma')
plt.title('Scatter Plot of Sigma')

plt.show()

#%% Test pred on train data

# Train final model with best kernel and hyperparameters
final_gp = GaussianProcessRegressor(kernel=best_kernel)
final_gp.fit(X_train, vp_train_noi)

vp_pred,sigma = final_gp.predict(X_train,return_std=True)
mse = mean_squared_error(vp_train, vp_pred)
print("mse between vp_test and vp_predict(on X_train) is" ,mse)

fig, axs = plt.subplots(1, 2, figsize=(12, 6))
# Scatter plot for vp
axs[0].scatter(range(len(vp_train)), vp_train, color='blue', marker='o')
axs[0].set_xlabel('Index')
axs[0].set_ylabel('vp')
axs[0].set_title('Scatter Plot of vp_train')

# Scatter plot for vp_pred
axs[1].scatter(range(len(vp_pred)), vp_pred, color='blue', marker='o')
axs[1].set_xlabel('Index')
axs[1].set_ylabel('vp_pred')
axs[1].set_title('Scatter Plot of vp_pred')

# Add a main title
plt.suptitle(f"Test pred on train data", fontsize=16)

plt.tight_layout()
plt.show()