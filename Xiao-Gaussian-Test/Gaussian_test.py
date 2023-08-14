#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 13:14:16 2023

@author: ryan
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.colors as colors
import pandas as pd
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF


#Implementation plot

df = pd.read_excel('/Users/ryan/Documents/GitHub/Shenzhen Bay lab/spike_protein_3structures_markerd_with_variants_from_cncb/6vxx_variants.xls')

sub_df = df.iloc[:972, 5:10]


# Extract the x, y, z, and m columns from the subset DataFrame
x = sub_df.iloc[:972,0]
y = sub_df.iloc[:972,1]
z = sub_df.iloc[:972,2]
m = sub_df.iloc[:972,4]
col = np.arange(972)



x = np.array(x)
y = np.array(y)
z = np.array(z)
m = np.array(m)

# Create the input feature matrix
X = np.column_stack((x, y, z))

# Define the Gaussian process regression model with the RBF kernel
kernel = RBF()
gp = GaussianProcessRegressor(kernel=kernel)

# Fit the model to your data
gp.fit(X, m)

# Predict viruspercentage for new input points
new_x = np.array([100, 110, 120])  # Replace with your new input values
predicted_viruspercentage, std = gp.predict([new_x], return_std=True)

# Print the predicted viruspercentage and its standard deviation
print("Predicted viruspercentage:", predicted_viruspercentage[0])
print("Standard deviation:", std[0])

#Visulisation

# Generate a meshgrid for plotting the surface
x_range = np.linspace(min(x), max(x), 100)
y_range = np.linspace(min(y), max(y), 100)
z_range = np.linspace(min(z), max(z), 100)
X_mesh, Y_mesh, Z_mesh = np.meshgrid(x_range, y_range, z_range)


# Flatten the meshgrid for prediction
X_pred = np.column_stack((X_mesh.flatten(), Y_mesh.flatten(), Z_mesh.flatten()))

# Predict the mean and standard deviation for the meshgrid points
m_pred, std_pred = gp.predict(X_pred, return_std=True)

# Reshape the predictions to match the meshgrid shape
M_pred = m_pred.reshape(X_mesh.shape)
Std_pred = std_pred.reshape(X_mesh.shape)

# Plot the mean function
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X_mesh, Y_mesh, Z_mesh, M_pred, cmap='jet', alpha=0.5)

# Plot confidence intervals as transparent surfaces
lower_bound = M_pred - 2 * Std_pred
upper_bound = M_pred + 2 * Std_pred
ax.plot_surface(X_mesh, Y_mesh, Z_mesh, lower_bound, color='gray', alpha=0.3)
ax.plot_surface(X_mesh, Y_mesh, Z_mesh, upper_bound, color='gray', alpha=0.3)

# Customize the plot
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Gaussian Process Model')

# Show the plot
plt.show()