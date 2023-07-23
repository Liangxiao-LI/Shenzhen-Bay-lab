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

#Sample plot

cm = plt.get_cmap("RdYlGn")

x = np.random.rand(30)
y = np.random.rand(30)
z = np.random.rand(30)
col = np.arange(30)

# 2D Plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(x, y, s=10, c=col, marker='o')  

# 3D Plot
fig = plt.figure()
ax3D = fig.add_subplot(111, projection='3d')
p3d = ax3D.scatter(x, y, z, s=30, c=col, marker='o')                                                                                

fig.colorbar(p3d)


plt.show()


#Implementation plot

df = pd.read_excel('/Users/ryan/Documents/GitHub/Shenzhen Bay lab/spike_protein_3structures_markerd_with_variants_from_cncb/6vxx_variants.xls')

sub_df = df.iloc[:972, 5:10]


# Extract the x, y, z, and m columns from the subset DataFrame
x = sub_df.iloc[:972,0]
y = sub_df.iloc[:972,1]
z = sub_df.iloc[:972,2]
m = sub_df.iloc[:972,4]
col = np.arange(972)

print(len(x))

"""
# Define the threshold value for m
m_min = 0.001

# Create a Boolean mask that selects only the rows of sub_df where m is greater than m_min
mask = sub_df['virusPercent'] > m_min

# Use the mask to select the corresponding values of x, y, z, and m to plot
x = sub_df.loc[mask, 'X']
y = sub_df.loc[mask, 'Y']
z = sub_df.loc[mask, 'Z']
m = sub_df.loc[mask, 'virusPercent']

print(len(x))


# Define the range of values for the m column
m_min = m.min()
m_max = m.max()

print("m_min is",m_min)  
print("m_max is",m_max)  

# Create a custom color map that maps low values to light blue and high values to dark blue
#cmap = colors.ListedColormap(['lightblue', 'darkblue'])
#norm = colors.Normalize(vmin=m_min, vmax=m_max)

# Map the values of m to colors using the custom color map
#col = cmap(norm(m))

# Create a color map based on the values in the m column
cmap = plt.cm.get_cmap('viridis')
col = cmap(m)


# Create a 3D scatter plot with color based on m
fig = plt.figure()
ax3D = fig.add_subplot(111, projection='3d')
p3d = ax3D.scatter(x, y, z, s=30, c=col, marker='o')

# Add a title to the plot
ax3D.set_title('virusPercent >= 0.001')

# Add a color bar to the plot
fig.colorbar(p3d)

# Show the plot
plt.show()

ax3D.view_init(elev=30, azim=45)
"""

dfs = {}

m_Min = [0,0.001,0.01,0.02,0.05,0.25,0.5]

for m_min in m_Min:
    
    
    # Create a Boolean mask that selects only the rows of sub_df where m is greater than m_min
    mask = sub_df['virusPercent'] > m_min
    
    # Use the mask to select the corresponding values of x, y, z, and m to plot
    x = sub_df.loc[mask, 'X']
    y = sub_df.loc[mask, 'Y']
    z = sub_df.loc[mask, 'Z']
    m = sub_df.loc[mask, 'virusPercent']
    
    print(len(x))
    
    # Create a data frame from the selected data points and store it in the dictionary
    dfs[m_min] = pd.DataFrame({'X': x, 'Y': y, 'Z': z, 'virusPercent': m})
    
    # Define the range of values for the m column
    m_min = m.min()
    m_max = m.max()
    
    print("m_min is",m_min)  
    print("m_max is",m_max)  
    
    # Create a custom color map that maps low values to light blue and high values to dark blue
    #cmap = colors.ListedColormap(['lightblue', 'darkblue'])
    #norm = colors.Normalize(vmin=m_min, vmax=m_max)
    
    # Map the values of m to colors using the custom color map
    #col = cmap(norm(m))
    
    # Create a color map based on the values in the m column
    cmap = plt.cm.get_cmap('viridis')
    col = cmap(m)
    
    
    # Create a 3D scatter plot with color based on m
    fig = plt.figure()
    ax3D = fig.add_subplot(111, projection='3d')
    p3d = ax3D.scatter(x, y, z, s=30, c=col, marker='o')
    
    # Add a title to the plot
    ax3D.set_title(f'virusPercent >= {m_min:.3f}')
    
    # Add a color bar to the plot
    fig.colorbar(p3d)
    
    # Show the plot
    plt.show()
    
    # Access the data frames by their keys
    for m_min, df in dfs.items():
        print(f'Data frame for virusPercent >= {m_min}:')
        print(df)
        print(f'the number of observation is {len(df)}:')