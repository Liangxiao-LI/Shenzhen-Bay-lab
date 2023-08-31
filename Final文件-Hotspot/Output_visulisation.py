
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
# %%



df = pd.read_excel('/Users/ryan/Documents/GitHub/Shenzhen Bay lab/spike_protein_3structures_markerd_with_variants_from_cncb/6vxx_variants.xls')

for i in range(1, 21):
    filename = f'/Users/ryan/Documents/GitHub/Shenzhen Bay lab/Final文件-Hotspot/output/output_{i}.csv'
    df_output = pd.read_csv(filename)
    merged_df = pd.merge(df, df_output, on='atomSerialNumber', how='outer')


    sub_df = merged_df.iloc[:972, 5:10]
    cluster = merged_df.iloc[:972, 15]
    sub_df_CYS = merged_df.iloc[:972, 5:10]

    # Extract the x, y, z, and m columns from the subset DataFrame
    x = sub_df.iloc[:972,0]
    y = sub_df.iloc[:972,1]
    z = sub_df.iloc[:972,2]
    vn = sub_df.iloc[:972,3]
    vp = sub_df.iloc[:972,4]
    vp = sub_df.iloc[:972,]
    is_clust = cluster.iloc[:972,]


    # Create a 3D scatter plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot the points with blue color
    #ax.scatter(x, y, z, c='blue')
    ax.scatter(x, y, z, c = is_clust)
    # Set the figure title
    ax.set_title(f'Output File: output_{i}.csv')

    # Show the plot
    plt.show()

# %%


# %%
