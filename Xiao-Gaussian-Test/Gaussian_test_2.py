#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  6 13:57:18 2023

@author: ryan
"""

from sklearn.gaussian_process import GaussianProcessClassifier

import numpy as np
import matplotlib.pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
from sklearn.gaussian_process.kernels import ConstantKernel as C



# Generate data
rng = np.random.default_rng(0)
X = rng.uniform(-5, 5, 25)[:, np.newaxis]
y = np.sin(X).ravel()
y += rng.normal(0, 0.1, y.size)

# Specify Gaussian Processes with fixed and optimized hyperparameters
gp_fix = GaussianProcessRegressor(kernel=1.0 * RBF(length_scale=1.0),optimizer=None)
gp_fix.fit(X, y)  
gp_opt = GaussianProcessRegressor(kernel=1.0 * RBF(length_scale=1.0))
gp_opt.fit(X, y)


# Make predictions
X_ = np.linspace(-5, 5, 100)[:, np.newaxis]
y_fix, sigma_fix = gp_fix.predict(X_, return_std=True)
y_opt, sigma_opt = gp_opt.predict(X_, return_std=True)


# Plot results
plt.figure()
plt.plot(X, y, "r.", markersize=10, label="Observations")
plt.plot(X_, y_fix, "b-", label="Prediction (fixed kernel)")
plt.plot(X_, y_opt, "g-", label="Prediction (optimized kernel)")
plt.fill(
    np.concatenate([X_, X_[::-1]]),
    np.concatenate([y_fix - 1.9600 * sigma_fix,
                    (y_fix + 1.9600 * sigma_fix)[::-1]]),
    alpha=.5,
    fc="b",
    ec="None",
    label="95% confidence interval (fixed kernel)",
)
plt.fill(
    np.concatenate([X_, X_[::-1]]),
    np.concatenate([y_opt - 1.9600 * sigma_opt,
                    (y_opt + 1.9600 * sigma_opt)[::-1]]),
    alpha=.5,
    fc="g",
    ec="None",
    label="95% confidence interval (optimized kernel)",
)
plt.xlabel("Feature")
plt.ylabel("Target")
plt.xlim(-5, 5)
plt.ylim(-3, 3)
plt.legend(loc="best")










# Plot LML landscape
plt.figure()
theta0 = np.logspace(-1, 3, 30)
theta1 = np.logspace(-1, 3, 29)
Theta0, Theta1 = np.meshgrid(theta0, theta1)
LML = [
    [
        gp_opt.log_marginal_likelihood(np.log([Theta0[i, j],
                                               Theta1[i, j]]))
        for i in range(Theta0.shape[0])
    ]
    for j in range(Theta0.shape[1])
]
LML = np.array(LML).T
plt.contour(Theta0, Theta1, LML)
plt.scatter(
    gp_opt.kernel_.theta[0],
    gp_opt.kernel_.theta[1],
    c="r",
    s=50,
    zorder=10,
    edgecolors=(0, 0, 0),
)
plt.plot(
    [gp_fix.kernel_.theta[0]],
    [gp_fix.kernel_.theta[1]],
    "bo",
    ms=10,
)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Length scale")
plt.ylabel("Noise level")
plt.title("Log-marginal-likelihood")






# Plot LML as a function of length scale
plt.figure()
plt.plot(
    gp_opt.kernel_.theta[0],
    gp_opt.log_marginal_likelihood(gp_opt.kernel_.theta),
    "bo",
    ms=10,
)
plt.plot(
    gp_fix.kernel_.theta[0],
    gp_fix.log_marginal_likelihood(gp_fix.kernel_.theta),
    "ro",
    ms=10,
)
plt.xlabel("Length scale")
plt.ylabel("Log-marginal-likelihood")
plt.title("Log-marginal-likelihood as\
a function of length scale")
  
plt.show()