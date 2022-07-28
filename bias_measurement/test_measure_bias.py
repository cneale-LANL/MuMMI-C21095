import numpy as np
import os
import pickle as pickle
import sys

import matplotlib.pyplot as plt
from matplotlib import rcParams,cm
rcParams.update({'figure.autolayout': True})



'''
The following code measures bias by fitting a regression line between two variables X and Y.

We assume that Y = AX + B

If Y is an unbiased estimator of X, then A = 1 and B = 0.

If either A =/= 1 or B =/= 0 is a significant result, then we predict that Y is
a biased estimator of X.

The bias calculation utilities can be found in bias_measurement_utility file

'''
from bias_measurement_utility import p_test


'''
Load data

The npz archive below contains three files:

Required
1. X : m x m x n numpy array
2. Y : m x m x n numpy array
Optional
3. labels: a list containing n variable labels

'''
f = np.load("data_for_bias_calculation.npz")

X = f['X']
Y = f['Y']
labels = f['labels']
nVar = X.shape[2]

'''

Measure p-values using functions from bias_measurement_utility

'''
p_th = 0.001

print("\n")
print("="*20)
print("\n")

print ("Using p = %0.1e to test significance."%p_th)
print("\n")

for l in range(nVar):
    x = X[:,:,l].flatten()
    y = Y[:,:,l].flatten()
    label = labels[l]
    p_test(x,y,label,p_th)

print("\n")
print("="*20)
print("\n")
