# -*- coding: utf-8 -*-

############################
#
# Quantitative Marketing PS2
# Yu Ushioda
#
############################

import numpy as np
import math
from random import uniform
from evaluateChebyshev import approx
from calculateChebyshevCoefficients import coef
from calculateChebyshevPolynomials import chebypol

s_size = 10
N = 3
M = 7
a = np.array([-5, -2, -3])
b = np.array([2, 4, 3])
D = len(a)

# prepare samples

X = np.zeros(D * s_size)
X.shape = (s_size, D)

for i in range(s_size):
    sample = np.zeros(D)
    for j in range(D):
        sample[j] = uniform(a[j], b[j])
    X[i] = sample

# define functions

def f(z):
    y = z[0] * pow(z[2],3) + z[1] * z[2] + pow(z[0],2) * z[1] * pow(z[2], 2)
    return y
    
def g(z):
    y = z[0] * math.log(5 + z[1] * z[2])
    return y
    
def h(z):
    y = pow(z[0],2) * math.cos(z[1]) * math.exp(z[2])
    return y
    
# calculate actual values

y_f = np.zeros(s_size)
for i in range(s_size):
    y_f[i] = f(X[i])
    
y_g = np.zeros(s_size)
for i in range(s_size):
    y_g[i] = g(X[i])

y_h = np.zeros(s_size)
for i in range(s_size):
    y_h[i] = h(X[i])
    
# calculate Chebyshev-approximated values
# since we don't want to calculate theta multiple times,
# we won't use evaluateChebyshev.py

"""

theta_f = coef(f,D,N,M,a,b)
    
y_f_hat = np.zeros(s_size)

for i in range(s_size):
    norm_x = - 1 + 2 * (X[i] - a) / (b - a)
    B = chebypol(N, norm_x)
    
    T = B
    for i in range(D - 1):
        T = np.kron(T, B)

    y_f_hat[i] = np.dot(T, theta_f)

"""