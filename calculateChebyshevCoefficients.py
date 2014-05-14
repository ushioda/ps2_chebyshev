# -*- coding: utf-8 -*-

# Chebyshev polynomial

import math
import numpy as np
from initializeChebyshevApproximator import init
import statsmodels.api as sm

def coef(f, D, N, M, a, b):
    X = init(D,N,M)[0]
    T = init(D,N,M)[1]
    y = np.empty(pow(M,D))
        
    for xi, yi in zip(X,y):
        ksi = a + (xi + 1) * (b - a) / 2
        yi = f(ksi)
    
    T = T.astype(float32)   
    y = y.astype(float32) 
    
    theta = np.linalg.lstsq(T,y)[0]
        
    return theta