# -*- coding: utf-8 -*-

# Chebyshev polynomial

import math
import numpy as np
from initializeChebyshevApproximator import init

def coef(f, D, N, M, a, b):
    X = init(D,N,M)[0]
    T = init(D,N,M)[1]
    y = np.zeros(pow(M,D))
        
    for i in range(len(y)):
        ksi = a + (X[i] + 1) * (b - a) / 2
        y[i] = f(ksi)
    
    theta = np.linalg.lstsq(T,y)[0]
        
    return theta