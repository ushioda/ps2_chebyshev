# -*- coding: utf-8 -*-

import math
import numpy as np
from calculateChebyshevCoefficients import coef
from calculateChebyshevPolynomials import chebypol

def approx(x,f,D,N,M,a,b):
    
    theta = coef(f,D,N,M,a,b)
    norm_x = - 1 + 2 * (x - a) / (b - a)
    B = chebypol(N, norm_x)
    
    T = B
    for i in range(D - 1):
        T = np.kron(B, T)
    
    if D == 1:
        j = 1
    else:
        k = 1
        A = []    
        while k <= D - 2:
            A.append(k * pow(D, D - k - 1))
            k = k + 1
        j = sum(A) + D
        
    t = T[j - 1]
        
    p =  np.dot(t, theta)
    return p