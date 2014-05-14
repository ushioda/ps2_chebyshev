# -*- coding: utf-8 -*-

import math
import numpy as np
from calculateChebyshevCoefficients import coef
from calculateChebyshevPolynomials import chebypol

def approx(x,f,D,N,M,a,b):
    
    theta = coef(f,D,N,M,a,b)
    norm_x = - 1 + 2 * (x - a) / (b - a)
    B = chebypol(N, x)
    
    T = B
    for i in range(D - 1):
        T = np.kron(T, B)
        
    p =  np.dot(T, theta)
    return p