# -*- coding: utf-8 -*-

import math
import itertools
import numpy as np
from calculateChebyshevPolynomials import chebypol
    
def init(D,N,M):
    
    # D-dimensional nodes for interval [-1,1]    
    
    nu = np.empty(M)
    for k, nuk in enumerate(nu):
        nuk = - math.cos((2 * k - 1) *  math.pi / (2 * M))
    
    # create a matrix that contains all possible length D product of nodes in each row
    
    prod = itertools.product(nu, repeat = D)
    
    X = np.empty(pow(M, D) * D)
    X.shape = (pow(M, D), D)
    for i, row in enumerate(prod):
        X[i] = row
        
    # create a matrix B, which contains vector of polynomials in each row
        
    B = chebypol(N, nu)
    
    # create a matrix T, where each row represents all possible tensor products
    
    T = B
    for i in range(D - 1):
        T = np.kron(T, B)
            
    return [X, T]