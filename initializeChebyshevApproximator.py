############################
#
# Quantitative Marketing PS2
# Yu Ushioda
#
############################

# Chebyshev polynomial

import math
import itertools
import numpy as np

def Tn(N,x):
    T = cos(N*acos(x))
    return T
    
def init(D,N,M,a,b):
    
    # D-dimensional nodes for interval [-1,1]    
    
    nu = np.empty(M)
    for k, nuk in enumerate(nu):
        nuk = - math.cos((2 * k - 1) *  math.pi / (2 * M))
    
    # create a matrix that contains all possible length D product of nodes in each row
    
    perm = itertools.product(nu, repeat = D)
    
    X = np.empty(pow(M, D) * D)
    X.shape = (pow(M, D), D)
    for i, row in enumerate(perm):
        X[i] = row