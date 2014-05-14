# -*- coding: utf-8 -*-

# Chebyshev polynomial

import math
import numpy as np

#########################################################
#
#    The following function returns 
#    N + 1 Chebyshev polynomials for some vector x 
#    (i, k) element of the output is 
#    the k-th polynomial of the i-th element of x
#
#########################################################

def chebypol(N, x):  
    T = np.ones(len(x))
    T1 = x
    T = np.vstack((T, T1))
    for i in range(N - 1):           
        t = 2 * x * T[-1] - T[-2]
        T = np.vstack((T, t))
    T = T.transpose()
    return T