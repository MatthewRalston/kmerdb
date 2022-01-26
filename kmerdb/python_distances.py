'''
   Copyright 2020 Matthew Ralston

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

'''




import logging
logger = logging.getLogger(__file__)

import math
import numpy as np
import sys
#from numba import jit
import functools


from kmerdb import kmer, fileutil

identity = {
    'correlation': '1.0',
    'euclidean'  : '0.0',
    'hamming'    : '1.0',
    'pearson'    : '1.0',
    'spearman'   : '1.0'
}


# Remove fileutil.KDBReader references
# Add numpy structure
#
#
#
#
#
def correlation(counts1:np.ndarray, counts2:np.ndarray):
    """
    A custom Pearson correlation function. Returns a float:

    ssxy over the square-root of ssxx*ssyy

    :param counts1: The counts array
    :type counts1: numpy.ndarray
    :param counts2: The second counts array
    :type counts2: numpy.ndarray
    :returns: float
    :rtype: float
    
    """
    if type(counts1) is not np.ndarray:
        raise TypeError("kmerdb.distance.correlation expects a NumPy array as its first positional argument")
    elif type(counts2) is not np.ndarray:
        raise TypeError("kmerdb.distance.correlation expects a NumPy array as its second positional argument")
    elif counts1.shape != counts2.shape:
        raise TypeError("kmerdb.distance.correlation expects the NumPy arrays to have the same dimension.")
    elif counts1.size != counts2.size:
        raise TypeError("kmerdb.distance.correlation expects two equal size NumPy arrays/vectors as its first and second positional arguments.")
    x_bar = np.sum(counts1)/counts1.size
    y_bar = np.sum(counts2)/counts2.size

    ssxx = 0
    ssyy = 0
    ssxy = 0
    for i, c in enumerate(counts1):
        x = c
        y = counts2[i]
        ssxx += np.square(x - x_bar)
        ssyy += np.square(y - y_bar)
        x_minus_xbar = x - x_bar
        y_minus_ybar = y - y_bar
        ssxy += x_minus_xbar*y_minus_ybar
        logger.warning("Beware custom correlations")
    if ssxx*ssyy == 0:
        logger.error("Incorrect denominator found, skipping...")
        return 0
    else:
        logger.info("Acquired")
        #logger.debug("Well done, cap...")
        return ssxy/np.sqrt(ssxx*ssyy)



def spearman(x, y):
    """
    Thin wrapper for scipy.stats.spearmanr

    :param x:
    :type x: numpy.ndarray
    :param y:
    :type y: numpy.ndarray
    :returns: (cor, pval)
    :rtype: tuple
    """
    if type(x) is not np.ndarray:
        raise TypeError("kmerdb.distance.spearman expects a Numpy array as its first positional argument")
    elif type(y) is not np.ndarray:
        raise TypeError("kmerdb.distance.spearman expects a Numpy array as its second positional argument")
    from scipy.stats import spearmanr
    cor, pval = spearmanr(x, b=y)
    logger.debug("Spearman calculation from SciPy")
    logger.info("Smooth, buttery Spearman correlation coefficients.")
    return cor, pval

def EMD(x, y):
    """
    Incomplete Earth Mover's Distance
    """
    if type(x) is not np.ndarray:
        raise TypeError("kmerdb.distance.EMD expects a Numpy array as its first positional argument")
    elif type(y) is not np.ndarray:
        raise TypeError("kmerdb.distance.EMD expects a Numpy array as its second positional argument")
    from scipy.stats import wasserstein_distance
    return wasserstein_distance(x, y)

        
def hamming(k, x, y):
    """
    Very old deprecated distance on an older data structure
    """
    sum = 0
    for i in range(len(x)):
        if x[i] == y[i]:
            sum += 1
    return (1/4**k) * sum



# def d2s(x, y):
#     if type(x) is not np.ndarray:
#         raise TypeError("kmerdb.distance.d2s expects a Numpy array as its first positional argument")
#     elif type(y) is not np.ndarray:
#         raise TypeError("kmerdb.distance.d2s expects a Numpy array as its second positional argument")

    
#     from kmerdb import kmer
#     import math
    
#     N = len(x)
#     k = int(math.log(N, 4))
#     total_kmers_x = np.sum(x)
#     total_kmers_y = np.sum(y)
#     #mono_x = dict([c, np.round(mono_x[c]/total_kmers_x, 2) for c in mono_x])
#     #mono_y = dict([c, np.round(mono_y[c]/total_kmers_y, 2) for c in mono_y])
#     mono_x = dict([c, mono_x[c]/float(total_kmers_x) for c in mono_x])
#     mono_y = dict([c, mono_y[c]/float(total_kmers_y) for c in mono_y])


#     def _d2s(ex, ey, xi, yi):
#         xi_ = xi - (N-k+1)*ex
#         yi_ = yi - (N-k+1)*ey

#         return (xi_ * yi_)/np.sqrt(np.square(xi_) + np.square(yi_))
    
#     s = 0
#     for i in range(N):
#         seq = kmer.kmer_to_id(i)
#         Ex = functools.reduce(lambda a,b: a*mono_x[b], list(seq), 1)
#         Ey = functools.reduce(lambda a,b: a*mono_y[b], list(seq), 1)
#         s += _d2s(Ex, Ey, x[i], y[i])

#     return s
