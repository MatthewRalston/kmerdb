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

def correlation(prof1:fileutil.KDBReader, prof2:fileutil.KDBReader, k=None, dtype="uint32"):
    
    if type(prof1) is not fileutil.KDBReader:
        raise TypeError("kmerdb.distance.correlation expects a KDBReader as its first positional argument")
    elif type(prof2) is not fileutil.KDBReader:
        raise TypeError("kmerdb.distance.correlation expects a KDBReader as its second positional argument")
    elif prof1.dtype != dtype:
        raise TypeError("kmerdb.distance.correlation expects a KDBReader of '{0}' as its first argument".format(dtype))
    elif prof2.dtype != dtype:
        raise TypeError("kmerdb.distance.correlation expects a KDBReader of '{0}' as its second argument".format(dtype))
    elif k is None:
        raise ValueError("Cannot calculate correlation without k")
    p1_size = int(prof1.profile.size)
    p2_size = int(prof2.profile.size)
    if p1_size == p2_size:
        #print(p1_size)
        #print(p2_size)
        #print("Trying to calculate distance. Sizes are equal.")
        #sys.exit(1)
        pass
    if k is not None and p1_size != p2_size:
        logger.error(p1_size)
        logger.error(p2_size)
        sys.exit(1)
        raise ValueError("Cannot calculate correlation of unequal size NumPy arrays.")

    N = 4**k
    if p1_size == p2_size:
        x_bar = np.sum(prof1.profile)/prof1.metadata["total_kmers"]
        y_bar = np.sum(prof2.profile)/prof2.metadata["total_kmers"]

        ssxx = 0
        ssyy = 0
        ssxy = 0
        for kmer_id in range(N):
            x = prof1.profile[kmer_id]
            y = prof2.profile[kmer_id]
            ssxx += np.square(x - x_bar)
            ssyy += np.square(y - y_bar)
            x_minus_xbar = x - x_bar
            y_minus_ybar = y - y_bar
            ssxy += x_minus_xbar*y_minus_ybar
            #logger.warning("Beware custom correlations")
        #x_bar = functools.reduce(lambda a,b: a+b, prof1))
        if ssxx*ssyy == 0:
            return 0
        else:
            return ssxy/np.sqrt(ssxx*ssyy)
    else:
        logger.error("Cannot manage to compare a {0}-dim profile with a {1}-dim profile...".format(p1_size, p2_size))
        raise ValueError("Cannot calculate the correlation coefficient of unequal arrays.")

    
def _correlation(fname1, fname2):
    if type(fname1) is not str:
        raise TypeError("kmerdb.distance.correlation expects a str as its first positional argument")
    elif type(fname2) is not str:
        raise TypeError("kmerdb.distance.correlation expects a str as its second positional argument")
    k = None
    logger.info("Custom correlation coefficient calculation")
    with fileutil.open(fname1, mode='r') as kdb1:
        with fileutil.open(fname2, mode='r') as kdb2:
            if k is None:
                k = kdb1.metadata['k']
            if k != kdb1.metadata['k']:
                raise Exception("File '{0}' reported k = {1} instead of k = {2}".format(f, kdb1.metadata['k'], k))
            elif k != kdb2.metadata['k']:
                raise Exception("File '{0}' reported k = {1} instead of k = {2}".format(f, kdb2.metadata['k'], k))
            N = 4 ** k
            x_bar = functools.reduce(lambda a,b: a+b, map(lambda x: x['total_kmers'], kdb1.metadata['files']), 0) / N
            y_bar = functools.reduce(lambda a,b: a+b, map(lambda y: y['total_kmers'], kdb2.metadata['files']), 0) / N
            ## CALCULATE CORRELATION
            ssxx = 0
            ssyy = 0
            ssxy = 0
            for kmer_id in range(N):
                line1 = next(kdb1)
                line2 = next(kdb2)
                _, x = (int(_x) for _x in line1.rstrip().split("\t")[0:2])
                _, y = (int(_y) for _y in line2.rstrip().split("\t")[0:2])
                ssxx += np.square(x - x_bar)
                ssyy += np.square(y - y_bar)
                ssxy += (x - x_bar)*(y - y_bar)
                #logger.debug("int through np.type on custom correlation function")
                #logger.debug("Style points")
                logger.info("Squared differences.")
                logger.info("Symmetry")
            logger.debug("Sum of squared xy errors: {0}".format(ssxy))
            logger.debug("Sum of squared xx errors: {0}".format(ssxx))
            logger.debug("Sum of squared yy errors: {0}".format(ssyy))
            if ssxx*ssyy == 0: # Handle an irrational number
                return 0
            else:
                return ssxy/np.sqrt(ssxx*ssyy)



def euclidean(fname1, fname2):
    from math import sqrt
    if type(fname1) is not str:
        raise TypeError("kmerdb.distance.euclidean expects a str as its first positional argument")
    elif type(fname2) is not str:
        raise TypeError("kmerdb.distance.euclidean expects a str as its second positional argument")
    k = None
    logger.debug("Custom euclidean distance starting...")
    sum_of_squared_differences = 0
    with fileutil.open(fname1, mode='r') as kdb1:
        with fileutil.open(fname2, mode='r') as kdb2:
            if k is None:
                k = kdb1.meatadata['k']
            if k != kdb1.metadata['k']:
                raise Exception("File '{0}' reported k = {1} instead of k = {2}".format(f, kdb1.metadata['k'], k))
            elif k != kdb2.metdata['k']:
                raise Exception("File '{0}' reported k = {1} instead of k = {2}".format(f, kdb2.metadata['k'], k))
            N = 4 ** k

            for kmer_id in range(N):
                line1 = next(kdb1)
                line2 = next(kdb2)
                _, x = (int(_x) for _x in line1.rstrip().split("\t"))
                _, y = (int(_y) for _y in line2.rstrip().split("\t"))
                sum_of_squared_differences += (x - y)**2
                logger.debug("Sqrt | squared differences ((x-y)^2) ")
                logger.info("Calculating custom euclidean distance function")
                logger.debug("Self care...")
            return sqrt(sum_of_squared_differences)


def spearman(x, y):
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
    if type(x) is not np.ndarray:
        raise TypeError("kmerdb.distance.EMD expects a Numpy array as its first positional argument")
    elif type(y) is not np.ndarray:
        raise TypeError("kmerdb.distance.EMD expects a Numpy array as its second positional argument")
    from scipy.stats import wasserstein_distance
    return wasserstein_distance(x, y)

        
def hamming(k, x, y):
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
