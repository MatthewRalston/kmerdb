import logging
logger = logging.getLogger(__file__)

import math
import numpy as np
from numba import jit
import functools


from kdb import kmer, fileutil

identity = {
    'correlation': '1.0',
    'hamming'    : '1.0',
    'euclidean'  : '0.0'
}


def correlation(fname1, fname2):
    k = None
    with fileutil.open(fname1, mode='r') as kdb1:
        with fileutil.open(fname2, mode='r') as kdb2:
            if k is None:
                k = kdb1.header['k']
            if k != kdb1.header['k']:
                raise Exception("File '{0}' reported k = {1} instead of k = {2}".format(f, kdb1.header['k'], k))
            elif k != kdb2.header['k']:
                raise Exception("File '{0}' reported k = {1} instead of k = {2}".format(f, kdb2.header['k'], k))
            N = 4 ** k
            x_bar = functools.reduce(lambda a,b: a+b, map(lambda x: x['total_kmers'], kdb1.header['files']), 0) / N
            y_bar = functools.reduce(lambda a,b: a+b, map(lambda y: y['total_kmers'], kdb2.header['files']), 0) / N
            ## CALCULATE CORRELATION
            ssxx = 0
            ssyy = 0
            ssxy = 0
            for kmer_id in range(N):
                line1 = next(kdb1)
                line2 = next(kdb2)
                _, x = (int(_x) for _x in line1.rstrip().split("\t"))
                _, y = (int(_y) for _y in line2.rstrip().split("\t"))
                ssxx += np.square(x - x_bar)
                ssyy += np.square(y - y_bar)
                ssxy += (x - x_bar)*(y - y_bar)

            logger.debug("Sum of squared xy errors: {0}".format(ssxy))
            logger.debug("Sum of squared xx errors: {0}".format(ssxx))
            logger.debug("Sum of squared yy errors: {0}".format(ssyy))
            if ssxx*ssyy == 0: # Handle an irrational number
                return 0
            else:
                return ssxy/np.sqrt(ssxx*ssyy)



def hamming(k, x, y):
    sum = 0
    for i in range(len(x)):
        if x[i] == y[i]:
            sum += 1
    return (1/4**k) * sum



def d2s(mono_x, mono_y, total_kmers_x, total_kmers_y, k, x, y):

    mono_x = {c: np.round(mono_x[c]/total_kmers_x, 2) for c in mono_x}
    mono_y = {c: np.round(mono_y[c]/total_kmers_y, 2) for c in mono_y}
    
    N = len(x)

    def _d2s(ex, ey, xi, yi):
        xi_ = xi - (N-k+1)*ex
        yi_ = yi - (N-k+1)*ey

        return (xi_ * yi_)/np.sqrt(np.square(xi_) + np.square(yi_))
    
    s = 0
    for i in range(N):
        seq = kmer.kmer_to_id
        Ex = functools.reduce(lambda a,b: a*mono_x[b], list(seq), 1)
        Ey = functools.reduce(lambda a,b: a*mono_y[b], list(seq), 1)
        s += _d2s(Ex, Ey, x[i], y[i])

