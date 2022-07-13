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
import multiprocessing as mp
cimport numpy as cnp
cimport cython
import sys
#from numba import jit
#import functools


def pearson_correlation(a, b, total_kmers, sharedr):
    r = correlation(a, b, total_kmers)
    sharedr.value = r

#cpdef double correlation(cnp.uint64_t[:] a, cnp.uint64_t[:] b, int total_kmers):
cpdef double correlation(long[:] a, long[:] b, int total_kmers):

    cdef int i
    cdef long double ssxx = 0
    cdef long double ssyy = 0
    cdef long double ssxy = 0
    cdef long double xx = 0
    cdef long double yy = 0
    cdef long double xy = 0

    cdef float x_bar = np.sum(a)/total_kmers
    cdef float y_bar = np.sum(b)/total_kmers
    if total_kmers != len(a) or total_kmers != len(b):
        raise ValueError("NumPy kmer count array total does not match length of arrays")
    for i in range(total_kmers):
        try:
            xx = (a[i] - x_bar)**2
            yy = (b[i] - y_bar)**2
            xy = (a[i] - x_bar)*(b[i] - y_bar)
        except OverflowError as e:
            logger.error("{0}\t{1}   |   {2}\t{3}\t{4}".format(a[i], b[i], xx, yy, xy))
            logger.error(e)
            raise e
        

        ssxx += xx
        ssyy += yy
        if xy > 0: # Add a log of the negative of the negative number instead
            logger.debug("{0}\t{1}   |   {2}\t{3}\t{4} | ssxy = {5}".format(a[i], b[i], xx, yy, xy, ssxy))
            ssxy += xy
        else: # For large positive xy (not supposed to happen), we will interpret this as a inf/ssxx*ssyy

            logger.debug("{0}\t{1}   |   {2}\t{3}\t{4} | ssxy = {5}".format(a[i], b[i], xx, yy, xy, ssxy))
            ssxy += xy
            if ssxy < 0:
                logger.info("Sum of squared residuals is less than 0")
            elif ssxy == 0:
                logger.info("Sum of squared residuals is 0")
            elif ssxy > 0:
                logger.info("Calculating next ssxy...")
                continue
    logger.info("Custom Pearson correlation acquired")
    logger.info("{0}/sqrt({1}*{2})".format(ssxy, ssxx, ssyy))
    r = ssxy/(np.sqrt(ssxx*ssyy))
    return r
    
