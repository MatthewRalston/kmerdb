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
cimport numpy as cnp
cimport cython
import sys
#from numba import jit
#import functools

#cpdef double correlation(cnp.uint64_t[:] a, cnp.uint64_t[:] b, int total_kmers):
cpdef double correlation(long[:] a, long[:] b, int total_kmers):

    cdef int i
    cdef float ssxx = 0
    cdef float ssyy = 0
    cdef float ssxy = 0

    cdef float x_bar = np.sum(a)/len(a)
    cdef float y_bar = np.sum(b)/len(b)
    
    for i in range(total_kmers):
        ssxx += np.square(a[i] - x_bar)
        ssyy += np.square(b[i] - y_bar)
        ssxy += (a[i] - x_bar)*(b[i] - y_bar)
    if (ssxx*ssyy) == 0:
        logger.error("Incorrect denominator found, skipping")
        return 0
    else:
        logger.info("Custom Pearson correlation acquired")
        return ssxy/np.sqrt(ssxx*ssyy)
    
