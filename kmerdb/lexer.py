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
# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False



import numpy as np

def max(X):
    length: cython.int
    m: cython.int
    l: cython.int

    if type(X) is list:
        length = len(X)
    elif type(X) is np.ndarray:
        length = X.shape[0]
    l = 0
    m = 0
    if (length, ) != (len(X), ):
        raise ValueError("kmerdb.lexer.max expects a unidimensional array")
    for i in range(length):
        #print("The ith value in the histo: {0} <=> the current max: {1}".format(X[i], m))

        
        if X[i] > m:
            l = i
            m = int(X[i])
    #print("The index i: {0} <=> the current max {1}".format(l, m))
    return l, m



