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



import numpy as np
cimport numpy as cnp
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy

# Ensure appropriate imports for types
from libc.stdint cimport int64_t

# A simple utility function to add two matrices
cdef cnp.ndarray[double, ndim=2] add(cnp.ndarray[double, ndim=2] A, 
                                     cnp.ndarray[double, ndim=2] B):
    cdef int i, j
    cdef int rows = A.shape[0]
    cdef int cols = A.shape[1]
    cdef cnp.ndarray[double, ndim=2] C = np.empty((rows, cols), dtype=cnp.double)

    for i in range(rows):
        for j in range(cols):
            C[i, j] = A[i, j] + B[i, j]
    
    return C

# A simple utility function to subtract two matrices
cdef cnp.ndarray[double, ndim=2] subtract(cnp.ndarray[double, ndim=2] A, 
                                          cnp.ndarray[double, ndim=2] B):
    cdef int i, j
    cdef int rows = A.shape[0]
    cdef int cols = A.shape[1]
    cdef cnp.ndarray[double, ndim=2] C = np.empty((rows, cols), dtype=cnp.double)

    for i in range(rows):
        for j in range(cols):
            C[i, j] = A[i, j] - B[i, j]
    
    return C

# Strassen's algorithm for matrix multiplication
cdef cnp.ndarray[double, ndim=2] strassen(cnp.ndarray[double, ndim=2] A,
                                          cnp.ndarray[double, ndim=2] B):
    cdef int n
    n = A.shape[0]

    if n <= 2:  # Base case
        return np.dot(A, B)

    # Splitting the matrices into quadrants
    cdef int half = n // 2
    cdef cnp.ndarray[double, ndim=2] A11 = A[:half, :half]
    cdef cnp.ndarray[double, ndim=2] A12 = A[:half, half:]
    cdef cnp.ndarray[double, ndim=2] A21 = A[half:, :half]
    cdef cnp.ndarray[double, ndim=2] A22 = A[half:, half:]

    cdef cnp.ndarray[double, ndim=2] B11 = B[:half, :half]
    cdef cnp.ndarray[double, ndim=2] B12 = B[:half, half:]
    cdef cnp.ndarray[double, ndim=2] B21 = B[half:, :half]
    cdef cnp.ndarray[double, ndim=2] B22 = B[half:, half:]

    # Compute the 7 products
    P1 = strassen(A11, subtract(B12, B22))
    P2 = strassen(add(A11, A12), B22)
    P3 = strassen(add(A21, A22), B11)
    P4 = strassen(A22, subtract(B21, B11))
    P5 = strassen(add(A11, A22), add(B11, B22))
    P6 = strassen(subtract(A12, A22), add(B21, B22))
    P7 = strassen(subtract(A11, A21), add(B11, B12))

    # Combine the products into a single result matrix
    C = np.empty((n, n), dtype=cnp.double)
    C[:half, :half] = add(subtract(add(P5, P4), P2), P6)
    C[:half, half:] = add(P1, P2)
    C[half:, :half] = add(P3, P4)
    C[half:, half:] = subtract(add(subtract(add(P5, P1), P3), P7), P6)

    return C

# Batch processing function for Strassen's algorithm
cpdef batch_strassen(cnp.ndarray[double, ndim=2] A_batch,
                     cnp.ndarray[double, ndim=2] B_batch):
    cdef int batch_size, i
    batch_size = A_batch.shape[0]

    # Result batch initialization
    result_shape = (batch_size, A_batch.shape[1], A_batch.shape[2])
    cdef cnp.ndarray[double, ndim=3] C_batch = cnp.empty(result_shape, dtype=cnp.double)

    for i in range(batch_size):
        C_batch[i] = strassen(A_batch[i], B_batch[i])

    return C_batch



