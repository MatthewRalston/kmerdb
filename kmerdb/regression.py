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

import sys



import logging
logger = logging.getLogger(__file__)



import numpy as np
#cimport numpy as cnp
#cimport cython


#from kmerdb cimport cublass_gemm

#from cpython.ref cimport PyObject


# def linear_regression(a, b):
#     x = regression(a, b, method="numpy")

# @cython.boundscheck(False)
# @cython.wraparound(False)
# def regression(double[:, :] A, double[:] b, char[:] method):
#     """
#     Linear regression wth two inputs.
#     A 2D array 'A' with two dimensions, and non-necessarily square.

#     and a vector of data determining the linear constants in the regression.
    

#     least squares problems occur when b is not in the column space of A.
    
#     They are solvable if and only if the matrix A:
#       * is non-singular
#       * has independent rows/columns
#       * a non-zero determinant (i.e. it is *not* row identical to the identity matrix)
#       * and is full-rank. 

    
#     ############
#     # definition
#     ############
#     At(b-Ax) = 0 OR
#     AtAx = Atb

#     but, less developed...
#     Ax = b

#     These are the normal equations.
    
#     """
#     cdef int n, m
#     n = A.shape[0]  # number of samples
#     m = A.shape[1]  # number of constituents


#     if n != m:
#         raise ValueError("kmerdb.regression expects a square matrix A of independent variables for regression against the dependent variable.")

    
#     cdef double[:, :] At = np.empty((m, n), dtype=np.double)
#     cdef double[:, :] At_test = np.empty((m, n), dtype=np.double)

#     # AtA: square symmetric positive definite required next if A is non-singular.
#     cdef double[:, :] AtA = np.empty((m, m), dtype=np.double)

    
#     cdef double[:] Atb = np.empty(m, dtype=np.double)
#     cdef double[:] x = np.empty(m, dtype=np.double)

#     # Calculate A^T
#     #At = A.T
#     for i in range(n):
#         for j in range(m):
#             At[i, j] = A[j][i]
#     At_test = A.T
#     assert np.array_equal(At_test, At), "kmerdb.regression expects the transpose of A to be equivalent to NumPy's transpose function."

#     # Calculate A^T * A

#     ## [MULTIPLY]
#     if method == "numpy": # Delegated to BLAS if configured.
#         AtA = np.dot(At, A)
#     elif method == "strassen":
#         raise ValueError("NOTE: not yet implemented. testing a cublas delegated matrix multiply.")

#         #AtA = strassen_multiplication(At, A)
#     elif method == "gemm":

#         #cdef cnp.ndarray[int64_t, ndim=2] a = A
#         #cdef cnp.ndarray[int64_t, ndim=2] b = B
#         #cdef cnp.ndarray[int64_t, ndim=2] c = C



#         raise ValueError("NOTE: not implemented. issues during Cythonize. ")

    
        
#         #cublas_gemm.matrix_mult_cuda_int64(a, b, c) # Uses cublas_gemm, specifically Dgemm

        
#         #AtA = c
#         #AtA = np.ndarray(c, dtype="uint64")

#         raise ValueError("NOTE: not yet implemented. testing a cublas delegated matrix multiply.")
        
#     # Ensure the solution can be found by checking for non-singularity.

#     detAtA = np.linalg.det(AtA)

#     if detAtA == 0:
#         raise ValueError("kmerdb.regression encountered a non-invertible matrix 'AtA', which is used in the solution of the least-squares problem.")
#     else:
#         try:
#             # TODO: Need an inverse function 
#             AtA_inv = np.linalg.inv(AtA)
#         except Exception as e:
#             raise e


#     # Calculate A^T * b
#     Atb = np.dot(At, b)
    

#     # Inverse A^T A and left-multiply with Atb. This yields the estimates for x that minimize error from the projection of b, p on the left nullspace of AtA. At => multiply with b. x is the Euclidean/Newtonian minimum of the convex space spanned by the left 
#     x = AtA_inv.dot(Atb)

#     return x


