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
#from libc.stdlib cimport malloc, free
#from libc.string cimport memcpy

# A simple utility function to add two matrices
def add(A, B):
    i: cython.int
    j: cython.int
    rows: cython.int
    cols: cython.int

    rows = A.shape[0]
    cols = A.shape[1]
    
    C = np.empty((rows, cols), dtype=np.double)
    for i in range(rows):
        for j in range(cols):
            C[i, j] = A[i, j] + B[i, j]

    return C

# cdef cnp.ndarray[double, ndim=2] add(cnp.ndarray[double, ndim=2] A, 
#                                      cnp.ndarray[double, ndim=2] B):
#     cdef int i, j
#     cdef int rows = A.shape[0]
#     cdef int cols = A.shape[1]
#     cdef cnp.ndarray[double, ndim=2] C = np.empty((rows, cols), dtype=np.double)

#     for i in range(rows):
#         for j in range(cols):
#             C[i, j] = A[i, j] + B[i, j]
    
#     return C

# A simple utility function to subtract two matrices
# cdef cnp.ndarray[double, ndim=2] subtract(cnp.ndarray[double, ndim=2] A, 
#                                           cnp.ndarray[double, ndim=2] B):
#     cdef int i, j
#     cdef int rows = A.shape[0]
#     cdef int cols = A.shape[1]
#     cdef cnp.ndarray[double, ndim=2] C = np.empty((rows, cols), dtype=np.double)

#     for i in range(rows):
#         for j in range(cols):
#             C[i, j] = A[i, j] - B[i, j]
    
#     return C

# # Strassen's algorithm for matrix multiplication
# cdef cnp.ndarray[double, ndim=2] strassen(cnp.ndarray[double, ndim=2] A,
#                                           cnp.ndarray[double, ndim=2] B, str method):
#     cdef int n

    
#     n = A.shape[0]
#     if method is None:
#         method = "numpy"
    
    
#     if method == "numpy":
#         return np.dot(A, B)
#     if n <= 4:#: and base == "numpy"  # Base case
#         return np.dot(A, B)
#     elif n > 4 and method == "strassen":
    
#         # Splitting the matrices into quadrants
#         half = n // 2
#     cdef cnp.ndarray[double, ndim=2] A11 = A[:half, :half]
#     cdef cnp.ndarray[double, ndim=2] A12 = A[:half, half:]
#     cdef cnp.ndarray[double, ndim=2] A21 = A[half:, :half]
#     cdef cnp.ndarray[double, ndim=2] A22 = A[half:, half:]
    
#     cdef cnp.ndarray[double, ndim=2] B11 = B[:half, :half]
#     cdef cnp.ndarray[double, ndim=2] B12 = B[:half, half:]
#     cdef cnp.ndarray[double, ndim=2] B21 = B[half:, :half]
#     cdef cnp.ndarray[double, ndim=2] B22 = B[half:, half:]
    
#     cdef cnp.ndarray[double, ndim=2] P1
#     cdef cnp.ndarray[double, ndim=2] P2
#     cdef cnp.ndarray[double, ndim=2] P3
#     cdef cnp.ndarray[double, ndim=2] P4
#     cdef cnp.ndarray[double, ndim=2] P5
#     cdef cnp.ndarray[double, ndim=2] P6
#     cdef cnp.ndarray[double, ndim=2] P7
    
#     cdef cnp.ndarray[double, ndim=2] temp1 = B12
#     cdef cnp.ndarray[double, ndim=2] temp2 = A11
#     cdef cnp.ndarray[double, ndim=2] temp3 = A21
#     cdef cnp.ndarray[double, ndim=2] temp4 = B11
    
#     cdef cnp.ndarray[double, ndim=2] temp5 = A11
#     cdef cnp.ndarray[double, ndim=2] temp6 = B11
#     cdef cnp.ndarray[double, ndim=2] temp7 = A22
#     cdef cnp.ndarray[double, ndim=2] temp8 = B22
    
#     cdef cnp.ndarray[double, ndim=2] temp9 = A11
#     cdef cnp.ndarray[double, ndim=2] temp10 = B11
    
#     # Compute the 7 products
    
#     temp1 = subtract(B12, B22)
#     P1 = strassen(A11, temp1, method="numpy")
    
#     temp2 = add(A11, A12)
#     P2 = strassen(temp2, B22, method="numpy")
    
#     temp3 = add(A21, A22)
#     P3 = strassen(temp3, B11, method="numpy")
    
#     temp4 = subtract(B21, B11)
#     P4 = strassen(A22, temp4, method="numpy")
    
#     temp5 = add(A11, A22)
#     temp6 = add(B11, B22)
#     P5 = strassen(temp5, temp6, method="numpy")
    
#     temp7 = subtract(A12, A22)
#     temp8 = add(B21, B22)
#     P6 = strassen(temp7, temp8, method="numpy")
    
    
#     temp9 = subtract(A11, A21)
#     temp10 = add(B11, B12)
#     P7 = strassen(temp9, temp10, method="numpy")
    
#     # Combine the products into a single result matrix
#     C = np.empty((n, n), dtype=np.double)
    
#     # C1 =  ((P5 + P4) - P2) + P6) (2 additions, 1 subtraction times exactly proportional to the matrix sizes.

#     # temp1_2 = (P6 + (P2 + (P4 + P5))) #
#     # temp2 = (P6 + (P5+P4))
#     temp1_1 = np.add(P5, P4) # P5(A12 + A22  x B11 + B22)  ADD+  P4(A22 x B21 - B11)
#     temp1_2 = np.subtract(temp1, P2)
    
#     C1 = np.add(temp1_2, P6)
    

    
#     # t7 = (A12 - A22) # 9/20/24
#     # t31 = (B21 + B22) # 9/20
#     # Parts of p5. guess which operator combines
    
#     # C2 = (A11 + A22)  x (B11 + B22)

#     # temp2_1 = add(temp5, temp6) # 9/20
#     # should be multipliy, which makes these useless
#     # And it's just p5 then
#     temp2_3 = np.add(P5, P6) #  ((A12 - A22) x (B11 + B22)) + P6 (A12 - A22) x (B21 + B22)
    
    
    
    
#     C2 = np.add(P1, P2)
#     C3 = np.add(P3, P4)
    

        
        
        
        
#     C[:half, :half] = add(subtract(P1, P2), P6)
#     C[:half, half:] = add(P1, P2)
#     C[half:, :half] = add(P3, P4)
#     C[half:, half:] = subtract(add(subtract(add(P5, P1), P3), P7), P6)
        
#     return C

#     # End Strassen



