
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
# cublas_gemm.pyx
import numpy as np
cimport numpy as cnp
from libc.stdlib cimport malloc, free
from libc.stdlib cimport memcpy
from cuda import cudaMalloc, cudaMemcpy, cudaFree
from cpython cimport PyObject

# Declare C interface to cuBLAS
cdef extern from "cublas_v2.h":
    ctypedef void* cublasHandle_t
    
    int cublasCreate(cublasHandle_t* handle) nogil
    int cublasDestroy(cublasHandle_t handle) nogil
    
    int cublasDgemm(cublasHandle_t handle, char transa, char transb,
                    int m, int n, int k,
                    double alpha, double *A, int lda,
                    double *B, int ldb,
                    double beta, double *C, int ldc)

# Function to perform matrix multiplication using cuBLAS
def matrix_mult_cuda_int64(cnp.ndarray[cnp.int64_t, ndim=2] A, 
                            cnp.ndarray[cnp.int64_t, ndim=2] B, 
                            cnp.ndarray[cnp.int64_t, ndim=2] C):
    cdef int m, n, k
    cdef int lda, ldb, ldc
    cdef double alpha = 1.0
    cdef double beta = 0.0

    # Set matrix dimensions
    m = A.shape[0]
    k = A.shape[1]
    n = B.shape[1]

    lda = m
    ldb = k
    ldc = m

    # Create the cuBLAS context
    cdef cublasHandle_t handle = NULL

    if cublasCreate(&handle) != 0:
        raise RuntimeError("Failed to create cuBLAS handle")
    
    #cublasCreate(&handle)  # Properly create the cuBLAS handle

    # Allocate device memory for inputs and result
    cdef double *d_A
    cdef double *d_B
    cdef double *d_C

    # Allocate device memory
    cudaMalloc(&d_A, m * k * sizeof(double))
    cudaMalloc(&d_B, k * n * sizeof(double))
    cudaMalloc(&d_C, m * n * sizeof(double))

    # Create temporary host arrays for conversion from int64 to double
    cdef double *h_A = <double *> malloc(m * k * sizeof(double))
    cdef double *h_B = <double *> malloc(k * n * sizeof(double))

    # Copy data from numpy arrays (int64) to temporary host arrays (double)
    for i in range(m * k):
        h_A[i] = <double> A[i // A.shape[1], i % A.shape[1]]  # Convert to double

    for i in range(k * n):
        h_B[i] = <double> B[i // B.shape[1], i % B.shape[1]]  # Convert to double

    # Copy from host arrays to device
    cudaMemcpy(d_A, h_A, m * k * sizeof(double), cudaMemcpyHostToDevice)
    cudaMemcpy(d_B, h_B, k * n * sizeof(double), cudaMemcpyHostToDevice)

    # Perform the matrix multiplication using cuBLAS
    cublasDgemm(handle, 'N', 'N', m, n, k, alpha, d_A, lda, 
                d_B, ldb, beta, d_C, ldc)

    # Copy result back from device to host (output C)
    #cudaMemcpy(C.data, d_C, m * n * sizeof(double), cudaMemcpyDeviceToHost)

    # Clean up device memory and host memory
    free(h_A)
    free(h_B)
    #cudaFree(d_A)
    #cudaFree(d_B)
    #cudaFree(d_C)
    cublasDestroy(handle)  # Properly destroy the cuBLAS handle
