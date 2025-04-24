#include <iostream>
#include <cuda_runtime.h>

__global__ void VecAdd(float* A, float* B, float* C) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if(i < N)
        C[i] = A[i] + B[i];
}

const int N = 10;

int main() {
    float h_A[N], h_B[N], h_C[N];
    float *d_A, *d_B, *d_C;
    size_t size = N * sizeof(float);

    // Initialize host arrays
    for (int i = 0; i < N; i++) {
        h_A[i] = 1.0f * i;
        h_B[i] = 10.0f + i;
    }

    // Allocate device memory
    cudaMalloc(&d_A, size);
    cudaMalloc(&d_B, size);
    cudaMalloc(&d_C, size);

    // Copy data from host to device
    cudaMemcpy(d_A, h_A, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_B, size, cudaMemcpyHostToDevice);

    // Launch kernel
    VecAdd<<<1, N>>>(d_A, d_B, d_C);

    // Copy result back to host
    cudaMemcpy(h_C, d_C, size, cudaMemcpyDeviceToHost);

    // Print results
    for (int i = 0; i < N; i++) {
        std::cout << h_A[i] << " + " << h_B[i] << " = " << h_C[i] << std::endl;
    }

    // Free device memory
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);

    return 0;
}
