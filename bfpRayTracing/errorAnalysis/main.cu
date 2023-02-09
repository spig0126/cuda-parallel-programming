#include "bfp.cuh"
#include "mkCuda.h"
#include <iostream>
#include <limits>
#include <fstream>
#include <string>

using namespace bfp;

const unsigned long GRID_SIZE = 256;        // 2^8
const unsigned long long BLOCK_SIZE = 1024; // 2^10
const unsigned long long N_SIZE = (unsigned long long)GRID_SIZE * (unsigned long long)BLOCK_SIZE * 4;

float h_n[N_SIZE];
unsigned long long cnt = 0;

float epsilon = std::numeric_limits<float>::epsilon();
typedef numeric_limits<float> flt;

__global__ void gpu_two_num_add(float *d_n, int a, int b_front)
{
    int b = (b_front << 18) + blockIdx.x * blockDim.x + threadIdx.x;
    unsigned long long idx = (blockIdx.x * blockDim.x + threadIdx.x) * 4;

    float *a_f = reinterpret_cast<float *>(&a);
    float *b_f = reinterpret_cast<float *>(&b);

    bfpNum a_b = float_to_bfpNum(*a_f);
    bfpNum b_b = float_to_bfpNum(*b_f);

    float res_f = *a_f + *b_f;
    float res_b = bfpNum_to_float(a_b + b_b);

    d_n[idx] = *a_f;
    d_n[idx + 1] = *b_f;
    d_n[idx + 2] = res_f;
    d_n[idx + 3] = res_b;
}

__global__ void gpu_two_num_mult(float *d_n, int a, int b_front)
{
    int b = (b_front << 18) + blockIdx.x * blockDim.x + threadIdx.x;
    unsigned long long idx = (blockIdx.x * blockDim.x + threadIdx.x) * 4;

    float *a_f = reinterpret_cast<float *>(&a);
    float *b_f = reinterpret_cast<float *>(&b);

    bfpNum a_b = float_to_bfpNum(*a_f);
    bfpNum b_b = float_to_bfpNum(*b_f);

    float res_f = *a_f * *b_f;
    float res_b = bfpNum_to_float(a_b * b_b);

    d_n[idx] = *a_f;
    d_n[idx + 1] = *b_f;
    d_n[idx + 2] = res_f;
    d_n[idx + 3] = res_b;
}

__global__ void gpu_two_num_div(float *d_n, int a, int b_front)
{
    int b = (b_front << 18) + blockIdx.x * blockDim.x + threadIdx.x;
    unsigned long long idx = (blockIdx.x * blockDim.x + threadIdx.x) * 4;

    float *a_f = reinterpret_cast<float *>(&a);
    float *b_f = reinterpret_cast<float *>(&b);

    bfpNum a_b = float_to_bfpNum(*a_f);
    bfpNum b_b = float_to_bfpNum(*b_f);

    float res_f = *a_f / *b_f;
    float res_b = bfpNum_to_float(a_b / b_b);

    d_n[idx] = *a_f;
    d_n[idx + 1] = *b_f;
    d_n[idx + 2] = res_f;
    d_n[idx + 3] = res_b;
}
__host__ void write_error_cases(float *h_n, ofstream &writeFile, string op)
{
    for (int i = 0; i < N_SIZE; i += 4)
    {
        if (isinf(h_n[i + 2]) || isnan(h_n[i + 2]))
            continue;
        if (!(fabs(h_n[i + 2] - h_n[i + 3]) < epsilon))
        {
            cnt++;
            if (cnt < 100)
            {
                writeFile << h_n[i] << op << h_n[i + 1] << " = " << h_n[i + 2] << " -> "
                          << h_n[i + 3] << endl;
            }
            else
                break;
        }
    }
}

__host__ void cnt_error(float *h_n, int a, int b_front, ofstream &cntFile, ofstream &caseFile, string op, bool writeCases)
{
    for (int i = 0; i < N_SIZE; i += 4)
    {
        if (isinf(h_n[i + 2]) || isnan(h_n[i + 2]))
            continue;
        if (!(fabs(h_n[i + 2] - h_n[i + 3]) < epsilon))
        {
            cnt++;
            if (writeCases && cnt < 100)
            {
                caseFile << h_n[i] << op << h_n[i + 1] << " = " << h_n[i + 2] << " -> "
                          << h_n[i + 3] << endl;
            }
        }
    }

    if(b_front == 0x3fff){
        if(cnt == 0){
            cntFile << a << endl;
        }
        else{
            cntFile << "a: " << a << ", errors: " << cnt << endl;
        }
    }
}


int main(int argc, char **argv)
{
    // int a = 33;
    int a = std::stoi(argv[1]);
    ofstream cntFile, caseFile;
    // cntFile.open("addErrorCnt.txt", std::ios_base::app);
    cntFile.open("multErrorCnt.txt", std::ios_base::app);
    // cntFile.open("divErrorCnt.txt", std::ios_base::app);

    // caseFile.open("addErrorCase.txt", std::ios_base::out);
    caseFile.open("multErrorCase.txt", std::ios_base::out);
    // caseFile.open("divErrorCase.txt", std::ios_base::out);
    caseFile.precision(flt::max_digits10 + 2);

    // GPU memory
    float *d_n;
    size_t N_BYTE_SIZE = N_SIZE * sizeof(float);

    // allocate memory in device
    cudaError_t err = cudaMalloc((void **)&d_n, N_BYTE_SIZE);
    checkCudaError(err);

    // grid setttings
    dim3 gridSize(GRID_SIZE, 1, 1);
    dim3 blockSize(BLOCK_SIZE, 1, 1);

    for (int b_front = 0; b_front <= 0x3fff; b_front++)
    {
        gpu_two_num_mult<<<gridSize, blockSize>>>(d_n, a, b_front);
        // gpu_two_num_add<<<gridSize, blockSize>>>(d_n, a, b_front);
        // gpu_two_num_div<<<gridSize, blockSize>>>(d_n, a, b_front);
        err = cudaMemcpy(h_n, d_n, N_BYTE_SIZE, cudaMemcpyDeviceToHost);
        checkCudaError(err);

        cnt_error(h_n, a, b_front, cntFile, caseFile, " * ", true);

    }

    // std::cout << "a: " << a << ", errors: " << cnt << endl;
    
    cudaFree(d_n);
    cntFile.close();
    caseFile.close();

    return 0;
}