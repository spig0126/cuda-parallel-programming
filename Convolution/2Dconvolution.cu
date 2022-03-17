/*
 * =====================================================================================
 *
 *       Filename:  main.cu
 *
 *    Description: 	Convolution
 *
 *        Version:  1.0
 *        Created:  02/21/2022 10:07:38
 *       Revision:  none
 *       Compiler:  nvcc
 *
 *         Author:  Yeonhee Jung
 *
 * =====================================================================================
 */

#include <assert.h>
#include "mkCuda.h"
#include "mkClockMeasure.h"

const int N_W = 256;
const int N_H = 256;
const int M_W = 5;
const int M_H = 5;
const int M_INTERVAL = M_W / 2;
const int MAX_NUM = 10;
const int MAX_ITER = 100;

unsigned int N[N_W * N_H];
unsigned int M[M_W * M_H];

const int tbSize = N_W * N_H;
const int threads = 16;
const int TILE_SIZE = threads;

unsigned int cpuOut[N_W * N_H];
unsigned int gpuOut_2D_basic[N_W * N_H];
unsigned int gpuOUt_2D_CM[N_W * N_H];
unsigned int gpuOUt_2D_tiled[N_W * N_H];
unsigned int gpuOUt_2D_tiled_internal[N_W * N_H];
__constant__ unsigned int d_m_CM[M_W * M_H];

void generateRandomValuesFor2D(unsigned int *input, const int row, const int col)
{
    for (int r = 0; r < row; r++)
    {
        for (int c = 0; c < col; c++)
        {
            input[r * col + c] = (unsigned int)float(rand()) / float(RAND_MAX) * MAX_NUM;
        }
    }
}

void print2DMask(const unsigned int *input, const int row, const int col)
{
    printf("Print 2D Mask \n -----------\n");
    for (int r = 0; r < row; r++)
    {
        for (int c = 0; c < col; c++)
        {
            printf("%u\t", input[r * col + c]);
        }
        printf("\n");
    }
    printf("--------\n");
}

void print2DMatrix(const unsigned int *input, const int row, const int col)
{
    printf("Print 2D Matrix \n -------------------\n");
    for (int r = 0; r < row; r++)
    {
        for (int c = 0; c < col; c++)
        {
            printf("%u\t", input[r * col + c]);
        }
        printf("\n");
    }
    printf("-----------------\n");
}

bool compare2DMatrix(const unsigned int *inputA, const unsigned int *inputB, const int row, const int col)
{
    bool ret = true;
    for (int i = 0; i < row * col; i++)
    {
        if (inputA[i] != inputB[i])
        {
            ret = false;
            break;
        }
    }
    return ret;
}

void cpu_2DConvolution_basic(const unsigned int *n, unsigned int *m, unsigned int *p, const int row, const int col, const int m_row, const int m_col, const int interval)
{
    for (int r = 0; r < row; r++)
    {
        for (int c = 0; c < col; c++)
        {
            int pVal = 0;
            int nRowStartIndex = r - interval;
            int nColStartIndex = c - interval;
            for (int r_m = 0; r_m < m_row; r_m++)
            {
                for (int c_m = 0; c_m < m_col; c_m++)
                {
                    int currentRow = nRowStartIndex + r_m;
                    int currentCol = nColStartIndex + c_m;
                    if(currentRow >= 0 && currentRow <row && currentCol>=0 && currentCol <col ){
                        int nIndex = currentRow * col + currentCol;
                        pVal += n[nIndex] * m[r_m * m_col + c_m];
                    }
                }
            }
            p[r * col + c] = pVal;
        }
    }
}

__global__ void gpu_2DConvolution_basic(const unsigned int *n, unsigned int *m, unsigned int *p, const int width, const int height, const int m_width, const int m_height, const int interval)
{
    int row = blockDim.y * blockIdx.y + threadIdx.y;
    int col = blockDim.x * blockIdx.x + threadIdx.x;
    int pVal = 0;
    int row_start_point = row - interval;
    int col_start_point = col - interval;

    for (int r_m = 0; r_m < m_height; r_m++)
    {
        for (int c_m = 0; c_m < m_width; c_m++)
        {
            int currentRow = row_start_point + r_m;
            int currentCol = col_start_point + c_m;
            if (currentRow >= 0 && currentRow < height && currentCol >= 0 && currentCol < width)
            {
                int nIndex = currentRow * width + currentCol;
                pVal += n[nIndex] * m[r_m * m_width + c_m];
            }
        }
    }

    p[row * width + col] = pVal;
}

__global__ void gpu_2DConvolution_cm(const unsigned int *n, unsigned int *p, const int width, const int height, const int m_width, const int m_height, const int interval)
{
    int row = blockDim.y * blockIdx.y + threadIdx.y;
    int col = blockDim.x * blockIdx.x + threadIdx.x;
    int pVal = 0;
    int row_start_point = row - interval;
    int col_start_point = col - interval;

    for (int r_m = 0; r_m < m_height; r_m++)
    {
        for (int c_m = 0; c_m < m_width; c_m++)
        {
            int currentRow = row_start_point + r_m;
            int currentCol = col_start_point + c_m;
            if (currentRow >= 0 && currentRow < height && currentCol >= 0 && currentCol < width)
            {
                int nIndex = currentRow * width + currentCol;
                pVal += n[nIndex] * d_m_CM[r_m * m_width + c_m];
            }
        }
    }

    p[row * width + col] = pVal;
}

__global__ void gpu_2DConvolution_tiled(const unsigned int *n, unsigned int *p, const int width, const int height, const int m_width, const int m_height, const int interval, const int TILE_SIZE)
{
    __shared__ unsigned int N_cm[TILE_SIZE + m_width - 1][TILE_SIZE + m_height - 1];
    int row = blockDim.y * blockIdx.y + threadIdx.y;
    int col = blockDim.x * blockIdx.x + threadIdx.x;
    int tx = threadIdx.x;
    int ty = threadIdx.y;
    int row_start_point = row - interval;
    int col_start_point = col - interval;

    if(row_start_point >= 0 && row_start_point < height && col_start_point >= 0 && col_start_point<width){
        N_cm[ty][tx] = n[row_start_point * width + col_start_point];
    }
    else{
        N_cm[ty][tx] = 0;
    }

    int pVal = 0;

    if(ty < )

    for (int r_m = 0; r_m < m_height; r_m++)
    {
        for (int c_m = 0; c_m < m_width; c_m++)
        {
            int currentRow = row_start_point + r_m;
            int currentCol = col_start_point + c_m;
            if (currentRow >= 0 && currentRow < height && currentCol >= 0 && currentCol < width)
            {
                int nIndex = currentRow * width + currentCol;
                pVal += n[nIndex] * d_m_CM[r_m * m_width + c_m];
            }
        }
    }

    p[row * width + col] = pVal;
}
int main()
{
    srand((unsigned int)time(NULL));
    generateRandomValuesFor2D(N, N_W, N_H);
    generateRandomValuesFor2D(M, M_W, M_H);

    // MK: GPU Memory
    unsigned int *d_n, *d_m, *d_p;

    int n_byteSize = N_W * N_H * sizeof(unsigned int);
    int m_byteSize = M_W * M_H * sizeof(unsigned int);

    // allocate memory in device
    cudaError_t err = cudaMalloc((void **)&d_n, n_byteSize);
    checkCudaError(err);
    err = cudaMalloc((void **)&d_m, m_byteSize);
    checkCudaError(err);
    err = cudaMalloc((void **)&d_p, n_byteSize);
    checkCudaError(err);

    err = cudaMemcpy(d_n, N, n_byteSize, cudaMemcpyHostToDevice);
    checkCudaError(err);
    err = cudaMemcpy(d_m, M, m_byteSize, cudaMemcpyHostToDevice);
    checkCudaError(err);
    err = cudaMemcpyToSymbol(d_m_CM, M, m_byteSize);
    checkCudaError(err);

    // MK: Thread Num
    // const int tbSize = N_W * N_H;
    const int threads = 16;
    dim3 gridSize(ceil((float)N_W / threads), ceil((float)N_H / threads), 1);
    dim3 blockSize(threads, threads, 1);

    // MK: Time Measurement
    mkClockMeasure *ckCpu = new mkClockMeasure("CPU CODE");
    ckCpu->clockReset();

    mkClockMeasure *ckGpu = new mkClockMeasure("GPU CODE");
    ckGpu->clockReset();
    mkClockMeasure *ckGpu_CM = new mkClockMeasure("GPU CODE - CM");
    ckGpu_CM->clockReset();

    for (int i = 0; i < MAX_ITER; i++)
    {
        ckCpu->clockResume();
        cpu_2DConvolution_basic(N, M, cpuOut, N_W, N_H, M_W, M_H, M_INTERVAL);
        ckCpu->clockPause();

        ckGpu->clockResume();
        gpu_2DConvolution_basic<<<gridSize, blockSize>>>(d_n, d_m, d_p, N_W, N_H, M_W, M_H, M_INTERVAL);
        err = cudaDeviceSynchronize();
        ckGpu->clockPause();

        ckGpu_CM->clockResume();
        gpu_2DConvolution_cm<<<gridSize, blockSize>>>(d_n, d_p, N_W, N_H, M_W, M_H, M_INTERVAL);
        err = cudaDeviceSynchronize();
        ckGpu_CM->clockPause();
        checkCudaError(err);
    }

    err = cudaMemcpy(gpuOut_2D_basic, d_p, n_byteSize, cudaMemcpyDeviceToHost);
    checkCudaError(err);

    cudaFree(d_n);
    cudaFree(d_m);
    cudaFree(d_p);

    if (compare2DMatrix(cpuOut, gpuOut_2D_basic, N_W, N_H))
    {
        printf("---------- [CPU] 2D convolution basic-------------\n");
        ckCpu->clockPrint();
        printf("---------- [GPU] 2D convolution basic-------------\n");
        ckGpu->clockPrint();
        printf("---------- [GPU] 2D convolution constant memory-------------\n");
        ckGpu_CM->clockPrint();
    }
    else
    {
        printf("ERROR: Two Matrices are not same\n");
    }
    // print2DMask(M, M_W, M_H);
    // print2DMatrix(N, N_W, N_H);
    // print2DMatrix(cpuOut, N_W, N_H);
    // print2DMatrix(gpuOut_2D_basic, N_W, N_H);

    // print1DMatrixValue(cpuOut, N_WIDTH);
    // print1DMatrixValue(gpuOut_1D_basic, N_WIDTH);
}
