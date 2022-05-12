/*
 * =====================================================================================
 *
 *       Filename:  main.cu
 *
 *    Description: 	2D-Convolution
 *
 *        Version:  1.0
 *        Created:  2022/02/21 10:07:38
 *       Revision:  none
 *       Compiler:  nvcc
 *
 *         Author:  Yeonhee Jung
 *   Organization:  Ewha Womans University
 *
 * =====================================================================================
 */

#include <assert.h>
#include "mkCuda.h"
#include "mkClockMeasure.h"

const int N_W = 512;
const int N_H = 512;
const int M_W = 5;
const int M_H = 5;
const int M_INTERVAL = M_W / 2;
const int MAX_NUM = 10;
const int MAX_ITER = 100;

unsigned int N[N_W * N_H];
unsigned int M[M_W * M_H];

const int threads = 32;
const int TILE_SIZE = threads;

unsigned int cpuOut[N_W * N_H];
unsigned int gpuOut_2D_basic[N_W * N_H];
unsigned int gpuOut_2D_CM[N_W * N_H];
unsigned int gpuOut_2D_tiled[N_W * N_H];
unsigned int gpuOut_2D_tiled_internal[N_W * N_H];
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

__global__ void gpu_2DConvolution_tiled(const unsigned int *n, unsigned int *p, const int width, const int height, const int m_width, const int m_height, const int interval, const int tile_size)
{
    __shared__ unsigned int SM[TILE_SIZE + M_W - 1][TILE_SIZE + M_H - 1];
 
    int tx = threadIdx.x;
    int ty = threadIdx.y;
    //output indexes
    int row_o = ty + blockIdx.y * tile_size;
    int col_o = tx + blockIdx.x * tile_size;
    //input indexes
    int row_i = row_o - interval;
    int col_i = col_o - interval;

    //load tile elements
    if(row_i >= 0 && row_i < height && col_i >= 0 && col_i<width){
        SM[ty][tx] = n[row_i * width + col_i];
    }
    else{
        SM[ty][tx] = 0;
    }

    //wait till all tile elements are loaded
    __syncthreads();

    //output value
    int pVal = 0;

    //only compute if it is an output tile element
    if(ty < tile_size && tx < tile_size){
        for (int m_row = 0; m_row < m_height; m_row++)
        {
            for (int m_col = 0; m_col < m_width; m_col++)
            {
                pVal +=  d_m_CM[m_row*m_width + m_col] * SM[m_row + ty][m_col + tx];
            }
        }

        //only save output value to output matrix when inside matrix bounds
        if(row_o<width && col_o<height){
            p[row_o*width + col_o] = pVal;
        }
    }
}

__global__ void gpu_2DConvolution_tiled_fixed(const unsigned int *n, unsigned int *p, const int width, const int height, const int m_width, const int m_height, const int interval, const int tile_size)
{
    __shared__ unsigned int SM[TILE_SIZE + M_W - 1][TILE_SIZE + M_H - 1];
 
    int tx = threadIdx.x;
    int ty = threadIdx.y;
    //output indexes
    int row_o = ty + blockIdx.y * tile_size;
    int col_o = tx + blockIdx.x * tile_size;
    //input indexes
    int row_i = row_o - interval;
    int col_i = col_o - interval;
    //tile indexes
    int row_t = ty + interval;
    int col_t = tx + interval;

    //LOAD TILE ELEMENTS
    //original tile
    SM[row_t][col_t] = n[row_o * width + col_o];
    //outside border
    if(ty == 0){    //top
        for(int outside_row=0; outside_row < interval; outside_row++){
            SM[outside_row][col_t] = row_i < 0? 0: n[(row_i + outside_row) * width + col_o];
        }
        if(tx == 0){    //top_left corner
            for(int outside_row = 0; outside_row < interval; outside_row ++ ){
                for(int outside_col = 0; outside_col < interval; outside_col ++ {
                    SM[outside_row][outside_col] = row_i < 0 ? 0 : n[(row_i + outside_row) * width + col_i + outside_col];
                }
            }
        }
    }
    else if(ty == tile_size - 1){ //bottom
        for(int outside_row=0; outside_row < interval; outside_row++){
            SM[row_t + outside_col][col_t] = row_o >= height? 0: n[(row_o + outside_row) * width + col_o];
        }
        if(tx == 0){    //bottom-right corner
            for(int outside_row = 0; outside_row < interval; outside_row ++ ){
                for(int outside_col = 0; outside_col < interval; outside_col ++ {
                    SM[row_t + outside_row][col_t + outside_col] = row_o >= height? 0: n[(row_o +outside_row) * width + col_o + outside_col];
                }
            }
        }
    }
    if(tx == 0){    //left
        for(int outside_col=0; outside_col < interval; outside_col++){
            SM[row_t][outside_col] = col_i < 0? n[row_o * width + col_i + outside_col];
        }
        if(tx == 0){    //bottom - left corner
            for(int outside_row = 0; outside_row < interval; outside_row ++ ){
                for(int outside_col = 0; outside_col < interval; outside_col ++ {
                    SM[row_t + outside_row][outside_col] = row_i<0 ? 0 : n[(row_o + outside_row) * width + col_i + outside_col];
                }
            }
        }
    }
    else if(tx == tile_size - 1){   //right
        for(int outside_col=0; outside_col < interval; outside_col++){
            SM[row_t][col_t + outside_col] = col_o >= height? n[row_o * width + col_o + outside_col];
        }
        if(ty == 0){    //top-right corner
            for(int outside_row = 0; outside_row < interval; outside_row ++ ){
                for(int outside_col = 0; outside_col < interval; outside_col ++ {
                    SM[outside_row][col_t + outside_col] = row_o>=height? 0: n[(row_i +outside_row) * width + col_o + outside_col];
                }
            }
        }
    }

    //wait till all tile elements are loaded
    __syncthreads();

    //output value
    int pVal = 0;

    //only compute if it is an output tile element
    if(ty < tile_size && tx < tile_size){
        for (int m_row = 0; m_row < m_height; m_row++)
        {
            for (int m_col = 0; m_col < m_width; m_col++)
            {
                pVal +=  d_m_CM[m_row*m_width + m_col] * SM[ty - interval + m_row][tx - interval + m_col];
            }
        }

        //only save output value to output matrix when inside matrix bounds
        if(row_o<width && col_o<height){
            p[row_o*width + col_o] = pVal;
        }
    }


}

// __global__ void gpu_2DConvolution_internal(const unsigned int *n, unsigned int *p, const int width, const int height, const int m_width, const int m_height, const int interval, const int tile_size)
// {
//     __shared__ unsigned int SM[TILE_SIZE][TILE_SIZE];

//     int tx = threadIdx.x;
//     int ty = threadIdx.y;
//     int row = ty + blockIdx.y * tile_size;
//     int col = tx + blockIdx.x * tile_size;
//     int row_start_point = row - interval;
//     int col_start_point = col - interval;

//     //load tile 
//     if(row < width && col < height){
//         SM[ty][tx] = n[row * width + col];
//     }

//     //wait till all tile elements are loaded
//     __syncthreads();

//     //output value
//     int pVal = 0;

//     //only compute if it is an output tile element
//     if(ty < tile_size && tx < tile_size){
//         for(int m_row=0; m_row<m_height; m_row++){
//             for(int m_col=0; m_col<m_width; m_col++){
//                 int row_i = row_start_point + m_row
//                 int col_i = col_start_point + m_col
//                 if(row_i >= 0 && row_i < height && col_i >= 0 && col_i<width){
//                     // if()
//                     // pVal += 
//                 }
//             }
//         }
//     }
// }
int main()
{
    srand((unsigned int)time(NULL));
    generateRandomValuesFor2D(N, N_W, N_H);
    generateRandomValuesFor2D(M, M_W, M_H);

    // MK: GPU Memory
    unsigned int *d_n, *d_m, *d_p_basic, *d_p_CM, *d_p_tiled;

    int n_byteSize = N_W * N_H * sizeof(unsigned int);
    int m_byteSize = M_W * M_H * sizeof(unsigned int);

    // allocate memory in device
    cudaError_t err = cudaMalloc((void **)&d_n, n_byteSize);
    checkCudaError(err);
    err = cudaMalloc((void **)&d_m, m_byteSize);
    checkCudaError(err);
    err = cudaMalloc((void **)&d_p_basic, n_byteSize);
    checkCudaError(err);
    err = cudaMalloc((void **)&d_p_CM, n_byteSize);
    checkCudaError(err);    
    err = cudaMalloc((void **)&d_p_tiled, n_byteSize);
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
    mkClockMeasure *ckGpu_tiled = new mkClockMeasure("GPU CODE - TILED");
    ckGpu_tiled->clockReset();

    for (int i = 0; i < MAX_ITER; i++)
    {
        ckCpu->clockResume();
        cpu_2DConvolution_basic(N, M, cpuOut, N_W, N_H, M_W, M_H, M_INTERVAL);
        ckCpu->clockPause();

        ckGpu->clockResume();
        gpu_2DConvolution_basic<<<gridSize, blockSize>>>(d_n, d_m, d_p_basic, N_W, N_H, M_W, M_H, M_INTERVAL);
        err = cudaDeviceSynchronize();
        ckGpu->clockPause();

        ckGpu_CM->clockResume();
        gpu_2DConvolution_cm<<<gridSize, blockSize>>>(d_n, d_p_CM, N_W, N_H, M_W, M_H, M_INTERVAL);
        err = cudaDeviceSynchronize();
        ckGpu_CM->clockPause();

        ckGpu_tiled->clockResume();
        gpu_2DConvolution_tiled<<<gridSize, blockSize>>>(d_n, d_p_tiled, N_W, N_H, M_W, M_H, M_INTERVAL, TILE_SIZE);
        err = cudaDeviceSynchronize();
        ckGpu_tiled->clockPause();

        checkCudaError(err);
    }

    err = cudaMemcpy(gpuOut_2D_basic, d_p_basic, n_byteSize, cudaMemcpyDeviceToHost);
    checkCudaError(err);
    err = cudaMemcpy(gpuOut_2D_CM, d_p_CM, n_byteSize, cudaMemcpyDeviceToHost);
    checkCudaError(err);
    err = cudaMemcpy(gpuOut_2D_tiled, d_p_tiled, n_byteSize, cudaMemcpyDeviceToHost);
    checkCudaError(err);

    cudaFree(d_n);
    cudaFree(d_m);
    cudaFree(d_p_basic);
    cudaFree(d_p_CM);
    cudaFree(d_p_tiled);


    if (compare2DMatrix(cpuOut, gpuOut_2D_basic, N_W, N_H))
    {
        printf("---------- [CPU] 2D convolution basic-------------\n");
        ckCpu->clockPrint();
        printf("\n---------- [GPU] 2D convolution basic-------------\n");
        ckGpu->clockPrint();
        printf("\n---------- [GPU] 2D convolution constant memory-------------\n");
        ckGpu_CM->clockPrint();
        printf("\n---------- [GPU] 2D convolution tiled-------------\n");
        ckGpu_tiled->clockPrint();
    }
    else
    {
        printf("ERROR: Two Matrices are not same\n");
    }

    // print2DMatrix(cpuOut, N_W, N_H);
    // print2DMatrix(gpuOut_2D_basic, N_W, N_H);
}
