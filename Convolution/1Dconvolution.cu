/*
 * =====================================================================================
 *
 *       Filename:  main.cu
 *
 *    Description: 1D Convolution
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

const int N_WIDTH = 65536;
const int M_WIDTH = 5;
const int M_INTERVAL = M_WIDTH / 2;
const int MAX_NUM = 10;
const int MAX_ITER = 100;

const int tbSize = N_WIDTH;
const int threads = 32;
const int TILE_SIZE = threads;

unsigned int n[N_WIDTH];
unsigned int m[M_WIDTH];

unsigned int cpuOut[N_WIDTH];
unsigned int gpuOut_1D_basic[N_WIDTH];
unsigned int gpuOUt_1D_CM[N_WIDTH];
unsigned int gpuOUt_1D_tiled[N_WIDTH];
unsigned int gpuOUt_1D_tiled_internal[N_WIDTH];
__constant__ unsigned int d_m_CM[M_WIDTH];

void generateRandomValuesFor1D(unsigned int *input, const int width)
{
	for (int i = 0; i < width; i++)
	{
		input[i] = (unsigned int)float(rand()) / float(RAND_MAX) * MAX_NUM;
	}
}

void print1DMask(const unsigned int *input, const int width)
{
	printf("Print 1D Mask \n -----------\n");
	for (int i = 0; i < width; i++)
	{
		printf("%u\t", input[i]);
	}

	printf("\n--------\n");
}

void print1DMatrix(const unsigned int *input, const int width)
{
	printf("Print Matrix \n -----------\n");
	for (int j = 0; j < 25; j++)
	{
		for (int i = 0; i < 16; i++)
		{
			printf("%u\t", input[j * 16 + i]);
		}
		printf("\n");
	}
	printf("\n--------\n");
}

bool compare1DMatrix(const unsigned int *inputA, const unsigned int *inputB, const int width)
{
	bool ret = true;
	for (int i = 0; i < width; i++)
	{
		if (inputA[i] != inputB[i])
		{
			ret = false;
			break;
		}
	}
	return ret;
}

void cpu_1DConvolution_basic(const unsigned int *n, unsigned int *m, unsigned int *p, const int width, const int m_width, const int interval)
{
	for (int i = 0; i < width; i++)
	{
		int pVal = 0;
		for (int j = i - interval, k = 0; k < m_width; j++, k++)
		{
			if (j >= 0 && j < width)
			{
				pVal += n[j] * m[k];
			}
		}
		p[i] = pVal;
	}
}

__global__ void gpu_1DConvolution_basic(const unsigned int *n, unsigned int *m, unsigned int *p, const int width, const int m_width, const int interval)
{
	int index = blockDim.x * blockIdx.x + threadIdx.x;

	int pVal = 0;
	for (int j = index - interval, i = 0; i < m_width; j++, i++)
	{
		if (j >= 0 && j < width)
		{
			pVal += n[j] * m[i];
		}
	}
	p[index] = pVal;
}

__global__ void gpu_1DConvolution_CM(const unsigned int *n, unsigned int *p, const int width, const int m_width, const int interval)
{
	int index = blockDim.x * blockIdx.x + threadIdx.x;

	int pVal = 0;
	for (int j = index - interval, i = 0; i < m_width; j++, i++)
	{
		if (j >= 0 && j < width)
		{
			pVal += n[j] * d_m_CM[i];
		}
	}
	p[index] = pVal;
}

__global__ void gpu_1DConvolution_tiled(const unsigned int *n, unsigned int *p, const int width, const int m_width, const int interval)
{
	__shared__ unsigned int n_sm[TILE_SIZE + M_WIDTH - 1];

	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int halo_index_left = blockDim.x * (blockIdx.x - 1) + threadIdx.x;
	if (threadIdx.x >= blockDim.x - interval)
	{
		n_sm[threadIdx.x - (blockDim.x - interval)] = halo_index_left < 0 ? 0 : n[halo_index_left];
	}

	n_sm[threadIdx.x + interval] = n[index];

	int halo_index_right = blockDim.x * (blockIdx.x + 1) + threadIdx.x;
	if (threadIdx.x < interval)
	{
		n_sm[interval + TILE_SIZE + threadIdx.x] = halo_index_right < width ? n[halo_index_right] : 0;
	}
	__syncthreads();

	int pVal = 0;
	for (int i = 0; i < m_width; i++)
	{
		pVal += n_sm[threadIdx.x + i] * d_m_CM[i];
	}
	p[index] = pVal;
}

__global__ void gpu_1DConvolution_tiled_internal(const unsigned int *n, unsigned int *p, const int width, const int m_width, const int interval)
{
	__shared__ unsigned int n_sm[TILE_SIZE];

	int index = blockDim.x * blockIdx.x + threadIdx.x;

	n_sm[threadIdx.x] = n[index];
	__syncthreads();

	int pVal = 0;
	int thisTileStartIndex = blockIdx.x * blockDim.x;
	int nextTileStartIndex = (blockIdx.x - 1) * blockDim.x;
	int nStartIndex = index - interval;
	for (int i = 0; i < m_width; i++)
	{
		int nIndex = nStartIndex + i;
		if (nIndex >= 0 && nIndex < width)
		{
			if (nIndex >= thisTileStartIndex && nIndex < nextTileStartIndex)
			{
				pVal += n_sm[threadIdx.x + i - interval] * d_m_CM[i];
			}
			else
			{
				pVal += n[nIndex] * d_m_CM[i];
			}
		}
	}
	p[index] = pVal;
}


int main()
{
	srand((unsigned int)time(NULL));
	generateRandomValuesFor1D(n, N_WIDTH);
	generateRandomValuesFor1D(m, M_WIDTH);

	// MK: GPU Memory
	unsigned int *d_n, *d_m, *d_p_basic, *d_p_CM, *d_p_tiled, *d_p_tiled_internal;

	int n_byteSize = N_WIDTH * sizeof(unsigned int);
	int m_byteSize = M_WIDTH * sizeof(unsigned int);

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
	err = cudaMalloc((void **)&d_p_tiled_internal, n_byteSize);
	checkCudaError(err);

	err = cudaMemcpy(d_n, n, n_byteSize, cudaMemcpyHostToDevice);
	checkCudaError(err);
	err = cudaMemcpy(d_m, m, m_byteSize, cudaMemcpyHostToDevice);
	checkCudaError(err);
	err = cudaMemcpyToSymbol(d_m_CM, m, m_byteSize);
	checkCudaError(err);

	// MK: Thread Num
	dim3 gridSize(tbSize / threads, 1, 1);
	dim3 blockSize(threads, 1, 1);

	// MK: Time Measurement
	mkClockMeasure *ckCpu = new mkClockMeasure("CPU CODE");
	ckCpu->clockReset();

	mkClockMeasure *ckGpu_basic = new mkClockMeasure("GPU BASIC CODE");
	ckGpu_basic->clockReset();

	mkClockMeasure *ckGpu_CM = new mkClockMeasure("GPU CM CODE");
	ckGpu_CM->clockReset();

	mkClockMeasure *ckGpu_tiled = new mkClockMeasure("GPU tiled CODE");
	ckGpu_tiled->clockReset();

	mkClockMeasure *ckGpu_tiled_internal = new mkClockMeasure("GPU tiled internal CODE");
	ckGpu_tiled_internal->clockReset();

	for (int i = 0; i < MAX_ITER; i++)
	{
		ckCpu->clockResume();
		cpu_1DConvolution_basic(n, m, cpuOut, N_WIDTH, M_WIDTH, M_INTERVAL);
		ckCpu->clockPause();

		ckGpu_basic->clockResume();
		gpu_1DConvolution_basic<<<gridSize, blockSize>>>(d_n, d_m, d_p_basic, N_WIDTH, M_WIDTH, M_INTERVAL);
		err = cudaDeviceSynchronize();
		ckGpu_basic->clockPause();

		ckGpu_CM->clockResume();
		gpu_1DConvolution_CM<<<gridSize, blockSize>>>(d_n, d_p_CM, N_WIDTH, M_WIDTH, M_INTERVAL);
		err = cudaDeviceSynchronize();
		ckGpu_CM->clockPause();
		checkCudaError(err);

		ckGpu_tiled->clockResume();
		gpu_1DConvolution_tiled<<<gridSize, blockSize>>>(d_n, d_p_tiled, N_WIDTH, M_WIDTH, M_INTERVAL);
		err = cudaDeviceSynchronize();
		ckGpu_tiled->clockPause();
		checkCudaError(err);

		ckGpu_tiled_internal->clockResume();
		gpu_1DConvolution_tiled_internal<<<gridSize, blockSize>>>(d_n, d_p_tiled_internal, N_WIDTH, M_WIDTH, M_INTERVAL);
		err = cudaDeviceSynchronize();
		ckGpu_tiled_internal->clockPause();
		checkCudaError(err);
	}

	err = cudaMemcpy(gpuOut_1D_basic, d_p_basic, n_byteSize, cudaMemcpyDeviceToHost);
	checkCudaError(err);
	err = cudaMemcpy(gpuOUt_1D_CM, d_p_CM, n_byteSize, cudaMemcpyDeviceToHost);
	checkCudaError(err);
	err = cudaMemcpy(gpuOUt_1D_tiled, d_p_tiled, n_byteSize, cudaMemcpyDeviceToHost);
	checkCudaError(err);
	err = cudaMemcpy(gpuOUt_1D_tiled_internal, d_p_tiled_internal, n_byteSize, cudaMemcpyDeviceToHost);
	checkCudaError(err);

	cudaFree(d_n);
	cudaFree(d_m);
	cudaFree(d_p_basic);
	cudaFree(d_p_CM);
	cudaFree(d_p_tiled);
	cudaFree(d_p_tiled_internal);

	if (compare1DMatrix(cpuOut, gpuOUt_1D_tiled_internal, N_WIDTH))
	{
		printf("---------- [CPU] 1D convolution basic-------------\n");
		ckCpu->clockPrint();
		printf("\n---------- [GPU] 1D convolution basic-------------\n");
		ckGpu_basic->clockPrint();
		printf("\n---------- [GPU] 1D convolution with Constant Memory-------------\n");
		ckGpu_CM->clockPrint();
		printf("\n---------- [GPU] 1D convolution tiled-------------\n");
		ckGpu_tiled->clockPrint();
		printf("\n---------- [GPU] 1D convolution tiled internal-------------\n");
		ckGpu_tiled_internal->clockPrint();
	}
	else
	{
		printf("ERROR: Two Matrices are not same\n");
	}
	// print1DMatrix(n, N_WIDTH);
	// print1DMask(m, M_WIDTH);
	// print1DMatrix(cpuOut, N_WIDTH);

	// print1DMatrixValue(gpuOut_1D_basic, N_WIDTH);
}
