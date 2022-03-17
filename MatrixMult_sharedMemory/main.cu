/*
 * =====================================================================================
 *
 *       Filename:  main.cu
 *
 *    Description: 	Matrix Multiplication using shared memory
 *
 *        Version:  1.0
 *        Created:  2022/02/01 10:07:38
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

const int A_H = 500;
const int A_W = 500;
const int B_H = A_W;
const int B_W = 500;
const int WIDTH = 500;
const int TILE_WIDTH = 16;
const unsigned int MAX_NUM = 10;
const int MAX_ITER = 1;

unsigned int matrixA[A_H * A_W];
unsigned int matrixB[B_H * B_W];
unsigned int cpuOut[A_H * B_W];
unsigned int gpuOut[A_H * B_W];

void generateRandomValues(unsigned int *input, const int rowSize, const int colSize)
{
	for (int i = 0; i < rowSize; i++)
	{
		for (int j = 0; j < colSize; j++)
		{
			input[i * colSize + j] = (unsigned int)float(rand()) / float(RAND_MAX) * MAX_NUM;
		}
	}
}

void printMatrixValue(const unsigned int *input, const int rowSize, const int colSize)
{
	printf("Print Matrix \n -----------\n");
	for (int i = 0; i < rowSize; i++)
	{
		for (int j = 0; j < colSize; j++)
		{
			printf("%u\t", input[i * colSize + j]);
		}
		printf("\n");
	}
	printf("--------\n");
}

bool compareMatrix(const unsigned int *inputA, const unsigned int *inputB, const int rowSize, const int colSize)
{
	bool ret = true;
	for (int i = 0; i < rowSize * colSize; i++)
	{
		if (inputA[i] != inputB[i])
		{
			ret = false;
			break;
		}
	}
	return ret;
}

void cpuMatrixMul(const unsigned int *h_a, const unsigned int *h_b, unsigned int *h_c, const int width)
{
	for (int r = 0; r < width; r++)
	{
		for (int c = 0; c < width; c++)
		{
			int temp = 0;
			for (int i = 0; i < width; i++)
			{
				temp += h_a[r * width + i] * h_b[i * width + c];
			}
			h_c[r * width + c] = temp;
		}
	}
}

__global__ void gpuMatrixMul(unsigned int *d_a, unsigned int *d_b, unsigned int *d_c, const int width)
{
	__shared__ int Asm[TILE_WIDTH][TILE_WIDTH];
	__shared__ int Bsm[TILE_WIDTH][TILE_WIDTH];

	// create automatic variables for faster access
	int bx = blockIdx.x;
	int by = blockIdx.y;
	int tx = threadIdx.x;
	int ty = threadIdx.y;

	// result array index
	int row = by * TILE_WIDTH + ty;
	int col = bx * TILE_WIDTH + tx;

	int cVal = 0;
	for (int tile = 0; tile < int(ceil((float)width / TILE_WIDTH)); tile++)
	{
		if (row < width && (tile * TILE_WIDTH + tx) < width)
			Asm[ty][tx] = d_a[row * width + tile * TILE_WIDTH + tx];
		else
			Asm[ty][tx] = 0;
		if (col < width && (tile * TILE_WIDTH + ty) < width)
			Bsm[ty][tx] = d_b[(tile * TILE_WIDTH + ty) * width + col];
		else
			Bsm[ty][tx] = 0;

		__syncthreads();

		for (int i = 0; i < TILE_WIDTH; i++)
		{
			cVal += Asm[ty][i] * Bsm[i][tx];
		}
		__syncthreads();
	}
	if (row < width && col < width)
		d_c[row * width + col] = cVal;
}

int main()
{
	srand((unsigned int)time(NULL));
	generateRandomValues(matrixA, WIDTH, WIDTH);
	generateRandomValues(matrixB, WIDTH, WIDTH);

	// MK: GPU Memory
	unsigned int *d_a, *d_b, *d_c;
	size_t matrixSizeA = sizeof(unsigned int) * WIDTH * WIDTH;
	size_t matrixSizeB = sizeof(unsigned int) * WIDTH * WIDTH;
	size_t matrixSizeC = sizeof(unsigned int) * WIDTH * WIDTH;

	// allocate memory in device
	cudaError_t err = cudaMalloc((void **)&d_a, matrixSizeA);
	checkCudaError(err);
	err = cudaMalloc((void **)&d_b, matrixSizeB);
	checkCudaError(err);
	err = cudaMalloc((void **)&d_c, matrixSizeC);
	checkCudaError(err);

	err = cudaMemcpy(d_a, matrixA, matrixSizeA, cudaMemcpyHostToDevice);
	checkCudaError(err);
	err = cudaMemcpy(d_b, matrixB, matrixSizeB, cudaMemcpyHostToDevice);
	checkCudaError(err);

	// MK: Thread Num
	const int threads = 16;
	dim3 gridSize(ceil((float)(WIDTH) / (float)threads), ceil((float)(WIDTH) / (float)threads));
	dim3 blockSize(threads, threads);

	// MK: Time Measurement
	mkClockMeasure *ckCpu = new mkClockMeasure("CPU CODE");
	ckCpu->clockReset();

	mkClockMeasure *ckGpu = new mkClockMeasure("GPU CODE");
	ckGpu->clockReset();

	for (int i = 0; i < MAX_ITER; i++)
	{
		ckCpu->clockResume();
		cpuMatrixMul(matrixA, matrixB, cpuOut, WIDTH);
		ckCpu->clockPause();

		ckGpu->clockResume();
		gpuMatrixMul<<<gridSize, blockSize>>>(d_a, d_b, d_c, WIDTH);
		err = cudaDeviceSynchronize();
		ckGpu->clockPause();
		checkCudaError(err);
	}

	err = cudaMemcpy(gpuOut, d_c, matrixSizeC, cudaMemcpyDeviceToHost);
	checkCudaError(err);

	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_c);

	if (compareMatrix(cpuOut, gpuOut, A_H, B_W))
	{
		ckCpu->clockPrint();
		ckGpu->clockPrint();
	}
	else
	{
		printf("ERROR: Two Matrices are not same\n");
	}

	// printMatrixValue(matrixA, A_H, A_W);
	// printMatrixValue(matrixB, B_H, B_W);
	// printMatrixValue(cpuOut, A_H, B_W);
	// printMatrixValue(gpuOut, A_H, B_W);
}
