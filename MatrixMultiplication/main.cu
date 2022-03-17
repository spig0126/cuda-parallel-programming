/*
 * =====================================================================================
 *
 *       Filename:  main.cu
 *
 *    Description:  Matrix Multiplication
 *
 *        Version:  1.0
 *        Created:  01/18/2022
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

const int A_H = 200;
const int A_W = 3000;
const int B_H = A_W;
const int B_W = 200;
const unsigned int MAX_NUM = 10;
const int MAX_ITER = 10;

unsigned int matrixA[A_H * A_W];
unsigned int matrixB[B_H * B_W];
unsigned int cpuOut[A_H * B_W];
unsigned int gpuOut[A_H * B_W];

void generateRandomValues(unsigned int *input, const int rowSize, const int colSize){
	for(int i = 0; i < rowSize; i++){
		for(int j = 0; j < colSize; j++){
			input[i * colSize + j] = (unsigned int) float(rand())/float(RAND_MAX) * MAX_NUM;
		}
	}
}

void printMatrixValue(const unsigned int *input, const int rowSize, const int colSize){
	printf("Print Matrix \n -----------\n");
	for(int i = 0; i < rowSize; i++){
		for(int j = 0; j < colSize; j++){
			printf("%u\t", input[i * colSize + j]);
		}
		printf("\n");
	}
	printf("--------\n");
}

bool compareMatrix(const unsigned int *inputA, const unsigned int *inputB, const int rowSize, const int colSize){
	bool ret = true;
	for(int i = 0; i < rowSize * colSize; i++){
		if(inputA[i] != inputB[i]){
			ret = false;
			break;
		}
	}
	return ret;
}

void cpuMatrixMul(const unsigned int *h_a, const unsigned int *h_b, unsigned int *h_c, const int aRowSize, const int aColSize, const int bRowSize, const int bColSize){
	assert(aColSize == bRowSize);

	for (int r = 0; r < aRowSize; r++){
		for (int c = 0; c < bColSize; c++){
			int temp = 0;
			for (int i = 0; i < aColSize; i++)
			{
				temp += h_a[r*aColSize + i] * h_b[i*bColSize + c];
			}
			h_c[r * bColSize + c] = temp;
		}
	}
}

__global__
void gpuMatrixMul(unsigned int *d_a, unsigned int *d_b, unsigned int *d_c, const int aRowSize, const int aColSize, const int bRowSize, const int bColSize){
	assert(aColSize == bRowSize);

	int row = blockIdx.y * blockDim.y + threadIdx.y;
	int col = blockIdx.x * blockDim.x + threadIdx.x;

	if(row<aRowSize && col<bColSize){
		int index = row * bColSize + col;
		d_c[index] = 0;
		for (int i = 0; i < aColSize; i++){
			d_c[index] += d_a[row * aColSize + i] * d_b[i * bColSize + col];
			//intf("(%d/%d)=(%d/%d) %d\n", row, col, d_c[index], (row*aColSize+i), (i*bRowSize + col));
		}
	}
}

int main(){
	srand((unsigned int)time(NULL));
	generateRandomValues(matrixA, A_H, A_W);
	generateRandomValues(matrixB, B_H, B_W);

	//MK: GPU Memory 
	unsigned int *d_a, *d_b, *d_c;
	size_t matrixSizeA = sizeof(unsigned int) * A_H * A_W;
	size_t matrixSizeB = sizeof(unsigned int) * B_H * B_W;
	size_t matrixSizeC = sizeof(unsigned int) * A_H * B_W;
	
	//allocate memory in device
	cudaError_t err = cudaMalloc((void **) &d_a, matrixSizeA);
	checkCudaError(err);
	err = cudaMalloc((void **) &d_b, matrixSizeB);
	checkCudaError(err);
	err = cudaMalloc((void **) &d_c, matrixSizeC);
	checkCudaError(err);

	err = cudaMemcpy(d_a, matrixA, matrixSizeA, cudaMemcpyHostToDevice);
	checkCudaError(err);
	err = cudaMemcpy(d_b, matrixB, matrixSizeB, cudaMemcpyHostToDevice);
	checkCudaError(err);

	//MK: Thread Num
	//const int tbSize = 256;
	const int threads = 16;
	dim3 gridSize(ceil((float)(A_H) / (float)threads), ceil((float)(B_W) / (float)threads), 1);
	dim3 blockSize(threads, threads, 1);

	//MK: Time Measurement
	mkClockMeasure *ckCpu = new mkClockMeasure("CPU CODE");
	ckCpu->clockReset();
	
	mkClockMeasure *ckGpu = new mkClockMeasure("GPU CODE");
	ckGpu->clockReset();

	for(int i = 0; i < MAX_ITER; i++){
		ckCpu->clockResume();
		cpuMatrixMul(matrixA, matrixB, cpuOut, A_H, A_W, B_H, B_W);
		ckCpu->clockPause();

		ckGpu->clockResume();
		gpuMatrixMul<<<gridSize, blockSize>> >(d_a, d_b, d_c, A_H, A_W, B_H, B_W);
		err=cudaDeviceSynchronize();
		ckGpu->clockPause();
		checkCudaError(err);

	}

	err = cudaMemcpy(gpuOut, d_c, matrixSizeC, cudaMemcpyDeviceToHost);
	checkCudaError(err);

	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_c);

	if(compareMatrix(cpuOut, gpuOut, A_H, B_W)){
		ckCpu->clockPrint();
		ckGpu->clockPrint();
	}else{
		printf("ERROR: Two Matrices are not same\n");
	}

	/*printMatrixValue(matrixA, A_H, A_W);
	printMatrixValue(matrixB, B_H, B_W);
	printMatrixValue(cpuOut, A_H, B_W);
	printMatrixValue(gpuOut, A_H, B_W);*/
}

