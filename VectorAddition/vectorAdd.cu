/*
 * =====================================================================================
 *
 *       Filename:  vectorAdd.cu
 *
 *    Description:  CH02 Example Code
 *
 *        Version:  1.0
 *        Created:  01/11/2022
 *       Revision:  none
 *       Compiler:  nvcc
 *
 *         Author:  Yeonhee Jung
 *
 * =====================================================================================
 */

const int MAX_SIZE=1000;
const float MAX_NUM=100.0;
const int MAX_ITER= 1000;

#include <iostream>

float inputA[MAX_SIZE];
float inputB[MAX_SIZE];
float gpuAns[MAX_SIZE];
float cpuAns[MAX_SIZE];

using namespace std;

#define checkCudaError(error) 			\
	if(error != cudaSuccess){ 				\
		printf("%s in %s at line %d\n", \
				cudaGetErrorString(error), 	\
				__FILE__ ,__LINE__); 				\
		exit(EXIT_FAILURE);							\
	}

void generateRandomValues(float *array, float max, const int size){
	for(int i = 0; i < size; i++){
		array[i] = float(rand())/float(RAND_MAX) * max;
	}
}

void cpuVectorAddition(float *h_a, float *h_b, float *h_c, const int size){
	for(int i = 0; i < size; i++){
		h_c[i] = h_a[i] + h_b[i];
	}
}

__global__
void gpuVectorAddition(float *d_a, float *d_b, float *d_c, const int size){
	int tId = blockDim.x * blockIdx.x + threadIdx.x;
	if(tId < size){
		d_c[tId] = d_a[tId] + d_b[tId];
	}
}

void checkAnswer(float *h_a, float *d_a, const int size){
	bool isSame = true;
	for(int i = 0; i < size; i++){
		if(h_a[i] != d_a[i]){
			cout<<"-\tERROR: IDX - "<< i << " (" << h_a[i] << " != " << d_a[i] << " )" << endl;
			isSame = false;
		}
	}
	if(isSame)
		printf("All values are same\n");
	else
		printf("Some values are not same\n");
}

//MK: Main Function
int main(){
	srand((unsigned int)time(NULL));
	
	//MK: Random
	generateRandomValues(inputA, MAX_NUM, MAX_SIZE);
	generateRandomValues(inputB, MAX_NUM, MAX_SIZE);
	
	//MK: GPU Memory
	float *d_a, *d_b, *d_c;
	int arraySize = MAX_SIZE * sizeof(float);
	cudaError_t err = cudaMalloc((void **) &d_a, arraySize);
	checkCudaError(err);
	err = cudaMalloc((void **) &d_b, arraySize);
	checkCudaError(err);
	err = cudaMalloc((void **) &d_c, arraySize);
	checkCudaError(err);

	err = cudaMemcpy(d_a, inputA, arraySize, cudaMemcpyHostToDevice);
	checkCudaError(err);
	err = cudaMemcpy(d_b, inputB, arraySize, cudaMemcpyHostToDevice);
	checkCudaError(err);

	const int tSize = 256;
	dim3 gridSize(ceil((float)MAX_SIZE/(float)tSize), 1, 1);
	dim3 blockSize(tSize, 1, 1);
	
	for(int i = 0; i < MAX_ITER; i++){
		cpuVectorAddition(inputA, inputB, cpuAns, MAX_SIZE);

		gpuVectorAddition<<<gridSize, blockSize>>>(d_a, d_b, d_c, MAX_SIZE);
		err=cudaDeviceSynchronize();
		checkCudaError(err);
	}

	err = cudaMemcpy(gpuAns, d_c, arraySize, cudaMemcpyDeviceToHost);
	checkCudaError(err);

	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_c);

	checkAnswer(cpuAns, gpuAns, MAX_SIZE);

	return 0;
}
