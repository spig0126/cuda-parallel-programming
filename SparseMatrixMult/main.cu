/*
 * =====================================================================================
 *
 *       Filename:  main.cu
 *
 *    Description: 	Sparse Matrix Multiplication
 *
 *        Version:  1.0
 *        Created:  2022/04/01
 *       Revision:  none
 *       Compiler:  nvcc
 *
 *         Author:  Yeonhee Jung
 *   Organization:  Ewha Womans University
 *
 * =====================================================================================
 */

#include <assert.h>
#include "sparseMatrix.h"
#include "mkCuda.h"
#include "mkClockMeasure.h"

const unsigned int MAX_NUM = 10;
const int MAX_ITER = 100;

const unsigned int SIZE = 256;
unsigned int X[SIZE];
unsigned int Y[SIZE];

unsigned int cpuOut_basic[SIZE];
unsigned int cpuOut_CSR[SIZE];
unsigned int gpuOut_basic[SIZE];
unsigned int gpuOut_CSR[SIZE];
unsigned int gpuOut_ELL[SIZE];
unsigned int gpuOut_COO[SIZE];

int a_byteSize = SIZE * SIZE * sizeof(unsigned int);
int vector_byteSize = SIZE * sizeof(unsigned int);
int ELL_byteSize = SIZE * 4 * sizeof(unsigned int);

//Thread Num
const int tbSize = SIZE;
const int threads = 32;
dim3 gridSize(ceil((float)tbSize / (float)threads),1, 1);
dim3 blockSize(threads, 1, 1);

void generateRandomVector(unsigned int *v, const int size)
{
	for (int i = 0; i < size; i++)
	{
		v[i] = (unsigned int)float(rand()) / float(RAND_MAX) * MAX_NUM;
	}
}

void printVector(unsigned int *v, const int size, char *description){
	printf("Print vector: %s \n -----------\n", description);
	for (int i = 0; i < size; i++)
	{
		printf("%d\t", v[i]);
	}
	printf("\n--------\n");
}

bool compareMatrix(const unsigned int *inputA, const unsigned int *inputB, const int size){
	bool ret = true;
	for(int i = 0; i < size; i++){
		if(inputA[i] != inputB[i]){
			ret = false;
			break;
		}
	}
	return ret;
}

void cpuSparseMatrixMult_basic(unsigned int a[SIZE][SIZE], unsigned int *x, unsigned int *y, unsigned int *y_out, const int rowSize, const int colSize)
{
	for(int row=0; row < rowSize; row ++){
		int dot = 0;

		for(int col=0; col <colSize; col++){
			dot += a[row][col] * x[col];
		}
		y_out[row] = y[row] + dot;
	}
}

void cpuSparseMatrixMult_CSR(unsigned int *x, unsigned int *y, unsigned int *data, unsigned int *col_index, unsigned int *row_ptr, unsigned int *y_out, const int data_size, const int ptr_size, const int num_rows)
{
	for(int row=0; row<num_rows; row++){
		int dot =0;
		int row_start = row_ptr[row];
		int row_end = row_ptr[row+1];
		for(int elem = row_start; elem<row_end; elem++){
			dot += data[elem] * x[col_index[elem]];
		}
		y_out[row] = y[row] + dot;
	}
}

__global__ void gpuSparseMatrixMult_basic(const unsigned int *a, const unsigned int *x, unsigned int *y, unsigned int *y_out, const int colSize){
	int row = blockDim.x * blockIdx.x + threadIdx.x;
	int dot=0;
	for(int col=0; col<colSize; col++){
		dot += a[row*colSize + col] * x[col];
	}
	y_out[row] = y[row] + dot;
}

__global__ void gpuSparseMatrixMult_CSR(unsigned int *x, unsigned int *y, unsigned int *y_out, unsigned int *data, unsigned int *col_index, unsigned int *row_ptr, const int size, const int ptr_size){
	int row = blockIdx.x * blockDim.x + threadIdx.x;

	if(row < ptr_size){
		int dot =0;
		int row_start = row_ptr[row];
		int row_end = row_ptr[row+1];

		for(int elem=row_start; elem<row_end; elem++){
			dot += data[elem] * x[col_index[elem]];
		}

		y_out[row] = y[row] + dot;
	}
}

__global__ void gpuSparseMatrixMult_ELL(unsigned int *x, unsigned int *y, unsigned int *y_out, unsigned int *transposed, unsigned int *transposed_col_index, const int row_size, const int col_size)
{
	int row = blockIdx.x * blockDim.x + threadIdx.x;

	if(row < row_size){
		int dot = 0;

		for(int col = 0; col<col_size; col++){
			dot += transposed[col*row_size + row] * x[transposed_col_index[col*row_size + row]];
		}

		y_out[row] = y[row] + dot;
	}
}

int convertToCooFormat(unsigned int *m, unsigned int *data, unsigned int *col_index, unsigned int *row_index, const int num_col, const int num_row){
	int index = 0;
	for(int col=0; col<num_col; col++){
		for(int row=0; row<num_row; row++){
			if(m[row*num_row + col] != 0){
				data[index] = m[row*num_row + col];
				col_index[index] = col;
				row_index[index] = row;
				index++;
			}
		}
	}

	return index;
}

void convertToEllFormat(unsigned int m_data[SIZE][4], unsigned int m_col_index[SIZE][4], unsigned int *data, unsigned int *col_index, const int num_row, const int num_col){
	for(int row=0; row<num_row; row++){
		for(int col=0; col<num_col; col++){
			data[col * num_row + row] = m_data[row][col];
			col_index[col * num_row + row] = m_col_index[row][col];
		}
	}
}

void hybridSparseMatrixMult(unsigned int *x, unsigned int *y, unsigned int padding[SIZE][4], unsigned int padding_col_index[SIZE][4], unsigned int transposed[4][SIZE], unsigned int *d_y_out_COO ,int coo_num_col, int ell_num_col, int row_size, int col_size){
	//row의 non zero element 최고 개수 구하기
	int max = 0, cnt=0;
	for(int row=0; row<SIZE; row++){
		cnt = 0;
		for(int col=0; col<SIZE; col++){
			if(A[row][col] != 0){
				cnt++;
			}
		}
		if(cnt>max){
			max = cnt;
		}
	}

	//ELL과 COO로 처리하기 위한 데이터들
	unsigned int *ELL_data = (unsigned int*)malloc(sizeof(unsigned int)*ell_num_col*SIZE);
	unsigned int *ELL_col_index = (unsigned int*)malloc(sizeof(unsigned int)*ell_num_col*SIZE);

	unsigned int *COO_matrix = (unsigned int*)malloc(sizeof(unsigned int)*coo_num_col*SIZE);
	unsigned int *COO_data = (unsigned int*)malloc(sizeof(unsigned int)*coo_num_col*SIZE);
	unsigned int *COO_col_index = (unsigned int*)malloc(sizeof(unsigned int)*coo_num_col*SIZE);
	unsigned int *COO_row_index = (unsigned int*)malloc(sizeof(unsigned int)*coo_num_col*SIZE);


	//ELL로 계산할 matrix 앞부분을 따로 COO 형태로 변환
	convertToEllFormat(padding, padding_col_index, ELL_data, ELL_col_index, SIZE, ell_num_col);

	//COO로 계산할 matrix 뒷부분을 따로 COO 형태로 변환
	for(int col=max - coo_num_col, coo_col=0; col<max; col++, coo_col++){
		for(int row=0; row<max; row++){
			COO_matrix[coo_col * col_size + row] = padding[row][col];
		}
	}
	int num_elem = convertToCooFormat(COO_matrix, COO_data, COO_col_index, COO_row_index, coo_num_col, SIZE);


	//ELL format으로 앞부분 행렬 곱셈 구하기
	unsigned int *d_hybrid_data, *d_hybrid_col_index;
	int byteSize = SIZE * ell_num_col * sizeof(unsigned int);

	cudaError_t err = cudaMalloc((void **) &d_hybrid_data, byteSize);
	checkCudaError(err);
	err = cudaMalloc((void **) &d_hybrid_col_index, byteSize);
	checkCudaError(err);

	err = cudaMemcpy(d_hybrid_data, ELL_data, byteSize, cudaMemcpyHostToDevice);
	checkCudaError(err);
	err = cudaMemcpy(d_hybrid_col_index, ELL_col_index, byteSize, cudaMemcpyHostToDevice);
	checkCudaError(err);

	gpuSparseMatrixMult_ELL<<<gridSize, blockSize>>>(x, y, d_y_out_COO, d_hybrid_data, d_hybrid_col_index, SIZE, ell_num_col);

	err=cudaDeviceSynchronize();
	checkCudaError(err);

	cudaFree(d_hybrid_data); 
	cudaFree(d_hybrid_col_index);

	//COO format으로 뒷부분 행렬 곱셈 구하고 y에 더하기
	for(int i=0; i<num_elem; i++){
		gpuOut_COO[COO_row_index[i]] += COO_data[i] * x[COO_col_index[i]];
	}
}

int main(){
	srand((unsigned int)time(NULL));
	generateRandomVector(X, SIZE);
	generateRandomVector(Y, SIZE);

	//MK: GPU Memory 
	unsigned int *d_a, *d_x, *d_y, *d_y_out_basic, *d_y_out_CSR, *d_y_out_COO, *d_y_out_ELL, *d_data, *d_col_index, *d_row_ptr, *d_row_index, *d_padding, *d_padding_col_index, *d_transposed, *d_transposed_col_index;

	//allocate memory in device
	cudaError_t err = cudaMalloc((void **) &d_a, a_byteSize);
	checkCudaError(err);
	err = cudaMalloc((void **) &d_x, vector_byteSize);
	checkCudaError(err);
	err = cudaMalloc((void **) &d_y, vector_byteSize);
	checkCudaError(err);
	err = cudaMalloc((void **) &d_y_out_basic, vector_byteSize);
	checkCudaError(err);
	err = cudaMalloc((void **) &d_y_out_CSR, vector_byteSize);
	checkCudaError(err);
	err = cudaMalloc((void **) &d_y_out_COO, vector_byteSize);
	checkCudaError(err);
	err = cudaMalloc((void **) &d_y_out_ELL, vector_byteSize);
	checkCudaError(err);
	err = cudaMalloc((void **) &d_data, vector_byteSize);
	checkCudaError(err);
	err = cudaMalloc((void **) &d_col_index, vector_byteSize);
	checkCudaError(err);
	err = cudaMalloc((void **) &d_row_index, vector_byteSize);
	checkCudaError(err);
	err = cudaMalloc((void **) &d_row_ptr, vector_byteSize);
	checkCudaError(err);
	err = cudaMalloc((void **) &d_padding, ELL_byteSize);
	checkCudaError(err);
	err = cudaMalloc((void **) &d_transposed, ELL_byteSize);
	checkCudaError(err);
	err = cudaMalloc((void **) &d_padding_col_index, ELL_byteSize);
	checkCudaError(err);
	err = cudaMalloc((void **) &d_transposed_col_index, ELL_byteSize);
	checkCudaError(err);

	//copy host memory to device memory
	err = cudaMemcpy(d_a, A, a_byteSize, cudaMemcpyHostToDevice);
	checkCudaError(err);
	err = cudaMemcpy(d_x, X, vector_byteSize, cudaMemcpyHostToDevice);
	checkCudaError(err);
	err = cudaMemcpy(d_y, Y, vector_byteSize, cudaMemcpyHostToDevice);
	checkCudaError(err);
	err = cudaMemcpy(d_data, data, vector_byteSize, cudaMemcpyHostToDevice);
	checkCudaError(err);
	err = cudaMemcpy(d_col_index, col_index, vector_byteSize, cudaMemcpyHostToDevice);
	checkCudaError(err);
	err = cudaMemcpy(d_row_index, row_index, vector_byteSize, cudaMemcpyHostToDevice);
	checkCudaError(err);
	err = cudaMemcpy(d_row_ptr, row_ptr, vector_byteSize, cudaMemcpyHostToDevice);
	checkCudaError(err);
	err = cudaMemcpy(d_transposed, transposed, ELL_byteSize, cudaMemcpyHostToDevice);
	checkCudaError(err);
	err = cudaMemcpy(d_transposed_col_index, transposed_col_index, ELL_byteSize, cudaMemcpyHostToDevice);
	checkCudaError(err);


	//MK: Time Measurement
	mkClockMeasure *ckCpu_basic = new mkClockMeasure("CPU CODE - BASIC");
	ckCpu_basic->clockReset();
	mkClockMeasure *ckCpu_CSR = new mkClockMeasure("CPU CODE - CSR");
	ckCpu_CSR->clockReset();
	mkClockMeasure *ckGpu_basic = new mkClockMeasure("GPU CODE - BASIC");
	ckGpu_basic->clockReset();
	mkClockMeasure *ckGpu_CSR = new mkClockMeasure("GPU CODE - CSR");
	ckGpu_CSR->clockReset();
	mkClockMeasure *ckGpu_ELL = new mkClockMeasure("GPU CODE - ELL");
	ckGpu_ELL->clockReset();
	mkClockMeasure *ckHybrid = new mkClockMeasure("GPU CODE - hybrid(ELL + COO)");
	ckHybrid->clockReset();

	for(int i = 0; i < MAX_ITER; i++){
		ckCpu_basic->clockResume();
		cpuSparseMatrixMult_basic(A, X, Y, cpuOut_basic, SIZE, SIZE);
		ckCpu_basic->clockPause();

		ckCpu_CSR->clockResume();
		cpuSparseMatrixMult_CSR(X, Y, data, col_index, row_ptr, cpuOut_CSR, SIZE, SIZE + 1, SIZE);
		ckCpu_CSR->clockPause();

		ckGpu_basic->clockResume();
		gpuSparseMatrixMult_basic<<<gridSize, blockSize>>>(d_a, d_x, d_y, d_y_out_basic, SIZE);
		err=cudaDeviceSynchronize();
		ckGpu_basic->clockPause();

		ckGpu_CSR->clockResume();
		gpuSparseMatrixMult_CSR<<<gridSize, blockSize>>>(d_x, d_y, d_y_out_CSR, d_data, d_col_index, d_row_ptr, SIZE, SIZE+1);
		err=cudaDeviceSynchronize();
		ckGpu_CSR->clockPause();

		ckGpu_ELL->clockResume();
		gpuSparseMatrixMult_ELL<<<gridSize, blockSize>>>(d_x, d_y, d_y_out_ELL, d_transposed, d_transposed_col_index, SIZE, 4);
		err=cudaDeviceSynchronize();
		ckGpu_ELL->clockPause();

		checkCudaError(err);

		// ckHybrid->clockResume();
		// hybridSparseMatrixMult(X, Y, padding, padding_col_index, transposed, d_y_out_COO ,1, 3, SIZE, 4);
		// ckHybrid->clockPause();

	}

	err = cudaMemcpy(gpuOut_basic, d_y_out_basic, vector_byteSize, cudaMemcpyDeviceToHost);
	checkCudaError(err);
	err = cudaMemcpy(gpuOut_CSR, d_y_out_CSR, vector_byteSize, cudaMemcpyDeviceToHost);
	checkCudaError(err);
	err = cudaMemcpy(gpuOut_ELL, d_y_out_ELL, vector_byteSize, cudaMemcpyDeviceToHost);
	checkCudaError(err);
	err = cudaMemcpy(gpuOut_COO, d_y_out_COO, vector_byteSize, cudaMemcpyDeviceToHost);
	checkCudaError(err);

	cudaFree(d_a); 
	cudaFree(d_x); 
	cudaFree(d_y); 
	cudaFree(d_y_out_basic); 
	cudaFree(d_y_out_CSR); 
	cudaFree(d_y_out_COO); 
	cudaFree(d_y_out_ELL); 
	cudaFree(d_data); 
	cudaFree(d_col_index); 
	cudaFree(d_row_ptr); 
	cudaFree(d_row_index); 
	cudaFree(d_transposed);
	cudaFree(d_transposed_col_index);

	if(compareMatrix(cpuOut_basic, gpuOut_ELL, SIZE)){
		printf("---------- [CPU] basic-------------\n");
		ckCpu_basic->clockPrint();

		printf("\n---------- [CPU] CSR-------------\n");
		ckCpu_basic->clockPrint();

		printf("\n---------- [GPU] basic-------------\n");
		ckGpu_basic->clockPrint();

		printf("\n---------- [GPU] CSR-------------\n");
		ckGpu_CSR->clockPrint();

		printf("\n---------- [GPU] ELL-------------\n");
		ckGpu_ELL->clockPrint();
	}else{
		printf("ERROR: Two Matrices are not same\n");
	}

	// printVector(cpuOut_basic, SIZE, "CPU basic");
	// printVector(gpuOut_CSR, SIZE, "GPU basic");
}
