/*
 * =====================================================================================
 *
 *       Filename:  basicAsyncProgram.cu
 *
 *    Description: to see host and device asynchronously operating
 *
 *        Version:  1.0
 *        Created:  04/11/22
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yeonhee Jung, spig0126@ewha.ac.kr
 *   Organization:  Ewha Womans University
 *
 * =====================================================================================
 */
#include "mkCuda.h"
#include "mkClockMeasure.h"

const int BLOCK_SIZE = 256;
const int GRID_SIZE = 65536;
const int MAX_ITER = 100;

const int N_SIZE = BLOCK_SIZE * GRID_SIZE;
const int N_BYTE_SIZE = N_SIZE * sizeof(int);
int N[N_SIZE];
int M1_gpuOut[N_SIZE];
int M2_gpuOut[N_SIZE];
int cpuOut[N_SIZE];
int gpuOut[N_SIZE];

const int N_STREAM = 4;
const int STREAM_SIZE = N_SIZE / N_STREAM ;
const int STREAM_BYTES = sizeof(int) * STREAM_SIZE;

bool compare(int *inputA, int *inputB, const int width)
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

void printArray(int *array, const int width, const char *des){
    printf("-------------%s------------\n", des);
    for(int b=0; b<GRID_SIZE; b++){
        for(int i=0; i<BLOCK_SIZE; i++){
            printf("%d\t", array[b*BLOCK_SIZE + i]);
        }
        printf("\n");
    }
}

void initializeArray(int *array, const int width){
    for(int i=0; i<width; i++){
        array[i] = i;
    }
}

void cpu_calSquare(int *h_n, const int width){
    for(int i=0; i<width; i++){
        cpuOut[i] = h_n[i] * h_n[i];
    }
}

__global__ void kernel_calSquare(int *d_n, int *d_out, const int start_index){
    int i = blockIdx.x * blockDim.x + threadIdx.x + start_index;    //start_index 필요성: GPU는 그저 지정된 grid block 내용을 기반으로 thread, block id를 제공해줄 뿐, 어디서부터 데이터 접근해야 할 지 모름!! 

    d_out[i] = d_n[i] * d_n[i];
}

void cpuOperation(int *h_n, const int width, mkClockMeasure *clk){
    //start clock
    clk -> clockResume();

    cpu_calSquare(h_n, width);

    //stop clock
    clk -> clockPause();
}

void gpu_sequential(int *h_n, int *gpuOut,int *d_n, int *d_out, mkClockMeasure *clk){
    clk -> clockResume();

    cudaMemcpy(d_n, h_n, N_BYTE_SIZE, cudaMemcpyHostToDevice);
    kernel_calSquare<<<GRID_SIZE, BLOCK_SIZE>>>(d_n, d_out, 0);
    cudaMemcpy(gpuOut, d_out,N_BYTE_SIZE, cudaMemcpyDeviceToHost);

    clk -> clockPause();
}

/* 
loop over all opertions by each chunk of array
(CHUNK host-to-device transfer -> CHUNK kernel launch -> CHUNK device-to-host transfer -> repeat)
*/
void gpu_loop_chunk(int *h_n, int *gpuOut,int *d_n, int *d_out, mkClockMeasure *clk){
    //start clock
    clk -> clockResume();

    //declare stream array and error
    cudaStream_t s[N_STREAM];
    cudaError_t e; 

    //cretae streams
    for(int i=0; i<N_STREAM; i++){
        e = cudaStreamCreate(&s[i]);
        checkCudaError(e);
    }

    for(int i=0; i<N_STREAM; i++){
        int chunkIndex = i * STREAM_SIZE;

        //CHUNK host-to-device transfer
        cudaMemcpyAsync(&d_n[chunkIndex], &h_n[chunkIndex], STREAM_BYTES, cudaMemcpyHostToDevice, s[i]);

        //CHUNK kernel launch -> grid 크기를 줄여 kernel 호출
        kernel_calSquare<<<STREAM_SIZE/BLOCK_SIZE, BLOCK_SIZE, 0, s[i]>>>(d_n, d_out, chunkIndex);

        //CHUNK device-to-host transfer
        cudaMemcpyAsync(&gpuOut[chunkIndex], &d_out[chunkIndex], STREAM_BYTES, cudaMemcpyDeviceToHost, s[i]);
    }

    //destory streams
    for(int i=0; i<N_STREAM; i++){
        e = cudaStreamDestroy(s[i]);
        checkCudaError(e);
    }
            
    //stop clock
    clk -> clockPause();
}

/* 
batch similar operations together
(ALL host-to-device transfer -> ALL kernel launches -> ALL device-to-host transfer)
*/
void gpu_batching(int *h_n, int *gpuOut,int *d_n, int *d_out, mkClockMeasure *clk){
    //clock start    
    clk -> clockResume();

    //declare stream array
    cudaStream_t s2[N_STREAM];

    cudaError_t e; 
    //cretae streams
    for(int i=0; i<N_STREAM; i++){
        e = cudaStreamCreate(&s2[i]);
        checkCudaError(e);
    }

    //host-to-device transfer
    for(int i=0; i<N_STREAM; i++){
        int chunkIndex = i * STREAM_SIZE;
        cudaMemcpyAsync(&d_n[chunkIndex], &h_n[chunkIndex], STREAM_BYTES, cudaMemcpyHostToDevice, s2[i]);
    }

    //kernel launch 
    for(int i=0; i<N_STREAM; i++){
        int chunkIndex = i * STREAM_SIZE;
        kernel_calSquare<<<STREAM_SIZE/BLOCK_SIZE, BLOCK_SIZE, 0, s2[i]>>>(d_n, d_out, chunkIndex);
    }

    //device-to-host transfer
    for(int i=0; i<N_STREAM; i++){
        int chunkIndex = i * STREAM_SIZE;
        cudaMemcpyAsync(&gpuOut[chunkIndex], &d_out[chunkIndex], STREAM_BYTES, cudaMemcpyDeviceToHost, s2[i]);
    }

    //destory streams
    for(int i=0; i<N_STREAM; i++){
        e = cudaStreamDestroy(s2[i]);
        checkCudaError(e);
    }

    //clock stop
    clk -> clockPause();

}


int main(void) {
    //initialize memory variables
    int *d_n, *d_out, *d_out_m1, *d_out_m2, *h_n;

    //allocate devic memory
    cudaError_t e = cudaMalloc((void **)&d_n, N_BYTE_SIZE);
	checkCudaError(e);
    e = cudaMalloc((void **)&d_out, N_BYTE_SIZE);
	checkCudaError(e);
    e = cudaMalloc((void **)&d_out_m1, N_BYTE_SIZE);
    checkCudaError(e);
    e = cudaMalloc((void **)&d_out_m2, N_BYTE_SIZE);
    checkCudaError(e);

    //pinned memory
    e = cudaMallocHost((void **)&h_n, N_BYTE_SIZE);

    //initialize arrays
    initializeArray(N, N_SIZE);
    initializeArray(h_n, N_SIZE);

    //clocks
    mkClockMeasure *ckCpu = new mkClockMeasure("CPU CODE");
    ckCpu -> clockReset();
    mkClockMeasure *ckGpu = new mkClockMeasure("GPU CODE");
    ckGpu -> clockReset();
    mkClockMeasure *M1_ckGpu = new mkClockMeasure("GPU CODE - METHOD #1 loop chuncks");
    M1_ckGpu -> clockReset();
    mkClockMeasure *M2_ckGpu = new mkClockMeasure("GPU CODE - METHOD #2 batch");
    M2_ckGpu -> clockReset();


    for(int i=0; i<MAX_ITER; i++){
        cpuOperation(N, N_SIZE, ckCpu);

        gpu_sequential(N, gpuOut, d_n, d_out, ckGpu);
        e = cudaDeviceSynchronize();
        checkCudaError(e);

        gpu_loop_chunk(N, M1_gpuOut, d_n, d_out_m1, M1_ckGpu);
        e = cudaDeviceSynchronize();
        checkCudaError(e);

        gpu_batching(N, M2_gpuOut, d_n, d_out_m2, M2_ckGpu);
    }
    

    //deallocate device memory
    cudaFree(d_n);
    cudaFree(d_out);
    cudaFree(d_out_m1);
    cudaFree(d_out_m2);

    //compare cpu and gpu outputs
    if(compare(cpuOut, M1_gpuOut, N_SIZE)){
        printf("Asynchronous programming complete\n\n");
        printf("----------------- [CPU] SQUARE -----------------\n");
        ckCpu -> clockPrint();
        printf("\n----------------- [GPU] SQUARE -----------------\n");
        ckGpu -> clockPrint();
        printf("\n----------------- [GPU] LOOP CHUNK -----------------\n");
        M1_ckGpu -> clockPrint();
        printf("\n----------------- [GPU] BATCHING -----------------\n");
        M2_ckGpu -> clockPrint();
    }
    else{
        printf("Outputs are not equal!!\n");
        // printArray(cpuOut, N_SIZE, "CPU");
        // printf("\n");
        // printArray(gpuOut, N_SIZE, "GPU");
        // printf("\n");
        // printArray(M1_gpuOut, N_SIZE, "GPU method 1");
    }
}