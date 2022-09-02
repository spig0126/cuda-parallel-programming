/*
 * =====================================================================================
 *
 *       Filename:  DFT.cu
 *
 *    Description:  DFT
 *
 *        Version:  1.0
 *        Created:  07/20/2022 
 *       Revision:  none
 *       Compiler:  nvcc
 *
 *         Author:  Yeonhee Jung
 *   Organization:  EWHA Womans Unversity
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "comp.cuh"
#include "mkClockMeasure.h"
#include "mkCuda.h"


const int SAMPLING_RATE = 131072;
const int N = SAMPLING_RATE;
const int FREQ_NUM = 3;
const int MAX_ITER = 1;
int db_bytesize = N * sizeof(double);
int comp_bytesize = N * 2 * sizeof(double);

double freq_amp_ph[FREQ_NUM][3] = {{1, 3, 0}, {4, 1, 0}, {7, 0.5, 0}};
double sample_points[N], freq[N], sig[N];
//cpu
Comp x[N];
double idft_sig[N], amp[N];
//gpu
Comp gpu_x[N];
double gpu_idft_sig[N], gpu_amp[N];

// thread, grid size */ 
int thread = 64;
int tbSize = N/thread;
dim3 gridSize(tbSize, 1, 1);
dim3 blockSize(thread, 1, 1);

void create_sample_points(double* sample_points){
    for (int i = 0; i<N; i++){
        sample_points[i] = (double) i/SAMPLING_RATE;
        freq[i] = i;

        // printf("%lf\n", sample_points[i]);
    }
}

void generate_sig(double* sample_points, double* sig){
    for (int s_i = 0; s_i < N; s_i++){
        sig[s_i] = 0;
        for (int f_i = 0; f_i < FREQ_NUM; f_i++)
        {
            sig[s_i] += freq_amp_ph[f_i][1] * sin(2 * M_PI * freq_amp_ph[f_i][0] * sample_points[s_i]);
        }
    }
}

bool compareResult(double* a, double* b, int size){
    double epsilon = 0.000001;

    for(int i=0; i<size; i++){
        if(fabs(a[i] - b[i]) > epsilon){
            // printf("a; %lf, \t b: %lf\n", a[i], b[i]);
            return false;
        }
    }
    return true;
}
void printData(double* x, double* y, int size, const char* x_label, const char* y_label){
    for(int i=0; i<size; i++){
        printf("%s: %lf\t%s: %lf\n", x_label, x[i], y_label, y[i]);
    }
}

void save_data(double* x_values, double* y_values, const char* path){
    FILE *dataf = fopen(path, "w");
    for (int i=0; i<N; i++){
        fprintf(dataf, "%lf %lf\n", x_values[i], y_values[i]);
    }
}

void cpu_dft(mkClockMeasure* ck){
    ck->clockResume();

    double exp = 2 * M_PI / N;

    for (int k = 0; k<N; k++){
        x[k].i = 0;
        x[k].r = 0;
        for (int n = 0; n < N; n++)
        {
            double bn = exp * k * n;
            x[k].i -= sig[n] * sin(bn);
            x[k].r += sig[n] * cos(bn);
        }
        amp[k] = 2 * sqrt(pow(x[k].r, 2) + pow(x[k].i, 2)) / N;
        // printf("frequency: %f\tamplitude: %f\n", freq[k], x[k]);
        // printf("[CPU] k: %d\txi: %lf\t, xr: %lf\n", k, xi[k], xr[k]);

    }
    ck->clockPause();
}

__global__ void gpu_dft(Comp* d_x, double* d_amp, double* d_sig, int N){
    //1 freq = 1 thread
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    double xi=0, xr=0;
    double exp = 2 * M_PI * k / N;


    for(int n=0; n<N; n++){
        double bn = exp * n;
        xi -= d_sig[n] * sin(bn);
        xr += d_sig[n] * cos(bn);
    }
    d_x[k].i = xi;
    d_x[k].r = xr;
    d_amp[k] =  2 * sqrt(pow(xr, 2) + pow(xi, 2)) / N;

    // printf("frequency: %d\tamplitude: %lf\t sig: %f\n", k, d_x[k], d_sig[k]);
    // printf("k: %d\n", k);
}

void cal_gpu_dft(mkClockMeasure* ck_mem_transfer, mkClockMeasure* ck_kernel, mkClockMeasure* ck_exec){
    // allocate device memory 
    Comp *d_x;
    double *d_amp, *d_sig;
    cudaError_t err = cudaMallocHost((void**)&d_x, comp_bytesize);
    checkCudaError(err);
    err = cudaMallocHost((void**)&d_amp, db_bytesize);
    checkCudaError(err);
    err = cudaMallocHost((void**)&d_sig, db_bytesize);

    ck_exec->clockResume();

    // memory transfer: host -> device
    ck_mem_transfer->clockResume();
    err = cudaMemcpy(d_sig, sig, db_bytesize, cudaMemcpyHostToDevice);
    checkCudaError(err);
    ck_mem_transfer->clockPause();


    // launch kernel
    ck_kernel->clockResume();
    gpu_dft<<<gridSize, blockSize>>>(d_x, d_amp, d_sig, N); 
    ck_kernel->clockPause();

    // memory transfer: device -> host 
    ck_mem_transfer->clockResume();
    err = cudaMemcpy(gpu_amp, d_amp, db_bytesize, cudaMemcpyDeviceToHost);
    checkCudaError(err);
    err = cudaMemcpy(gpu_x, d_x, comp_bytesize, cudaMemcpyDeviceToHost);
    ck_mem_transfer->clockPause();

    ck_exec->clockPause();

    // free device memory
    cudaFree(d_amp);
    cudaFree(d_sig);
}

void cpu_idft(mkClockMeasure *ck){
    ck->clockResume();
    for (int n = 0; n < N; n++)
    {
        idft_sig[n] = 0;
        for (int k = 0; k < N; k++)
        {
            double bn = 2 * M_PI * k * n / N;
            idft_sig[n] += x[k].r * cos(bn) + x[k].i * sin(bn);
        }
        idft_sig[n] /= N;
        // printf("t: %d\tamplitude: %lf\n", n, idft_sig[n]);
    }
    ck->clockPause();

}

__global__ void gpu_idft(double* d_idft_sig, Comp* d_x, int N){
    // 1 signal at sample point n = 1 thread
    int n = blockIdx.x * blockDim.x + threadIdx.x;
    double exp = 2 * M_PI * n / N;
    double sig = 0;
    

    for(int k=0; k<N; k++){
        double bn = exp * k;
        sig += d_x[k].r * cos(bn) + d_x[k].i * sin(bn);
    }

    d_idft_sig[n] = sig / N;   
    // printf("t: %d\tamplitude: %lf\n", n, d_idft_sig[n]);
    // printf("[GPU] k: %d\txi: %lf\t, xr: %lf\n", n, d_xi[n], d_xr[n]);
}

void cal_gpu_idft(mkClockMeasure* ck_mem_transfer, mkClockMeasure* ck_kernel, mkClockMeasure* ck_exec){
    // allocate device memory 
    Comp *d_x;
    double *d_idft_sig;
    cudaError_t err = cudaMallocHost((void**)&d_idft_sig, db_bytesize);
    err = cudaMallocHost((void**)&d_x, comp_bytesize);
    checkCudaError(err);

    ck_exec->clockResume();

    // memory transfer: host -> device
    ck_mem_transfer->clockResume();
    err = cudaMemcpy(d_x, gpu_x, comp_bytesize, cudaMemcpyHostToDevice);
    ck_mem_transfer->clockPause();

    // launch kernel
    ck_kernel->clockResume();
    gpu_idft<<<gridSize, blockSize>>>(d_idft_sig, d_x, N); 
    ck_kernel->clockPause();

    // memory transfer: device -> host 
    ck_mem_transfer->clockResume();
    err = cudaMemcpy(gpu_idft_sig, d_idft_sig, db_bytesize, cudaMemcpyDeviceToHost);
    checkCudaError(err);
    ck_mem_transfer->clockPause();

    ck_exec->clockPause();

    // ree device memory 
    cudaFree(d_idft_sig);
}



int main(void){
    // int sampling_rates = [64, 256, 1024, 4096]
    printf("SAMPLING_RATE : %d\nMAX_ITERATION : %d\n\n", SAMPLING_RATE, MAX_ITER);

    /* set basic data */
    create_sample_points(sample_points);
    generate_sig(sample_points, sig);

    /* create clocks */
    mkClockMeasure *ckCpu_dft = new mkClockMeasure("CPU - DFT");
    mkClockMeasure *ckCpu_idft = new mkClockMeasure("CPU - IDFT");
    mkClockMeasure *ckGpu_dft_mem_transfer = new mkClockMeasure("MEMORY TRANSFER");
    mkClockMeasure *ckGpu_dft_kernels = new mkClockMeasure("KERNELS");
    mkClockMeasure *ckGpu_dft = new mkClockMeasure("GPU - DFT TOTAL");
    mkClockMeasure *ckGpu_idft_mem_transfer = new mkClockMeasure("MEMORY TRANSFER");
    mkClockMeasure *ckGpu_idft_kernels = new mkClockMeasure("KERNELS");
    mkClockMeasure *ckGpu_idft = new mkClockMeasure("GPU - IDFT TOTAL");
    ckCpu_dft->clockReset(), ckCpu_idft->clockReset();
    ckGpu_dft_mem_transfer->clockReset(), ckGpu_dft_kernels->clockReset(), ckGpu_dft->clockReset();
    ckGpu_idft_mem_transfer->clockReset(), ckGpu_idft_kernels->clockReset(), ckGpu_idft->clockReset();

    for(int i=0; i<MAX_ITER; i++){
        // /* CPU - DFT */
        // cpu_dft(ckCpu_dft);

        // /* CPU - IDFT */
        // cpu_idft(ckCpu_idft);

        /* GPU - DFT */
        cal_gpu_dft(ckGpu_dft_mem_transfer, ckGpu_dft_kernels, ckGpu_dft);

        /* GPU - IDFT */  
        cal_gpu_idft(ckGpu_idft_mem_transfer, ckGpu_idft_kernels, ckGpu_idft);
    }




    /* print DFT performance */
    // if(compareResult(amp, gpu_amp, N)){
    if(true){
        printf("-------------------[CPU] DFT ---------------------\n");
		ckCpu_dft->clockPrint();
        printf("\n-------------------[GPU] DFT ---------------------\n");
        ckGpu_dft_mem_transfer->clockPrint();
        ckGpu_dft_kernels->clockPrint();
		ckGpu_dft->clockPrint();
    }
    else{
        printf("ERROR: DFT results are not the same\n\n");
    }

    /* print Inverse DFT performance */
    // if(compareResult(idft_sig, gpu_idft_sig, N)){
    if(true){
        printf("\n\n-------------------[CPU] Inverse DFT ---------------------\n");
        ckCpu_idft->clockPrint();
        printf("\n-------------------[GPU] Inverse DFT ---------------------\n");
        ckGpu_idft_mem_transfer->clockPrint();
        ckGpu_idft_kernels->clockPrint();
		ckGpu_idft->clockPrint();
    }
    else{
        printf("ERROR: Inverse DFT results are not the same\n\n");
    }

    /* save data */
    save_data(sample_points, sig,  "data/original_signal.txt");
    save_data(freq, amp, "data/cpu_dft_frequencies.txt");
    save_data(freq, gpu_amp, "data/gpu_dft_frequencies.txt");
    save_data(sample_points, idft_sig, "data/idft_signal.txt");
    save_data(sample_points, gpu_idft_sig, "data/gpu_idft_signal.txt");
}