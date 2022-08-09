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
#include "comp.h"
#include "mkClockMeasure.h"
#include "mkCuda.h"


const int SAMPLING_RATE = 4096;
const int N = SAMPLING_RATE;
const int FREQ_NUM = 3;
const int MAX_ITER = 1;


double freq_amp_ph[FREQ_NUM][3] = {{1, 3, 0}, {4, 1, 0}, {7, 0.5, 0}};
double sample_points[N], freq[N], sig[N];
double idft_sig[N], xr[N], xi[N], x[N];
double gpu_idft_sig[N], gpu_xr[N], gpu_xi[N], gpu_x[N];

__constant__ double d_sig[N];

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
    double epsilon = 0.000001f;

    for(int i=0; i<size; i++){
        if(fabs(a[i] - b[i]) < epsilon){
            // printf("a; %lf, \t b: %lf\n", a[i], b[i]);
            return true;
        }
    }
    return false;
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
        xi[k] = 0;
        xr[k] = 0;
        for (int n = 0; n < N; n++)
        {
            double bn = exp * k * n;
            xi[k] -= sig[n] * sin(bn);
            xr[k] += sig[n] * cos(bn);
        }
        x[k] = 2 * sqrt(pow(xr[k], 2) + pow(xi[k], 2)) / N;
        // printf("frequency: %f\tamplitude: %f\n", freq[k], x[k]);
        // printf("[CPU] k: %d\txi: %lf\t, xr: %lf\n", k, xi[k], xr[k]);

    }
    ck->clockPause();
}

__global__ void gpu_dft(double* d_xi, double* d_xr, double* d_x, int N){
    //1 freq = 1 thread
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    double xi=0, xr=0;
    double exp = 2 * M_PI * k / N;


    for(int n=0; n<N; n++){
        double bn = exp * n;
        xi -= d_sig[n] * sin(bn);
        xr += d_sig[n] * cos(bn);
    }
    d_xi[k] = xi;
    d_xr[k] = xr;
    d_x[k] =  2 * sqrt(pow(xr, 2) + pow(xi, 2)) / N;
    // printf("frequency: %d\tamplitude: %lf\t sig: %f\n", k, d_x[k], d_sig[k]);
    // printf("k: %d\n", k);

}

void cpu_idft(mkClockMeasure *ck){
    ck->clockResume();
    for (int n = 0; n < N; n++)
    {
        idft_sig[n] = 0;
        for (int k = 0; k < N; k++)
        {
            double bn = 2 * M_PI * k * n / N;
            idft_sig[n] += xr[k] * cos(bn) + xi[k] * sin(bn);
        }
        idft_sig[n] /= N;
        // printf("t: %d\tamplitude: %lf\n", n, idft_sig[n]);
    }
    ck->clockPause();

}

__global__ void gpu_idft(double* d_idft_sig, double* d_xi, double* d_xr, int N){
    // 1 signal at sample point n = 1 thread
    int n = blockIdx.x * blockDim.x + threadIdx.x;
    double exp = 2 * M_PI * n / N;
    double sig = 0;
    

    for(int k=0; k<N; k++){
        double bn = exp * k;
        sig += d_xr[k] * cos(bn) + d_xi[k] * sin(bn);
    }

    d_idft_sig[n] = sig / N;

    // printf("t: %d\tamplitude: %lf\n", n, d_idft_sig[n]);
    // printf("[GPU] k: %d\txi: %lf\t, xr: %lf\n", n, d_xi[n], d_xr[n]);
}


int main(void){
    // int sampling_rates = [64, 256, 1024, 4096]
    printf("SAMPLING_RATE : %d\nMAX_ITERATION : %d\n\n", SAMPLING_RATE, MAX_ITER);

    /* set basic data */
    create_sample_points(sample_points);
    generate_sig(sample_points, sig);

    /* create clocks */
    mkClockMeasure *ckCpu_dft = new mkClockMeasure("CPU - DFT"),  *ckCpu_idft = new mkClockMeasure("CPU - IDFT"), *ckGpu_dft = new mkClockMeasure("GPU - DFT"), *ckGpu_idft = new mkClockMeasure("GPU - IDFT");
    ckCpu_dft->clockReset(), ckCpu_idft->clockReset(), ckGpu_dft->clockReset(), ckGpu_idft->clockReset();

    /* allocate device memory */ 
    double *d_idft_sig, *d_xi, *d_xr, *d_x;
    int bytesize = N * sizeof(double);

    cudaError_t err = cudaMalloc((void**)&d_idft_sig, bytesize);
    checkCudaError(err);
    err = cudaMalloc((void**)&d_xi, bytesize);
    checkCudaError(err);
    err = cudaMalloc((void**)&d_xr, bytesize);
    checkCudaError(err);
    err = cudaMalloc((void**)&d_x, bytesize);
    checkCudaError(err);
    err = cudaMalloc((void**)&d_idft_sig, bytesize);
    checkCudaError(err);

    /* set thread, grid size */ 
    int thread = 256;
    int tbSize = N/thread;
    dim3 gridSize(tbSize, 1, 1);
    dim3 blockSize(thread, 1, 1);


    for(int i=0; i<MAX_ITER; i++){
        /* CPU - DFT */
        cpu_dft(ckCpu_dft);

        /* GPU - DFT */
        ckGpu_dft->clockResume();

        // memory transfer: host -> device
        err = cudaMemcpyToSymbol(d_sig, sig, bytesize); //constant memory
        checkCudaError(err);

        gpu_dft<<<gridSize, blockSize>>>(d_xi, d_xr, d_x, N); 

        // memory transfer: device -> host 
        err = cudaMemcpy(gpu_xi, d_xi, bytesize, cudaMemcpyDeviceToHost);
        checkCudaError(err);
        err = cudaMemcpy(gpu_xr, d_xr, bytesize, cudaMemcpyDeviceToHost);
        checkCudaError(err);
        err = cudaMemcpy(gpu_x, d_x, bytesize, cudaMemcpyDeviceToHost);
        checkCudaError(err);

        ckGpu_dft->clockPause();


        /* CPU - IDFT */
        cpu_idft(ckCpu_idft);

        /* GPU - IDFT */
        ckGpu_idft->clockResume();

        gpu_idft<<<gridSize, blockSize>>>(d_idft_sig, d_xi, d_xr, N); 

        // memory transfer: device -> host 
        err = cudaMemcpy(gpu_idft_sig, d_idft_sig, bytesize, cudaMemcpyDeviceToHost);
	    checkCudaError(err);

        ckGpu_idft->clockPause();
    }


    /* memory: device -> host */
    err = cudaMemcpy(gpu_xi, d_xi, bytesize, cudaMemcpyDeviceToHost);
	checkCudaError(err);
    err = cudaMemcpy(gpu_xr, d_xr, bytesize, cudaMemcpyDeviceToHost);
	checkCudaError(err);
    err = cudaMemcpy(gpu_x, d_x, bytesize, cudaMemcpyDeviceToHost);
	checkCudaError(err);
    err = cudaMemcpy(gpu_idft_sig, d_idft_sig, bytesize, cudaMemcpyDeviceToHost);
	checkCudaError(err);

    cudaFree(d_sig);
    cudaFree(d_idft_sig);
    cudaFree(d_xi);
    cudaFree(d_sig);
    cudaFree(d_x);

    /* print DFT performance */
    if(compareResult(x, gpu_x, N)){
        printf("-------------------[CPU] DFT ---------------------\n");
		ckCpu_dft->clockPrint();
        printf("-------------------[GPU] DFT ---------------------\n");
		ckGpu_dft->clockPrint();
    }
    else{
        printf("ERROR: DFT results are not the same\n\n");
    }

    /* print Inverse DFT performance */
    if(compareResult(idft_sig, gpu_idft_sig, N)){
        printf("\n-------------------[CPU] Inverse DFT ---------------------\n");
        ckCpu_idft->clockPrint();
        printf("-------------------[GPU] Inverse DFT ---------------------\n");
		ckGpu_idft->clockPrint();
    }
    else{
        printf("ERROR: Inverse DFT results are not the same\n\n");
    }

    /* save data */
    save_data(sample_points, sig,  "data/original_signal.txt");
    save_data(freq, x, "data/cpu_dft_frequencies.txt");
    save_data(freq, gpu_x, "data/gpu_dft_frequencies.txt");
    save_data(sample_points, idft_sig, "data/idft_signal.txt");
    save_data(sample_points, gpu_idft_sig, "data/gpu_idft_signal.txt");
}