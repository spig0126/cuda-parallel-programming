/*
 * =====================================================================================
 *
 *       Filename:  FFT.cu
 *
 *    Description:  FFT
 *
 *        Version:  1.0
 *        Created:  07/28/2022 
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
#include "comp.cuh"
#include "mkCuda.h"
#include "mkClockMeasure.h"

const int SAMPLING_RATE = 4096;
const int N = SAMPLING_RATE;
const int FREQ_NUM = 3;
const int MAX_ITER = 1;

double freq_amp_ph[FREQ_NUM][3] = {{1, 3, 0}, {4, 1, 0}, {7, 0.5, 0}};
double sample_points[N], freq[N], sig[N];
__constant__ double d_sig[N];
/* cpu */
double fft_amp_rec[N], fft_amp_iter[N];
Comp fft_res_rec[N], fft_res_iter[N];
/* gpu */
double gpu_amp[N];
Comp gpu_res[N];
int db_bytesize = N * sizeof(double);
int comp_bytesize = N * 2 * sizeof(double);

/* thread, grid size */ 
int thread = 256;    // block size = half length of signal (시그널의 주기성 때문)
int tbSize = N/thread;
dim3 half_gridSize(tbSize/2, 1, 1);
dim3 total_gridSize(tbSize, 1, 1);
dim3 blockSize(thread, 1, 1);


/* basic functions */
void db_to_comp(double* a, Comp* b, int size){
    for(int i=0; i<size; i++){
        b[i].r = a[i];
        b[i].i = 0;
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

void save_data(double* x_values, double* y_values, const char* path){
    FILE *dataf = fopen(path, "w");
    for (int i=0; i<N; i++){
    fprintf(dataf, "%lf %lf\n", x_values[i], y_values[i]);
    }
}

/* data settings */
void create_sample_points(double* sample_points){
    for (int i = 0; i<N; i++){
        sample_points[i] = (double)i/SAMPLING_RATE;
        freq[i] = i;
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

void generate_sig_TEST(double* sig){
    for(int i=0; i< N; i++){
        sig[i] = i;
    }
}


//radix-2 cooley-tukey fft
void cpu_fft_recursive(int len, Comp* x){ //x = signal 값
    if(len == 1){
        return;
    }
    int half_len = len >> 1;

    //divide x into 2 subgroups: even, odd
    Comp even[half_len], odd[half_len];
    for(int i=0; i<half_len; i++){
        even[i] = x[2 * i];
        odd[i] = x[2 * i + 1];
    }

    //divde
    cpu_fft_recursive(half_len, even);
    cpu_fft_recursive(half_len, odd);

    //conquer & combine
    for(int k=0; k<half_len; k++){

        Comp w_k = cal_euler(-2 * M_PI * k / len);
        Comp t = comp_mult(w_k, odd[k]);

        x[k] = comp_add(even[k], t);
        x[k + half_len] = comp_sub(even[k], t);
    }
}

void cpu_cal_fft_recursive(double* sample_points, double* sig, mkClockMeasure* ck){
    //initialize fft as signal values
    db_to_comp(sig, fft_res_rec, N);
    
    //start clock
    ck -> clockResume();

    //calculate fft
    cpu_fft_recursive(N, fft_res_rec);

    //cal magnitude of each frequency
    for(int k=0; k<N; k++){
        fft_amp_rec[k] = 2 * sqrt(pow(fft_res_rec[k].r, 2) + pow(fft_res_rec[k].i, 2)) / N;
    }

    ck->clockPause();
}

void cpu_fft_iterative(int len, Comp* x){
    int depth = (int)log2(len);
    int half_len = len/2;
    Comp x_copy[N];

    // 1 iteration = 1 layer
    for(int l=2, d=1; l<=len; l<<=1, d<<=1){
        // update data from lower layer
        for(int i=0; i<N; i++){
            x_copy[i].r = x[i].r;
            x_copy[i].i = x[i].i;
        }

        int itvl = len / l;

        for(int i=0; i<itvl; i++){
            for(int k=0; k<d; k++){
                int temp = k * itvl;
                int data_even_i = i + temp*2;
                int data_odd_i = data_even_i + itvl;
                int res_i = i + temp;
                int res_half_i = res_i + half_len;

                Comp w_k = cal_euler(-2 * M_PI * k / l);
                Comp t = comp_mult(w_k, x_copy[data_odd_i]);
                Comp even = x_copy[data_even_i];
                x[res_i] = comp_add(even, t);
                x[res_half_i] = comp_sub(even, t);
            }
        }
    }
}

void cpu_cal_fft_iterative(double* sample_points, double* sig, mkClockMeasure* ck){
    //initialize fft as signal values
    db_to_comp(sig, fft_res_iter, N);
    
    //start clock
    ck -> clockResume();

    //calculate fft
    cpu_fft_iterative(N, fft_res_iter);

    //cal magnitude of each frequency
    for(int k=0; k<N; k++){
        fft_amp_iter[k] = 2 * sqrt(pow(fft_res_iter[k].r, 2) + pow(fft_res_iter[k].i, 2)) / N;
    }

    ck->clockPause();
}

__global__ void gpu_cal_amp(Comp* d_res, double* d_amp, int half_len){
    int k = blockDim.x * blockIdx.x + threadIdx.x;

    d_amp[k] = 2 * sqrt(pow(d_res[k].r, 2) + pow(d_res[k].i, 2)) / N;
}

__global__ void gpu_update_data(Comp* d_res, Comp* d_res_copy, int len, int temp){
    /* update data from lower layer */
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    // if(temp == 4){
    //     printf("i: %d => d_res_copy[i].r = %lf,  d_res_copy[i].i = %lf, d_res[i].r = %lf, d_res[i].i = %lf\n", i, d_res_copy[i].r, d_res_copy[i].i, d_res[i].r,  d_res[i].i);
    // }

    d_res_copy[i].r = d_res[i].r;
    d_res_copy[i].i = d_res[i].i;
}
__global__ void gpu_fft(Comp* d_res, Comp* d_res_copy, int len, int half_len, int itvl){
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    int k = (int)(i / itvl);
    int data_even_i = i % itvl + k * itvl * 2;
    int data_odd_i = data_even_i + itvl;
    int res_i = i;
    int res_half_i = i + half_len;
  
    Comp w_k = cal_euler(-2 * M_PI * k / len);
    Comp t = comp_mult(w_k, d_res_copy[data_odd_i]);
    Comp even = d_res_copy[data_even_i];
    d_res[res_i] = comp_add(even, t);
    d_res[res_half_i] = comp_sub(even, t);
 
    // if(len == 2){
    //     // printf("k: %d \n", k);
    //     // printf("\tbefore => \n");
    //     // printf("\t\tx[k]: %lf, x[k + half_len]: %lf\n",  x[i + k*itvl].r, x[i + k*itvl + half_len].r);
    //     // printf("\t\teven: %lf, odd: %lf \n", x_copy[i + k*itvl*2].r, x_copy[i + k*itvl*2 + itvl].r);

    //     printf("\t\ti: %d, k: %d\t=>\tx[k]: %lf, x[k + half_len]: %lf\n",  i, k, d_res[res_i].r, d_res[res_half_i].r);
    //     printf("\t\ti: %d, k: %d\t=>\teven: %lf, odd: %lf \n", i, k, d_res_copy[data_even_i].r, d_res_copy[data_odd_i].r);

    // }
}


void gpu_cal_fft(double* sample_points, double* sig, mkClockMeasure* ck_mem_transfer, mkClockMeasure* ck_kernels, mkClockMeasure* ck_total){
    /* initialize fft as signal values */
    db_to_comp(sig, gpu_res, N);

    /* allocate device memory */ 
    Comp *d_res, *d_res_copy;
    double *d_amp;

    cudaError_t e = cudaMalloc((void**)&d_res, comp_bytesize);
    checkCudaError(e);
    e = cudaMalloc((void**)&d_res_copy, comp_bytesize);
    checkCudaError(e);
    e = cudaMalloc((void**)&d_amp, db_bytesize);
    checkCudaError(e);

    ck_total -> clockResume();

    /* data transfer - host -> dev */
    ck_mem_transfer -> clockResume();

    // ck_mem_transfer -> clockResume();
    e = cudaMemcpy(d_res, gpu_res, comp_bytesize, cudaMemcpyHostToDevice);
    checkCudaError(e);
    ck_mem_transfer -> clockPause();

    /* cal Fourier efficients */
    ck_kernels -> clockResume();
    int depth = (int)log2(N);
    int half_len = N / 2;
    for(int l=2; l <= N; l<<=1){ 
        int itvl = N / l;

        gpu_update_data<<<total_gridSize, blockSize>>>(d_res, d_res_copy, N, l);
        e=cudaDeviceSynchronize();

        gpu_fft<<<half_gridSize, blockSize>>>(d_res, d_res_copy, l, half_len, itvl);
        e=cudaDeviceSynchronize();

	    checkCudaError(e);
    }

    gpu_cal_amp<<<total_gridSize, blockSize>>>(d_res, d_amp, half_len);
    ck_kernels -> clockPause();


    /* data transfer - dev -> host */
    ck_mem_transfer -> clockResume();
    e = cudaMemcpy(gpu_amp, d_amp, db_bytesize, cudaMemcpyDeviceToHost);
    checkCudaError(e);
    ck_mem_transfer -> clockPause();


    /* free device memory */
    cudaFree(d_res_copy);
    cudaFree(d_res);
    cudaFree(d_amp);

    ck_total -> clockPause();
}


int main(void){
    /* set basic data */
    create_sample_points(sample_points);
    generate_sig(sample_points, sig);
    // generate_sig_TEST(sig);

    /* create clocks */
    mkClockMeasure *ckCpu_fft_recur = new mkClockMeasure("CPU - FFT RECURSIVE");
    mkClockMeasure *ckCpu_fft_iter = new mkClockMeasure("CPU - FFT ITERATIVE");
    mkClockMeasure *ckGpu_mem_transfer = new mkClockMeasure("MEMORY TRANSFER");
    mkClockMeasure *ckGpu_exec = new mkClockMeasure("KERNELS");
    mkClockMeasure *ckGpu_fft = new mkClockMeasure("GPU - FFT (TOTAL)");
    ckCpu_fft_recur->clockReset();
    ckCpu_fft_iter->clockReset();
    ckGpu_mem_transfer->clockReset();
    ckGpu_exec->clockReset();
    ckGpu_fft->clockReset();


    for(int i=0; i<MAX_ITER; i++){
        /* CPU - FFT recurisve */
        cpu_cal_fft_recursive(sample_points, sig, ckCpu_fft_recur);

        /* CPU - FFT iterative */
        cpu_cal_fft_iterative(sample_points, sig, ckCpu_fft_iter);

        /* GPU - FFT */
        gpu_cal_fft(sample_points, sig, ckGpu_mem_transfer, ckGpu_exec, ckGpu_fft);
    }


    if(compareResult(fft_amp_rec, gpu_amp, N) && compareResult(fft_amp_iter, fft_amp_rec, N)){
        printf("SAMPLING_RATE : %d\nITERATION : %d\n\n", SAMPLING_RATE, MAX_ITER);
        printf("-------------------[CPU] FFT - recursive ---------------------\n");
        ckCpu_fft_recur->clockPrint();
        printf("\n-------------------[CPU] FFT - iterative---------------------\n");
        ckCpu_fft_iter->clockPrint();
        printf("\n-------------------[GPU] FFT --------------------------------\n");
        ckGpu_mem_transfer->clockPrint();
        ckGpu_exec->clockPrint();
        ckGpu_fft->clockPrint();


    }
    else{
        printf("Error: the two are not the same\n");
    }
    // save_data(sample_points, sig,  "data/original_signal.txt");
    save_data(freq, fft_amp_rec, "data/fft_freq_rec.txt");
    save_data(freq, fft_amp_iter, "data/fft_freq_iter.txt");
    save_data(freq, gpu_amp, "data/gpu_fft_freq.txt");
    save_data(sample_points, sig, "data/original_signal.txt");
}