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
#include "comp.h"
#include "mkClockMeasure.h"

const int SAMPLING_RATE = 4086;
const int N = SAMPLING_RATE;
const int FREQ_NUM = 3;
const int MAX_ITER = 1;

double freq_amp_ph[FREQ_NUM][3] = {{1, 3, 0}, {4, 1, 0}, {7, 0.5, 0}};
double sample_points[N], freq[N], sig[N];
Comp sig_comp[N];
__constant__ double d_sig[N];
//cpu
double fft_amp_rec[N], fft_amp_iter[N];
Comp fft_res_rec[N], fft_res_iter[N];
//gpu
double gpu_amp[N];
Comp gpu_x[N];


//data settings
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
    db_to_comp(sig, sig_comp);
}

void generate_sig_TEST(double* sig){
    for(int i=0; i< N; i++){
        sig[i] = i;
    }
}

void save_data(double* x_values, double* y_values, const char* path){
    FILE *dataf = fopen(path, "w");
    for (int i=0; i<N; i++){
    fprintf(dataf, "%lf %lf\n", x_values[i], y_values[i]);
    }
}

bool compareResult(double* a, double* b, int size){
    double epsilon = 0.000001;

    for(int i=0; i<size; i++){
        if(fabs(a[i] - b[i]) < epsilon){
            // printf("a; %lf, \t b: %lf\n", a[i], b[i]);
            return true;
        }
    }
    return true;
}

void db_to_comp(double* a, Comp* b, int size){
    for(int i=0; i<size; i++){
        b[i].r = a[i];
        b[i].i = 0;
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
    fft_recursive(half_len, even);
    fft_recursive(half_len, odd);

    //conquer
    for(int k=0; k<half_len; k++){

        Comp w_k = cal_euler(-2 * M_PI * k / len);
        Comp t = comp_mult(w_k, odd[k]);

        x[k] = comp_add(even[k], t);
        x[k + half_len] = comp_sub(even[k], t);
    }
}

void cpu_cal_fft_recursive(double* sample_points, double* sig, mkClockMeasure* ck){
    for(int j=0; j < MAX_ITER; j++){
        //initialize fft as signal values
        db_to_comp(sig, fft_res_rec);
        
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
}

void cpu_fft_iterative(int len, Comp* x){
    int depth = (int)log2(len);
    int half_len = len/2;
    Comp x_copy[N];

    for(int l=2, d=1; l<=len; l*=2, d*=2){
        int itvl = len / l;
        for(int i=0; i<N; i++){
            x_copy[i].r = x[i].r;
            x_copy[i].i = x[i].i;
        }

        for(int i=0; i<itvl; i++){
            for(int k=0; k<d; k++){
                Comp w_k = cal_euler(-2 * M_PI * k / l);
                Comp t = comp_mult(w_k, x_copy[i + k*itvl*2 + itvl]);
                Comp even = x_copy[i + k*itvl*2];
                x[i + k*itvl] = comp_add(even, t);
                x[i + k*itvl + half_len] = comp_sub(even, t);
            }
        }
    }
}

void cpu_cal_fft_iterative(double* sample_points, double* sig, mkClockMeasure* ck){
    for(int j=0; j < MAX_ITER; j++){
        //initialize fft as signal values
        db_to_comp(sig, fft_res_iter);
        
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
}

__global__ void gpu_fft(Comp* d_x, Comp* d_x_copy, int depth, int len, int half_len, int N){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i < N){
        comp_cpy(d_x[i], d_x_copy[i]);

        for()
    }
}


int main(void){
    /* set basic data */
    create_sample_points(sample_points);
    generate_sig(sample_points, sig);
    // generate_sig_TEST(sig);

    /* create clocks */
    mkClockMeasure *ckCpu_fft_recur = new mkClockMeasure("CPU - FFT RECURSIVE");
    mkClockMeasure *ckCpu_fft_iter = new mkClockMeasure("CPU - FFT ITERATIVE");
    mkClockMeasure *ckGpu_fft = new mkClockMeasure("GPU - FFT");
    ckCpu_fft_recur->clockReset();
    ckCpu_fft_iter->clockReset();
    ckGpu_fft->clockReset();

    /* allocate device memory */ 
    Comp *d_x, *d_x_copy;
    double *d_amp;
    int db_bytesize = N * sizeof(double);
    int comp_bytesize = N * 2 * sizeof(double);

    cudaError_t err = cudaMalloc((void**)&d_x, comp_bytesize);
    checkCudaError(err);
    err = cudaMalloc((void**)&d_xd_x_copy, comp_bytesize);
    checkCudaError(err);
    err = cudaMalloc((void**)&d_amp, db_bytesize);
    checkCudaError(err);

    /* set thread, grid size */ 
    int thread = 256;
    int tbSize = N/thread;
    dim3 gridSize(tbSize, 1, 1);
    dim3 blockSize(thread, 1, 1);


for(int i=0; i<MAX_ITER; i++){
        /* GPU - FFT */
        ckGpu_fft->clockResume();

        // memory transfer: host -> device
        err = cudaMemcpy(d_x, sig_comp, comp_bytesize);
        checkCudaError(err);
        err = cudaMemcpy(d_x_copy, sig_comp, comp_bytesize);
        checkCudaError(err);

        gpu_dft<<<gridSize, blockSize>>>(d_x, d_amp, N); 

        // memory transfer: device -> host 
        err = cudaMemcpy(gpu_x, d_x, comp_bytesize, cudaMemcpyDeviceToHost);
        checkCudaError(err);
        err = cudaMemcpy(gpu_amp, d_amp, db_bytesize, cudaMemcpyDeviceToHost);
        checkCudaError(err);

        ckGpu_fft->clockPause();
    }

    /* free device memory */
    cudaFree(d_x_copy);
    cudaFree(d_x);
    cudaFree(d_amp);

    /* fft recursive version */
    cpu_cal_fft_recursive(sample_points, sig, ckCpu_fft_recur);

    /* fft iterative version */
    cpu_cal_fft_iterative(sample_points, sig, ckCpu_fft_iter);



    if(compareResult(fft_amp_rec, gpu_amp, N)){
        printf("SAMPLING_RATE : %d\nITERATION : %d\n\n", SAMPLING_RATE, MAX_ITER);
        printf("-------------------[CPU] FFT - recursive ---------------------\n");
        ckCpu_fft_recur->clockPrint();
        printf("\n-------------------[CPU] FFT - iterative---------------------\n");
        ckCpu_fft_iter->clockPrint();
        printf("\n-------------------[GPU] FFT --------------------------------\n");
        ckGpu_fft->clockPrint();

    }
    else{
        printf("Error: the two are not the same\n");
    }
    // save_data(sample_points, sig,  "data/original_signal.txt");
    save_data(freq, fft_amp_rec, "data/fft_freq_rec.txt");
    save_data(freq, fft_amp_iter, "data/fft_freq_iter.txt");
    save_data(freq, gpu_amp, "data/gpu_fft_freq_iter.txt");
    save_data(sample_points, sig, "data/original_signal.txt");
}