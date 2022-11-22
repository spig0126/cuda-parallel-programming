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
#include "comp.cuh"
#include "mkCuda.h"
#include "mkClockMeasure.h"
#include "data.cuh"

const unsigned long N = 131072;
const unsigned long SAMPLING_RATE = N;
const unsigned long BLOCK_SIZE = 32;
const unsigned long GRID_SIZE = N/BLOCK_SIZE;
const unsigned long BATCH = 1; //batch는 동시에 수행하는 fourier transform 개수!!
const unsigned long NX = N;

const unsigned long FREQ_NUM = 3;
const unsigned long MAX_ITER = 1;

fp freq_amp_ph[FREQ_NUM][3] = {{1, 3, 0}, {4, 1, 0}, {7, 0.5, 0}};
fp sample_points[N], freq[N], sig[N];
/* cpu */
fp fft_amp_rec[N], fft_amp_iter[N];
fpComplex fft_res_rec[N], fft_res_iter[N];
/* gpu */
fp gpu_amp[N];
fpComplex gpu_res[N];
/* cuFFT */
fp cufft_amp[N];
int db_bytesize = N * sizeof(fp);
int comp_bytesize = N * sizeof(fpComplex);

/* thread, grid size */ 
dim3 half_gridSize(GRID_SIZE/2, 1, 1);
dim3 total_gridSize(GRID_SIZE, 1, 1);
dim3 blockSize(BLOCK_SIZE, 1, 1);

/* basic functions */
void db_to_comp(fp* a, fpComplex* b, int size){
    for(int i=0; i<size; i++){
        b[i].x = a[i];
        b[i].y = 0;
    }
}
bool compareResult(fp* a, fp* b, int size){
    fp epsilon = 0.000001;

    for(int i=0; i<size; i++){
        if(fabs(a[i] - b[i]) > epsilon){
            // printf("a; %lf, \t b: %lf\n", a[i], b[i]);
            return false;
        }
    }
    return true;
}
void save_data(fp* x_values, fp* y_values, const char* path){
    FILE *dataf = fopen(path, "w");
    for (int i=0; i<N; i++){
    fprintf(dataf, "%lf %lf\n", x_values[i], y_values[i]);
    }
}

/* data settings */
void create_sample_points(fp* sample_points){
    for (int i = 0; i<N; i++){
        sample_points[i] = (fp)i/SAMPLING_RATE;
        freq[i] = i;
    }
}
void generate_sig(fp* sample_points, fp* sig){
    for (int s_i = 0; s_i < N; s_i++){
        sig[s_i] = 0;
        for (int f_i = 0; f_i < FREQ_NUM; f_i++)
        {
            sig[s_i] += freq_amp_ph[f_i][1] * sin(2 * M_PI * freq_amp_ph[f_i][0] * sample_points[s_i]);
        }
    }
}
void generate_sig_TEST(fp* sig){
    for(int i=0; i< N; i++){
        sig[i] = i;
    }
}

//radix-2 cooley-tukey fft
void cpu_fft_recursive(int len, fpComplex* x){ //x = signal 값
    if(len == 1){
        return;
    }
    int half_len = len >> 1;

    //divide x into 2 subgroups: even, odd
    fpComplex even[half_len], odd[half_len];
    for(int i=0; i<half_len; i++){
        even[i] = x[2 * i];
        odd[i] = x[2 * i + 1];
    }

    //divde
    cpu_fft_recursive(half_len, even);
    cpu_fft_recursive(half_len, odd);

    //conquer & combine
    for(int k=0; k<half_len; k++){

        fpComplex w_k = cal_euler(-2 * M_PI * k / len);
        fpComplex t = cuCmulf(w_k, odd[k]);

        x[k] = cuCaddf(even[k], t);
        x[k + half_len] = cuCsubf(even[k], t);
        
        // if(len == 16){
        //     printf("\tfft even : ");
        //      for(int i=0; i<half_len; i++){
        //          printf("%d   ", int(even[i].x));
        //      }
        //      printf("\n");
        //      printf("\tfft odd : ");
        //      for(int i=0; i<half_len; i++){
        //          printf("%d   ", int(odd[i].x));
        //      }
        //      printf("\n");
        //     printf("\t\t\tx[k]: %lf, x[k + half_len]: %lf \n", x[k].x,  x[k+ half_len].x);
        //     printf("\t\t\teven: %lf, odd: %lf \n", even[k].x, odd[k].x);
        // }
    }
}

void cpu_cal_fft_recursive(fp* sample_points, fp* sig, mkClockMeasure* ck){
    //initialize fft as signal values
    db_to_comp(sig, fft_res_rec, N);
    
    //start clock
    ck -> clockResume();

    //calculate fft
    cpu_fft_recursive(N, fft_res_rec);

    //cal magnitude of each frequency
    for(int k=0; k<N; k++){
        fft_amp_rec[k] = 2 * sqrt(pow(cuCrealf(fft_res_rec[k]), 2) + pow(cuCimagf(fft_res_rec[k]), 2)) / N;
    }

    ck->clockPause();
}

void cpu_fft_iterative(int len, fpComplex* x){
    int depth = (int)log2(len);
    int half_len = len/2;
    fpComplex x_copy[N];

    // 1 iteration = 1 layer
    for(int l=2, d=1; l<=len; l<<=1, d<<=1){
        // update data from lower layer
        // printf("len: %d\n", l);
        for(int i=0; i<N; i++){
            x_copy[i].x = cuCrealf(x[i]);
            x_copy[i].y = cuCimagf(x[i]);
        }

        int itvl = len / l;

        for(int i=0; i<itvl; i++){
            for(int k=0; k<d; k++){
                int temp = k * itvl;
                int data_even_i = i + temp*2;
                int data_odd_i = data_even_i + itvl;
                int res_i = i + temp;
                int res_half_i = res_i + half_len;

                fpComplex w_k = cal_euler(-2 * M_PI * k / l);
                fpComplex t = cuCmulf(w_k, x_copy[data_odd_i]);
                fpComplex even = x_copy[data_even_i];
                x[res_i] = cuCaddf(even, t);
                x[res_half_i] = cuCsubf(even, t);

                // if(l == 4){
                //      printf("\t\tx[k]: %lf, x[k + half_len]: %lf\n",  cuCrealf(x[i + k*itvl]), cuCrealf(x[i + k*itvl + half_len]));
                //      printf("\t\teven: %lf, odd: %lf \n", cuCrealf(x_copy[i + k*itvl*2]), cuCrealf(x_copy[i + k*itvl*2 + itvl]));
                //  }
            }
        }
    }
}

void cpu_cal_fft_iterative(fp* sample_points, fp* sig, mkClockMeasure* ck){
    //initialize fft as signal values
    db_to_comp(sig, fft_res_iter, N);
    
    //start clock
    ck -> clockResume();

    //calculate fft
    cpu_fft_iterative(N, fft_res_iter);

    //cal magnitude of each frequency
    for(int k=0; k<N; k++){
        fft_amp_iter[k] = 2 * sqrt(pow(cuCrealf(fft_res_iter[k]), 2) + pow(cuCimagf(fft_res_iter[k]), 2)) / N;
    }

    ck->clockPause();
}

__global__ void gpu_cal_amp(fpComplex* d_res, fp* d_amp){
    int k = blockDim.x * blockIdx.x + threadIdx.x;

    d_amp[k] = 2 * sqrt(pow(cuCrealf(d_res[k]), 2) + pow(cuCimagf(d_res[k]), 2)) / N;
}

__global__ void gpu_update_data(fpComplex* d_res, fpComplex* d_res_copy){
    /* update data from lower layer */
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    d_res_copy[i].x = d_res[i].x;
    d_res_copy[i].y = d_res[i].y;
}

__global__ void gpu_fft(fpComplex* after, fpComplex* before, int len, int half_len, int itvl){
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    int k = (int)(i / itvl);
    int data_even_i = i % itvl + k * itvl * 2;
    int data_odd_i = data_even_i + itvl;
    int res_i = i;
    int res_half_i = i + half_len;
  
    fpComplex w_k = cal_euler(-2 * M_PI * k / len);
    fpComplex t = cuCmulf(w_k, before[data_odd_i]);
    fpComplex even = before[data_even_i];
    after[res_i] = cuCaddf(even, t);
    after[res_half_i] = cuCsubf(even, t);

    // if(len == 4){
    //     printf("k: %d \n", k);
    //     printf("\t\ti: %d, k: %d\t=>\tx[k]: %lf, x[k + half_len]: %lf\n",  i, k, after[res_i].x, after[res_half_i].x);
    //     printf("\t\ti: %d, k: %d\t=>\teven: %lf, odd: %lf \n", i, k, before[data_even_i].x, before[data_odd_i].x);
    // }
}

void gpu_cal_fft(fp* sample_points, fp* sig, mkClockMeasure* ck_mem_transfer, mkClockMeasure* ck_kernels, mkClockMeasure* ck_total){
    /* initialize fft as signal values */
    db_to_comp(sig, gpu_res, N);

    /* allocate device memory */ 
    fpComplex *d_res, *d_res_copy;
    fp *d_amp;

    cudaError_t e = cudaMalloc((void**)&d_res, comp_bytesize);
    checkCudaError(e);
    e = cudaMalloc((void**)&d_res_copy, comp_bytesize);
    checkCudaError(e);
    e = cudaMalloc((void**)&d_amp, db_bytesize);
    checkCudaError(e);

    ck_total -> clockResume();

    /* data transfer - host -> dev */
    ck_mem_transfer -> clockResume();
    e = cudaMemcpy(d_res_copy, gpu_res, comp_bytesize, cudaMemcpyHostToDevice);
    checkCudaError(e);
    ck_mem_transfer -> clockPause();

    /* cal Fourier efficients */
    ck_kernels -> clockResume();
    int depth = (int)log2(N);
    int half_len = N / 2;
    for(int l=2, d=1; l <= N; l<<=1, d+=1){ 
        int itvl = N / l;
        fpComplex *after = d%2== 1? d_res: d_res_copy;
        fpComplex *before = d%2 == 1? d_res_copy: d_res;
        
        gpu_fft<<<half_gridSize, blockSize>>>(after, before, l, half_len, itvl);
        e=cudaDeviceSynchronize();

	    checkCudaError(e);
    }
    ck_kernels -> clockPause();

    gpu_update_data<<<total_gridSize, blockSize>>>(d_res, d_res_copy);
    e=cudaDeviceSynchronize();
    gpu_cal_amp<<<total_gridSize, blockSize>>>(d_res, d_amp);

    /* data transfer - dev -> host */
    ck_mem_transfer -> clockResume();
    e = cudaMemcpy(gpu_amp, d_amp, db_bytesize, cudaMemcpyDeviceToHost);
    checkCudaError(e);
    ck_mem_transfer -> clockPause();

    ck_total -> clockPause();

    /* free device memory */
    cudaFree(d_res_copy);
    cudaFree(d_res);
    cudaFree(d_amp);
}

void cal_using_cufft(fp* signal, mkClockMeasure* ck_total, mkClockMeasure* ck_plan, mkClockMeasure* ck_execute){
    /* allocate device memory */ 
    //amplitude related
    fp *d_amp;
    cudaMalloc((void**)&d_amp, db_bytesize);
    //cufft related
    cufftHandle plan;
    fpComplex *data;
    cudaMalloc((void**)&data, comp_bytesize);
    fpComplex sig_comp[N];
    db_to_comp(signal, sig_comp, N);

    //copy signal values to device(cufft needs to have all the data in the device)
    cudaMemcpy(data, sig_comp, comp_bytesize, cudaMemcpyHostToDevice);

    /* create 1D FFT plan */
    ck_total->clockResume();
    ck_plan->clockResume();
    if(cufftPlan1d(&plan, NX, CUFFT_C2C, BATCH) != CUFFT_SUCCESS){
        printf("CUFFT error: paln creation failed\n");
    }
    ck_plan->clockPause();

    /* execute cuFFT plan for transformation */
    ck_execute->clockResume();
    if (cufftExecC2C(plan, data, data, CUFFT_FORWARD) != CUFFT_SUCCESS){
        printf("CUFFT error: ExecD2C forward failed\n");
    }
    ck_execute->clockPause();
    ck_total->clockPause();


    gpu_cal_amp<<<total_gridSize, blockSize>>>(data, d_amp);
    cudaMemcpy(cufft_amp, d_amp, db_bytesize, cudaMemcpyDeviceToHost);

    if(cudaDeviceSynchronize() != cudaSuccess){
        printf("Cuda error: failed to synchronize\n");
    }

    /* destroy plan */
    cudaFree(data);
    cudaFree(d_amp);
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
    mkClockMeasure *ck_cufft = new mkClockMeasure("cuFFT");
    mkClockMeasure *ck_cufft_plan = new mkClockMeasure("cuFFT - PLAN");
    mkClockMeasure *ck_cufft_exec = new mkClockMeasure("cuFFT - EXECUTE");
    ckCpu_fft_recur->clockReset();
    ckCpu_fft_iter->clockReset();
    ckGpu_mem_transfer->clockReset();
    ckGpu_exec->clockReset();
    ckGpu_fft->clockReset();
    ck_cufft->clockReset();
    ck_cufft_plan->clockReset();
    ck_cufft_exec->clockReset();


    for(int i=0; i<MAX_ITER; i++){
        /* CPU - FFT recurisve */
        // cpu_cal_fft_recursive(sample_points, sig, ckCpu_fft_recur);

        /* CPU - FFT iterative */
        cpu_cal_fft_iterative(sample_points, sig, ckCpu_fft_iter);

        /* GPU - FFT */
        gpu_cal_fft(sample_points, sig, ckGpu_mem_transfer, ckGpu_exec, ckGpu_fft);

        /* cufft */
        cal_using_cufft(sig, ck_cufft, ck_cufft_plan, ck_cufft_exec);
    }


    if(compareResult(fft_amp_iter, gpu_amp, N) && compareResult(fft_amp_iter, cufft_amp, N)){ 
        printf("SAMPLING_RATE : %ld\nITERATION : %ld\n\n", SAMPLING_RATE, MAX_ITER);
        printf("-------------------[CPU] FFT - recursive ---------------------\n");
        ckCpu_fft_recur->clockPrint();
        printf("\n-------------------[CPU] FFT - iterative---------------------\n");
        ckCpu_fft_iter->clockPrint();
        printf("\n-------------------[GPU] FFT --------------------------------\n");
        ckGpu_mem_transfer->clockPrint();
        ckGpu_exec->clockPrint();
        ckGpu_fft->clockPrint();
        printf("\n------------------- CUFFT --------------------------------\n");
        ck_cufft->clockPrint();
        ck_cufft_plan->clockPrint();
        ck_cufft_exec->clockPrint();

    }
    else{
        printf("Error: the two are not the same\n");
    }
    
    // save_data(sample_points, sig,  "data/original_signal.txt");
    save_data(freq, fft_amp_rec, "data/fft_freq_rec.txt");
    save_data(freq, fft_amp_iter, "data/fft_freq_iter.txt");
    save_data(freq, gpu_amp, "data/gpu_fft_freq.txt");
    save_data(freq, cufft_amp, "data/cufft_freq.txt");
    save_data(sample_points, sig, "data/original_signal.txt");
}