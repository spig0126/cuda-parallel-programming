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

const int SAMPLING_RATE = 256;
const int N = SAMPLING_RATE;
const int FREQ_NUM = 3;
const int MAX_ITER = 1;

double freq_amp_ph[FREQ_NUM][3] = {{1, 3, 0}, {4, 1, 0}, {7, 0.5, 0}};
double freq[N];
double fft_amp_rec[N], fft_amp_iter[N];
double sample_points[N];
double sig[N];
Comp fft_res_rec[N], fft_res_iter[N];

//Comp operation functions


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
    double epsilon = 0.0000001f;

    for(int i=0; i<size; i++){
        if(fabs(a[i] - b[i]) < epsilon){
            // printf("a; %lf, \t b: %lf\n", a[i], b[i]);
            return true;
        }
    }
    return true;
}

//radix-2 cooley-tukey fft
int cnt = 0;

void fft_recursive(int len, Comp* x){ //x = signal 값
    // printf("#%d\n", ++cnt);
    if(len == 1){
        // printf("\tlen: 1\n\t\tval: %d\n", int(x[0].r));
        return;
    }
    // printf("\tlen: %d\n", len);

    int half_len = len >> 1;

    //divide x into 2 subgroups: even, odd
    Comp even[half_len], odd[half_len];
    for(int i=0; i<half_len; i++){
        even[i] = x[2 * i];
        odd[i] = x[2 * i + 1];
    }

    //conquer
    fft_recursive(half_len, even);
    fft_recursive(half_len, odd);

    //combine
    for(int k=0; k<half_len; k++){
        Comp w_k = cal_euler(-2 * M_PI * k / len);
        Comp t = comp_mult(w_k, odd[k]);
        x[k] = comp_add(even[k], t);
        x[k + half_len] = comp_sub(even[k], t);

        // printf("\t\t\t\tsig[%d] = even %d * exp( %d * %d * %d)", int(x[k].r), N / len, );

        // printf("\t\t fft even : ");
        // for(int i=0; i<half_len; i++){
        //     printf("%d   ", int(even[i].r));
        // }
        // printf("\n");
        // printf("\t\tfft odd : ");
        // for(int i=0; i<half_len; i++){
        //     printf("%d   ", int(odd[i].r));
        // }
        // printf("\n");
  
        // printf("\t\t\tlen: %d, k: %d, x[k]: %lf, x[k + half_len]: %lf \n", len, k, x[k].r,  x[k+ half_len].r);

    
    }
}

void cal_fft_recursive(double* sample_points, double* sig, mkClockMeasure* ck){
    for(int j=0; j < MAX_ITER; j++){
        //initialize fft as signal values
        for(int i=0; i<N; i++){
            fft_res_rec[i].r = sig[i];
            fft_res_rec[i].i = 0;
        }
        
        //start clock
        ck -> clockResume();

        //calculate fft
        fft_recursive(N, fft_res_rec);

        //cal magnitude of each frequency
        for(int k=0; k<N; k++){
            fft_amp_rec[k] = 2 * sqrt(pow(fft_res_rec[k].r, 2) + pow(fft_res_rec[k].i, 2)) / N;
        }

        ck->clockPause();
    }
}

void fft_iterative(int len, Comp* x){
    int depth = (int)log2(len);
    int half_len = len>>1;

    for(int l=2; l<=len; l*=2){
        int itvl = len / l;
        for(int k=0; k<half_len; k++){
            for(int i=0; i<itvl; i++){
                Comp w_k = cal_euler(-2 * M_PI * k / l);
                Comp t = comp_mult(w_k, x[i + itvl]);
                x[k] = comp_add(x[i], t);
                x[k + itvl] = comp_sub(x[i], t);
            }
        }
    }
}

void cal_fft_iterative(double* sample_points, double* sig, mkClockMeasure* ck){
    for(int j=0; j < MAX_ITER; j++){
        //initialize fft as signal values
        for(int i=0; i<N; i++){
            fft_res_iter[i].r = sig[i];
            fft_res_iter[i].i = 0;
        }
        
        //start clock
        ck -> clockResume();

        //calculate fft
        fft_iterative(N, fft_res_iter);

        //cal magnitude of each frequency
        for(int k=0; k<N; k++){
            fft_amp_iter[k] = 2 * sqrt(pow(fft_res_iter[k].r, 2) + pow(fft_res_iter[k].i, 2)) / N;
        }

        ck->clockPause();
    }
}





int main(void){
    mkClockMeasure *ckCpu_fft_recur = new mkClockMeasure("CPU - FFT RECURSIVE");
    mkClockMeasure *ckCpu_fft_iter = new mkClockMeasure("CPU - FFT ITERATIVE");
    ckCpu_fft_recur->clockReset();
    ckCpu_fft_iter->clockReset();

    create_sample_points(sample_points);
    // generate_sig(sample_points, sig);
    generate_sig_TEST(sig);

    /* fft recursive version */
    cal_fft_recursive(sample_points, sig, ckCpu_fft_recur);

    /* fft iterative version */
    cal_fft_iterative(sample_points, sig, ckCpu_fft_iter);

    if(compareResult(fft_amp_rec, fft_amp_iter, N)){
        printf("SAMPLING_RATE : %d\nITERATION : %d\n\n", SAMPLING_RATE, MAX_ITER);
        printf("-------------------[CPU] FFT - recursive ---------------------\n");
        ckCpu_fft_recur->clockPrint();
        printf("\n-------------------[CPU] FFT - iterative---------------------\n");
        ckCpu_fft_iter->clockPrint();
    }
    else{
        printf("Error: the two are not the same\n");
    }
    // save_data(sample_points, sig,  "data/original_signal.txt");
    save_data(freq, fft_amp_rec, "data/fft_freq_rec.txt");
    save_data(freq, fft_amp_iter, "data/fft_freq_iter.txt");
    save_data(sample_points, sig, "data/original_signal.txt");
}