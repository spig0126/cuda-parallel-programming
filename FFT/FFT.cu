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

const int SAMPLING_RATE = 32768;
const int N = SAMPLING_RATE;
const int FREQ_NUM = 3;

double freq_amp_ph[FREQ_NUM][3] = {{1, 3, 0}, {4, 1, 0}, {7, 0.5, 0}};
double freq[N];
double fft_amp[N];
double sample_points[N];
double sig[N];
Comp fft_res[N];

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

void save_data(double* x_values, double* y_values, const char* path){
    FILE *dataf = fopen(path, "w");
    for (int i=0; i<N; i++){
    fprintf(dataf, "%lf %lf\n", x_values[i], y_values[i]);
    }
}

//radix-2 cooley-tukey fft
void fft_recursive(int len, Comp* x){ //x = signal 값
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

void fft_iterative(int len, Comp* x){
    int depth = (int)log2(len);
    int half_len == len >> 1;

    for(int l == len; l > 0; l/=2){
        int itvl = l >> 1;

        for(int k=0; k<half_len; k++){
            Comp w_k = cal_euler(-2 * M_PI * k / l);
            Comp odd = comp_mult(w_k, x[k + itvl]);
            Comp even = x[k];
            x[k] = comp_add(even, odd);
            x[k + half_len] = comp_sub(even, odd);
        }
    }
}

void cal_fft(double* sample_points, double* sig, mkClockMeasure* ck){
    //initialize fft as signal values
    for(int i=0; i<N; i++){
        fft_res[i].r = sig[i];
        fft_res[i].i = 0;
    }
    
    //start clock
    ck -> clockResume();

    //calculate fft
    fft_recursive(N, fft_res);

    //cal magnitude of each frequency
    for(int k=0; k<N; k++){
        fft_amp[k] = 2 * sqrt(pow(fft_res[k].r, 2) + pow(fft_res[k].i, 2)) / N;
    }

    //print clock record
    ck->clockPause();
    printf("SAMPLING_RATE : %d\n", SAMPLING_RATE);
    printf("-------------------[CPU] FFT ---------------------\n");
    ck->clockPrint();
}




int main(void){
    mkClockMeasure *ckCpu_dft = new mkClockMeasure("CPU - FFT CODE");
    ckCpu_dft->clockReset();

    create_sample_points(sample_points);
    generate_sig(sample_points, sig);
    // cal_fft(sample_points, sig, ckCpu_dft);

    /* fft iterative version */
    //initialize fft as signal values
    for(int i=0; i<N; i++){
        fft_res[i].r = sig[i];
        fft_res[i].i = 0;
    }
    fft_iterative(N, fft_res);

    save_data(sample_points, sig,  "data/original_signal.txt");
    // save_data(freq, fft_amp, "fft_frequencies.txt");
}