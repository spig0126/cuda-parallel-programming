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
#include <complex.h>
#include "mkClockMeasure.h"

const int SAMPLING_RATE = 32768;
const int N = SAMPLING_RATE;
const int FREQ_NUM = 3;

typedef struct Comp{
    double r, i;
} Comp;

double freq_amp_ph[FREQ_NUM][3] = {{1, 3, 0}, {4, 1, 0}, {7, 0.5, 0}};
double freq[N];
double fft_amp[N];
double sample_points[N];
double sig[N];
Comp fft_res[N];

//Comp operation functions
Comp cal_euler(double x){
    Comp res;
    res.r = cos(x);
    res.i = sin(x);
    return res;
}

Comp comp_mult(Comp a, Comp b){
    Comp res;
    res.r = a.r * b.r - a.i * b.i;
    res.i = a.r * b.i + a.i * b.r;
    return res;
}

Comp comp_add(Comp a, Comp b){
    Comp res;   
    res.r = a.r + b.r;
    res.i = a.i + b.i;
    return res;
}

Comp comp_sub(Comp a, Comp b){
    Comp res;   
    res.r = a.r - b.r;
    res.i = a.i - b.i;
    return res;
}

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

//algorithms
void fft(int len, Comp* x){ //x = signal 값
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
    fft(half_len, even);
    fft(half_len, odd);

    //conquer
    for(int k=0; k<half_len; k++){
        Comp w_k = cal_euler(-2 * M_PI * k / len);
        Comp t = comp_mult(w_k, odd[k]);
        x[k] = comp_add(even[k], t);
        x[k + half_len] = comp_sub(even[k], t);
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
    fft(N, fft_res);

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
    cal_fft(sample_points, sig, ckCpu_dft);

    save_data(sample_points, sig,  "original_signal.txt");
    save_data(freq, fft_amp, "fft_frequencies.txt");
}