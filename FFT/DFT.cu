/*
 * =====================================================================================
 *
 *       Filename:  DFT.cpp
 *
 *    Description:  Ch03 Samples
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
#include "mkClockMeasure.h"

const int SAMPLING_RATE = 100;
const int N = SAMPLING_RATE;
const int FREQ_NUM = 3;

double sample_points[N];
double freq[N];
double sig[N];
double idft_sig[N];
double xr[N];
double xi[N];
double x[N];
double freq_amp_ph[FREQ_NUM][3] = {{1, 3, 0}, {4, 1, 0}, {7, 0.5, 0}};

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

void cal_dft(double* sample_points, double* sig, mkClockMeasure* ck){
    ck->clockResume();

    for (int k = 0; k<N; k++){
        xi[k] = 0;
        xr[k] = 0;
        for (int n = 0; n < N; n++)
        {
            double bn = 2 * M_PI * k * n / N;
            xi[k] -= sig[n] * sin(bn);
            xr[k] += sig[n] * cos(bn);
        }
        x[k] = 2 * sqrt(pow(xr[k], 2) + pow(xi[k], 2)) / N;
        // printf("frequency: %f\tamplitude: %f\n", x[k], freq[k]);
    }
    ck->clockPause();
    printf("-------------------[CPU] DFT ---------------------\n");
    ck->clockPrint();
}

void cal_idft(mkClockMeasure *ck){
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
        // printf("t: %f\tamplitude: %f\n", n, idft_sig[n]);
    }
    ck->clockPause();
    printf("-------------------[CPU] Inverse DFT ---------------------\n");
    ck->clockPrint();
}

int main(void){
    mkClockMeasure *ckCpu_dft = new mkClockMeasure("CPU - DFT CODE");
    ckCpu_dft->clockReset();
    mkClockMeasure *ckCpu_idft = new mkClockMeasure("CPU - IDFT CODE");
    ckCpu_idft->clockReset();

    create_sample_points(sample_points);
    generate_sig(sample_points, sig);
    save_data(sample_points, sig,  "original_signal.txt");

    cal_dft(sample_points, sig, ckCpu_dft);
    save_data(freq, x, "dft_frequencies.txt");

    cal_idft(ckCpu_idft);
    save_data(sample_points, idft_sig, "idft_signal.txt");
}