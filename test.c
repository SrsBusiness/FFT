#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fft.h"

#define PI 3.14159265358979323846
#define NUM 8

int main() {
    printf("============ Fourier Tests ============\n");
    //double samples[NUM];
    //putchar('[');
    //for (i = 0; i < NUM; i++) {
    //    samples[i] = (double)(rand() % 100) + ((double)rand()) / ((double)rand());
    //    printf("%f, ", samples[i]);
    //}
    //printf("]\n");

    double samples[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    struct ZNum samples_z[] = {{1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0},
        {5.0, 0.0}, {6.0, 0.0}, {7.0, 0.0}, {8.0, 0.0}};
    int i;
    struct ZNum coef[NUM];

    printf("Fourier Transform\n");
    FFT_real_radix_2(samples, coef, NUM, FFT_FORWARD);
    for (i = 0; i < NUM; i++) {
        ZNum_print(&coef[i]);
    }
    
    printf("\n");
    memset(samples, 0, sizeof(samples));
    printf("Inverse Fourier Transform\n");
    FFT_real_radix_2(samples, coef, NUM, FFT_BACKWARD);
    for (i = 0; i < NUM; i++) {
        printf("%f\n", samples[i]);
    }
    
    printf("\n");
    printf("Fourier Transform\n");
    FFT_radix_2(samples_z, coef, NUM, FFT_FORWARD);
    for (i = 0; i < NUM; i++) {
        ZNum_print(&coef[i]);
    }
    
    printf("\n");
    printf("Inverse Fourier Transform\n");
    FFT_radix_2(samples_z, coef, NUM, FFT_BACKWARD);
    for (i = 0; i < NUM; i++) {
        ZNum_print(&samples_z[i]);
    }

    return 0;
}
