#include <stdlib.h>
#include <stdio.h>
#include "fft.h"

#define PI 3.14159265358979323846

struct vec {
    uint64_t length; 
    double *data;
};

struct vec_z{
    uint64_t length;
    struct ZNum *data;
};

static void _FFT_radix_2(struct vec_z *samples, struct vec_z *coef, struct vec_z *copy, uint64_t start, uint64_t num, uint64_t space, int8_t direction) {
    if (num == 1) {
        coef->data[start]= samples->data[start];
        return;
    }
    /* Split into Even/Odd */
    _FFT_radix_2(samples, copy, coef, start, num / 2, space * 2, direction);
    _FFT_radix_2(samples, copy, coef, start + space, num / 2, space * 2, direction);
    struct ZNum twiddle, tmp, result;
    uint64_t i, current, even, odd, half = num * space / 2, double_space = space * 2;
    for (i = 0; i < num / 2; i++) {
        current = i * space + start;
        even = i * double_space + start;
        odd = even + space;
        ZNum_root_of_unity(direction * 2.0 * i * PI / num, &twiddle);
        tmp = copy->data[even];
        /* i: even + twiddle * odd */ 
        ZNum_mul(&twiddle, &copy->data[odd], &result);
        ZNum_add(&tmp, &result, &coef->data[current]);
        /* i + N/2: even - twiddle * odd */
        ZNum_sub(&tmp, &result, &coef->data[current + half]);
    }
}

static void _FFT_real_radix_2(struct vec *samples, struct vec_z *coef, struct vec_z *copy, uint64_t start, uint64_t num, uint64_t space) {
    if (num == 1) {
        coef->data[start].real = samples->data[start];
        coef->data[start].imag = 0.0;
        return;
    }
    /* Split into Even/Odd */
    _FFT_real_radix_2(samples, copy, coef, start, num / 2, space * 2);
    _FFT_real_radix_2(samples, copy, coef, start + space, num / 2, space * 2);
    struct ZNum twiddle, tmp, result;
    uint64_t i, current, even, odd, half = num * space / 2, double_space = space * 2;
    for (i = 0; i < num / 2; i++) {
        current = i * space + start;
        even = i * double_space + start;
        odd = even + space;
        ZNum_root_of_unity(-2.0 * i * PI / num, &twiddle);
        tmp = copy->data[even];
        /* i: even + twiddle * odd */ 
        ZNum_mul(&twiddle, &copy->data[odd], &result);
        ZNum_add(&tmp, &result, &coef->data[current]);
        /* i + N/2: even - twiddle * odd */
        ZNum_sub(&tmp, &result, &coef->data[current + half]);
    }
}

/* Pre: samples and coeff are of length 2^n */
void FFT_radix_2(struct ZNum *samples, struct ZNum *coef, uint64_t length, int8_t direction) {
    struct vec_z samples_vec = {length, samples};
    struct vec_z coef_vec = {length, coef};
    struct ZNum *copy = malloc(sizeof(struct ZNum) * length);
    struct vec_z copy_vec = {length, copy};
    if (direction == FFT_BACKWARD) { /* Need to scale by 1/length */
        _FFT_radix_2(&coef_vec, &samples_vec, &copy_vec, 0, length, 1, FFT_BACKWARD); 
        uint64_t i;
        for (i = 0; i < length; i++) {
            samples[i].real /= (double)length;
            samples[i].imag /= (double)length;
        }
    } else {
        _FFT_radix_2(&samples_vec, &coef_vec, &copy_vec, 0, length, 1, FFT_FORWARD); 
    }
    free(copy);
}

/* Pre: samples and coeff are of length 2^n */
void FFT_real_radix_2(double *samples, struct ZNum *coef, uint64_t length, int8_t direction) {
    struct vec_z coef_vec = {length, coef};
    struct ZNum *copy = malloc(sizeof(struct ZNum) * length);
    struct vec_z copy_vec = {length, copy};
    if (direction == FFT_BACKWARD) {
        struct ZNum *to_complex = malloc(sizeof(struct ZNum) * length);
        struct vec_z samples_vec = {length, to_complex};
        _FFT_radix_2(&coef_vec, &samples_vec, &copy_vec, 0, length, 1, FFT_BACKWARD);
        uint64_t i;
        for (i = 0; i < length; i++)
            samples[i] = to_complex[i].real / (double)length;
        free(to_complex);
    } else {
        struct vec samples_vec = {length, samples};
        _FFT_real_radix_2(&samples_vec, &coef_vec, &copy_vec, 0, length, 1);
    }
    free(copy);
}
