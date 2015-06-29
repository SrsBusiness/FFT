#include <math.h>
#include <stdio.h>
#include "complex.h"

#define PI 3.14159265358979323846

void ZNum_add(struct ZNum *a, struct ZNum *b, struct ZNum *result) {
    result->real = a->real + b->real;
    result->imag = a->imag + b->imag;
}

void ZNum_sub(struct ZNum *a, struct ZNum *b, struct ZNum *result) {
    result->real = a->real - b->real;
    result->imag = a->imag - b->imag;
}

void ZNum_mul(struct ZNum *a, struct ZNum *b, struct ZNum *result) {
    result->real = a->real * b->real - a->imag * b->imag;
    result->imag = a->real * b->imag + a->imag * b->real;
}

void ZNum_div(struct ZNum *a, struct ZNum *b, struct ZNum *result) {
    /* Conjugate */
    struct ZNum conj = {b->real, -b->imag};
    double divisor = conj.real * conj.real + conj.imag * conj.imag;
    result->real = (a->real * conj.real - a->imag * conj.imag) / divisor;
    result->imag = (a->real * conj.imag + a->imag * conj.real) / divisor;
}

void ZNum_print(struct ZNum *z) {
    printf("(%f, %f)\n", z->real, z->imag);
}

void ZNum_root_of_unity(double radians, struct ZNum *result) {
    result->real = cos(radians);
    result->imag = sin(radians);
}
