#include <stdint.h>
#include "complex.h"

#define FFT_FORWARD -1
#define FFT_BACKWARD 1

/* Samples are always real, coefficients are complex */
void FFT_real_radix_2(double *samples, struct ZNum *coef, uint64_t length, int8_t direction);

/* Both samples and coefficients are complex */
void FFT_radix_2(struct ZNum *samples, struct ZNum *coef, uint64_t length, int8_t direction);

