# FFT
Cooley-Tukey implementation of DFT in C

Relevant functions are:

```c
void FFT_real_radix_2(double *samples, struct ZNum *coef, uint64_t length, int8_t direction);
void FFT_radix_2(struct ZNum *samples, struct ZNum *coef, uint64_t length, int8_t direction);
```

The first one can be used for purely real samples, and computes the DFT and saves it in an
array of complex numbers supplied by the caller.

The second one is for complex samples.

The "direction" of the transform determines whether it is a forward or inverse tranform.
Use the constants FFT_FORWARD and FFT_BACKWARD.
