struct ZNum {
    double real;
    double imag; /* In PI Radians */
};

void ZNum_add(struct ZNum *a, struct ZNum *b, struct ZNum *result);
void ZNum_sub(struct ZNum *a, struct ZNum *b, struct ZNum *result);
void ZNum_mul(struct ZNum *a, struct ZNum *b, struct ZNum *result);
void ZNum_div(struct ZNum *a, struct ZNum *b, struct ZNum *result);
void ZNum_print(struct ZNum *z);
void ZNum_root_of_unity(double radians, struct ZNum *result);
