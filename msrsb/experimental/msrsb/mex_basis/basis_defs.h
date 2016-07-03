typedef struct submatrix{
    int n;
    int *row;
    int *rowPos;
    double *diagonal;
    double *values;
}submatrix;

#ifndef max
 #define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
#endif
     
#ifndef min
 #define min(a,b) -max(-a, -b);
#endif

#ifndef __builtin_expect
  #define __builtin_expect(a, b) \
     a
#endif