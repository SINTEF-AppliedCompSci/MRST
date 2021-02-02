#include <mex.h>

#ifndef HAVE_OCTAVE
#include "blas.h"
#endif

#include <omp.h>
#include <chrono>
#include <iostream>

#if !defined(_WIN32)
#define dscal dscal_
#endif

// #define operation diagMultBLAS
#define operation diagProductMult

void diagProductMult(const size_t n, const size_t m,
                                     const double* v1, const double* D1,
                                     const double* v2, const double* D2, double* out) {
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        auto v1_i = v1[i];
        auto v2_i = v2[i];
        for (int j = m * i; j < m * (i + 1); j++) {
            out[j] = v1_i * D1[j] + v2_i * D2[j];
        }
    }
}

template <int m>
void diagProductMult(const size_t n,
    const double* v1, const double* D1,
    const double* v2, const double* D2, double* out) {
    diagProductMult(n, m, v1, D1, v2, D2, out);
}


void diagProductMultScalar(const size_t n, const double* v1, const double* D1,
                                           const double* v2, const double* D2, double* out) {
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        out[i] = v1[i]*D2[i] +  v1[i]*D2[i];
    }
}


void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[])

{
    if (nrhs == 0) {
        if (nlhs > 0) {
            mexErrMsgTxt("Cannot give outputs with no inputs.");
        }
        // We are being called through compilation testing. Just do nothing. 
        // If the binary was actually called, we are good to go.
        return;
    }
    // In: diagonal (m x nc), v (nc x 1)
    int n = mxGetN(prhs[0]);
    int m = mxGetM(prhs[0]);
    bool rowMajor = mxGetScalar(prhs[4]);

    if (!rowMajor) {
        mexErrMsgTxt("Column major not supported for this mex file.");
    }
    const double* v_ptr = mxGetPr(prhs[1]);
    const double* w_ptr = mxGetPr(prhs[3]);

    double* d_ptr = mxGetPr(prhs[0]);
    double* b_ptr = mxGetPr(prhs[2]);

    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    double* out_ptr = mxGetPr(plhs[0]);
    diagProductMult(n, m, v_ptr, d_ptr, w_ptr, b_ptr, out_ptr);
};

