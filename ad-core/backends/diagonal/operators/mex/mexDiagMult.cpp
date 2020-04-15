#include <mex.h>

#ifndef HAVE_OCTAVE
#include "blas.h"
#endif

#include <omp.h>
#include <chrono>
#include <iostream>


template <int m, class V_type>
void diagMult(const size_t n, const V_type* v, const double* D, double* out) {
    diagMult(n, m, v, D, out);
}

template <class V_type>
void diagMult(const size_t n, const size_t m, const V_type * v, const double* D, double* out) {
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        const double vi = v[i];
        for (int j = i * m; j < (i + 1) * m; j++) {
            out[j] = D[j] * vi;
        }
    }
}


/* MEX gateway */
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
    // auto t0 = std::chrono::high_resolution_clock::now();
    // In: diagonal (m x nc), v (nc x 1)
    int n = mxGetN(prhs[0]);
    int m = mxGetM(prhs[0]);
    bool rowMajor = mxGetScalar(prhs[2]);

    if (!rowMajor) {
        mexErrMsgTxt("Column major not supported for this mex file.");
    }

    // auto t1 = std::chrono::high_resolution_clock::now();
    

    // auto t2 = std::chrono::high_resolution_clock::now();
    // plhs[0] = mxDuplicateArray(prhs[0]);
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    double* out_ptr = mxGetPr(plhs[0]);
    double* d_ptr = mxGetPr(prhs[0]);

    auto is_logical_vector = mxIsLogical(prhs[1]);
    if (!is_logical_vector){
        const double* v_ptr = mxGetPr(prhs[1]);
        diagMult(n, m, v_ptr, d_ptr, out_ptr);
    }
    else {
        const mxLogical* v_ptr = mxGetLogicals(prhs[1]);
        diagMult(n, m, v_ptr, d_ptr, out_ptr);
    }
}
