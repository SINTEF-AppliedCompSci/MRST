#include <omp.h>
#include <chrono>
#include <iostream>

#ifdef MRST_OCTEXT
    #include <octave/oct.h>
    #include <octave/dMatrix.h>
#else
    #include <mex.h>
    #ifndef HAVE_OCTAVE
        #include "blas.h"
    #endif
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

#ifdef MRST_OCTEXT
    /* OCT gateway */
    DEFUN_DLD (mexDiagProductMult, args, nargout,
               "mexDiagProductMult (see ADI.diagProductMult)")
    {
        const int nrhs = args.length();
        const int nlhs = nargout;
        if (nrhs == 0) {
            if (nlhs > 0) {
                error("Cannot give outputs with no inputs.");
            }
            // We are being called through compilation testing. Just do nothing. 
            // If the binary was actually called, we are good to go.
            return octave_value_list();
        }
        // D, v, B, w, rowMajor

        // Matrices
        const NDArray D_nd = args(0).array_value();
        const NDArray B_nd = args(2).array_value();
        // Vectors
        const NDArray v_nd = args(1).array_value();
        const NDArray w_nd = args(3).array_value();

        // D
        int m = D_nd.rows();
        int n = D_nd.cols();
        // D+B ptrs
        const double * d_ptr = D_nd.data();
        const double * b_ptr = B_nd.data();
        // v+w ptrs
        const double* v_ptr = v_nd.data();
        const double* w_ptr = w_nd.data();

        bool rowMajor = args(4).scalar_value();

        if (!rowMajor) {
            error("Column major not supported for this mex file.");
        }

        NDArray output({m, n});
        double * out_ptr = output.fortran_vec();

        diagProductMult(n, m, v_ptr, d_ptr, w_ptr, b_ptr, out_ptr);
        return octave_value (output);
    }
#else
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

        plhs[0] = mxCreateUninitNumericMatrix(m, n, mxDOUBLE_CLASS, mxREAL);
        double* out_ptr = mxGetPr(plhs[0]);
        diagProductMult(n, m, v_ptr, d_ptr, w_ptr, b_ptr, out_ptr);
    };
#endif

