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

template <int m, class V_type>
void diagMult(const size_t n, const V_type* v, const double* D, double* out) {
    diagMult(n, m, v, D, out);
}

template <class V_type>
void diagMult(const size_t n, const size_t m, const V_type * v, const double* D, double* out) {
#pragma omp parallel for
    for (int i = 0; i < n; i++) {
        const double vi = v[i];
        for (int j = i * m; j < (i + 1) * m; j++) {
            out[j] = D[j] * vi;
        }
    }
}

#ifdef MRST_OCTEXT
    /* OCT gateway */
    DEFUN_DLD (mexDiagMult, args, nargout,
               "mexDiagMult (see ADI.diagMult)")
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
        const NDArray diagonal_nd = args(0).array_value();
        const NDArray v_nd = args(1).array_value();

        // D
        int m = diagonal_nd.rows();
        int n = diagonal_nd.cols();
        const double * d_ptr = diagonal_nd.data();
        // v
        const double* v_ptr = v_nd.data();
        bool rowMajor = args(2).scalar_value();

        if (!rowMajor) {
            error("Column major not supported for this mex file.");
        }

        NDArray output({m, n});
        double * out_ptr = output.fortran_vec();

        diagMult(n, m, v_ptr, d_ptr, out_ptr);
        return octave_value (output);
    }
#else
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
        int n = mxGetN(prhs[0]);
        int m = mxGetM(prhs[0]);
        bool rowMajor = mxGetScalar(prhs[2]);
        const double* d_ptr = mxGetPr(prhs[0]);

        if (!rowMajor) {
            mexErrMsgTxt("Column major not supported for this mex file.");
        }

        plhs[0] = mxCreateUninitNumericMatrix(m, n, mxDOUBLE_CLASS, mxREAL);
        double* out_ptr = mxGetPr(plhs[0]);

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
#endif
