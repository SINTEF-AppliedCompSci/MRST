//
// include necessary system headers
//
#include <cmath>
#include <array>
#ifdef _OPENMP
    #include <omp.h>
#endif
#include <iostream>
#ifdef MRST_OCTEXT
    #include <octave/oct.h>
    #include <octave/dMatrix.h>
#else
    #include <mex.h>
#endif

void faceGradient(const int nf, const int nc, const int dim, const double* value, const double* N, double* result) {
     #pragma omp parallel for
     for(int i=0;i<nf;i++){
         int left = N[i] - 1;
         int right = N[i + nf] - 1;
         for (int j = 0; j < dim; j++) {
            result[i + j * nf] = value[right + j * nc] - value[left + j * nc];
        }
    }
}

const char* inputCheck(const int nin, const int nout, int & status_code){
    if (nin == 0) {
        if (nout > 0) {
            status_code = -1;
            return "Cannot give outputs with no inputs.";
        }
        // We are being called through compilation testing. Just do nothing.
        // If the binary was actually called, we are good to go.
        status_code = 1;
        return "";
    } else if (nin != 2) {
        status_code = -2;
        return "2 input arguments required.";
    } else if (nout > 1) {
        status_code = -3;
        return "Wrong number of output arguments.";
    } else {
        // All ok.
        status_code = 0;
        return "";
    }
}


#ifdef MRST_OCTEXT
    /* OCT gateway */
    DEFUN_DLD (mexTwoPointGradientVal, args, nargout,
               "Face gradient operator for MRST - value only.")
    {
        const int nrhs = args.length();
        const int nlhs = nargout;
        int status_code = 0;
        auto msg = inputCheck(nrhs, nlhs, status_code);
        if(status_code < 0){
            // Some kind of error
            error(msg);
        }else if (status_code == 1){
            // Early return
            return octave_value_list();
        }

        const NDArray value_nd = args(0).array_value();
        const NDArray N_nd = args(1).array_value();
        
        const double * value = value_nd.data();
        const double * N = N_nd.data();
        const dim_vector sz = value_nd.dims();
        const dim_vector nsz = N_nd.dims();
        
        int dim = sz(1);
        int nc = sz(0);
        
        int nf = nsz(0);

        NDArray output({nf, dim});
        double * result = output.fortran_vec();

       faceGradient(nf, nc, dim, value, N, result);
       return octave_value (output);
    }
#else
    /* MEX gateway */
    void mexFunction( int nlhs, mxArray *plhs[], 
              int nrhs, const mxArray *prhs[] )

    { 
        // In: Cell value (nc x d), N (nf x 2)
        // Out: Face value of (nf x d)
        int status_code = 0;
        auto msg = inputCheck(nrhs, nlhs, status_code);
        if(status_code < 0){
            // Some kind of error
            mexErrMsgTxt(msg);
        }else if (status_code == 1){
            // Early return
            return;
        }

        double * value = mxGetPr(prhs[0]);
        double * N = mxGetPr(prhs[1]);

        int dim = mxGetN(prhs[0]);
        int nc = mxGetM(prhs[0]);
        int nf = mxGetM(prhs[1]);
        // plhs[0] = mxCreateDoubleMatrix(nf, dim, mxREAL);
        plhs[0] = mxCreateUninitNumericMatrix(nf, dim, mxDOUBLE_CLASS, mxREAL);
        double * result = mxGetPr(plhs[0]);
        faceGradient(nf, nc, dim, value, N, result);

        return;
    }
#endif

