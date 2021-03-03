//
// include necessary system headers
//
#include <cmath>
#include <array>
#ifdef _OPENMP
    #include <omp.h>
#endif
#include <iostream>
#include <chrono>
#ifdef MRST_OCTEXT
    #include <octave/oct.h>
#else
    #include <mex.h>
#endif

// In: 
// acc (nc x 1) or empty
// flux (nf x 1)
// N (nf x 2)
// nc (scalar)
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
    } else if (nin != 6) {
        status_code = -2;
        return "6 input arguments required.";
    } else if (nout > 1) {
        status_code = -3;
        return "Too many outputs requested. Function has a single output argument.";
    } else {
        // All ok.
        status_code = 0;
        return "";
    }
}

template <bool has_accumulation> // nc, accumulation, flux, faces, facePos, N, result
void divergenceVal(const int nc, const double * accumulation, const double * flux, const double * faces,
                   const double * facePos, const double * N, double * result){
    #pragma omp parallel for
    for (int cell = 0; cell < (int)nc; cell++) {
        // Each cell has number of connections equal to the number of half-
        // faces for that cell plus itself multiplied by the block size
        double v;
        if(has_accumulation){
            v = accumulation[cell];
        }else{
            v = 0;
        }
        for (int i = facePos[cell]; i < facePos[cell+1]; i++){
            int f = faces[i];
            double f_v = flux[f];
            if(N[f] == cell+1){
                // Positive flux is out
                v += f_v;
            }else{
                v -= f_v;
            }
        }
        result[cell] = v;
    }
}
#ifdef MRST_OCTEXT
    /* OCT gateway */
    DEFUN_DLD (mexDiscreteDivergenceVal, args, nargout,
               "Discrete divergence for MRST - diagonal Jacobian.")
    {
        int status_code = 0;
        auto msg = inputCheck(args.length(), nargout, status_code);
        
        if(status_code < 0){
            // Some kind of error
            error(msg);
        } else if (status_code == 1){
            // Early return
            return octave_value_list();
        }
        // Accumulation term (may be empty)
        const NDArray acc_nd = args(0).array_value();
        const double * acc = acc_nd.data();
        int n_acc = acc_nd.numel();
        // Fluxes
        const NDArray v_nd = args(1).array_value();
        const double * flux = v_nd.data();
        // Dimensions etc
        const NDArray N_nd = args(2).array_value();
        const double * N = N_nd.data();

        const double nc = args(3).scalar_value();

        const NDArray facePos_nd = args(4).array_value();
        const double * facePos = facePos_nd.data();
        const NDArray faces_nd = args(5).array_value();
        const double * faces = faces_nd.data();

        NDArray output({(int)nc, 1});
        double * result = output.fortran_vec();

        if(n_acc > 0){
            divergenceVal<true>(nc, acc, flux, faces, facePos, N, result);
        }else{
            divergenceVal<false>(nc, acc, flux, faces, facePos, N, result);
        }
        return octave_value (output);
    }
#else
    /* MEX gateway */
    void mexFunction( int nlhs, mxArray *plhs[], 
            int nrhs, const mxArray *prhs[] )
        
    { 
        int status_code = 0;
        auto msg = inputCheck(nrhs, nlhs, status_code);
        if(status_code < 0){
            // Some kind of error
            mexErrMsgTxt(msg);
        } else if (status_code == 1){
            // Early return
            return;
        }
        // Accumulation term (may be empty)
        const mxArray * acc = prhs[0];
        const double * accumulation = mxGetPr(prhs[0]);
        int n_acc = mxGetNumberOfElements(acc);
        // Fluxes
        const mxArray * v = prhs[1];
        const double * flux = mxGetPr(v);
        // Dimensions etc
        const double * N = mxGetPr(prhs[2]);
        const double nc = mxGetScalar(prhs[3]);
        const double * facePos = mxGetPr(prhs[4]);
        const double * faces   = mxGetPr(prhs[5]);

        // int nf = mxGetM(prhs[1]);
        
        // plhs[0] = mxCreateDoubleMatrix(nc, 1, mxREAL);
        plhs[0] = mxCreateUninitNumericMatrix(nc, 1, mxDOUBLE_CLASS, mxREAL);
        double * result = mxGetPr(plhs[0]);
        if(n_acc > 0){
            divergenceVal<true>(nc, accumulation, flux, faces, facePos, N, result);
        }else{
            divergenceVal<false>(nc, accumulation, flux, faces, facePos, N, result);
        }
    }
#endif


