//
// include necessary system headers
//
#include <cmath>
#include <mex.h>
#include <array>
#ifdef _OPENMP
    #include <omp.h>
#endif
#include <iostream>
// In: Cell value (nc x d), N (nf x 2), flag (nf x 1) bool
// Out: Face value of (nf x d)

template <class logic_type> 
void upwind(const int nf, const int nc, const int dim, const double* value, const double* N, const logic_type * flag, double* result) {
     #pragma omp parallel for
    for(int i=0;i<nf;i++){
        int fpos = i + nf;
        if(flag[i]){
            fpos = i;
        }
        else{
            fpos = i + nf;
        }
        int cell_inx = N[fpos] - 1;
        for(int j =0; j<dim; j++){
            result[i+nf*j] = value[cell_inx + nc*j];
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
    } else if (nin != 3) {
        status_code = -2;
        return "3 input arguments required. In: Cell value (nc x d), N (nf x 2), flag (nf x 1) bool";
    } else if (nout > 1) {
        status_code = -3;
        return "Wrong number of output arguments. Out: Face value of (nf x d)";
    } else {
        // All ok.
        status_code = 0;
        return "";
    }
}

/* MEX gateway */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
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
    mxLogical * flag = mxGetLogicals(prhs[2]);

    int dim = mxGetN(prhs[0]);
    int nc = mxGetM(prhs[0]);
    int nf = mxGetM(prhs[1]);
        
    plhs[0] = mxCreateUninitNumericMatrix(nf, dim, mxDOUBLE_CLASS, mxREAL);
    double * result = mxGetPr(plhs[0]);
    
    upwind<mxLogical>(nf, nc, dim, value, N, flag, result);
    return;
}


