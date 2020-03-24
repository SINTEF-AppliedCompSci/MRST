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

/* MEX gateway */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
    // In: Cell value (nc x d), N (nf x 2)
    // Out: Face value of (nf x d)
    if (nrhs == 0) {
        if (nlhs > 0) {
            mexErrMsgTxt("Cannot give outputs with no inputs.");
        }
        // We are being called through compilation testing. Just do nothing. 
        // If the binary was actually called, we are good to go.
        return;
    } else if (nrhs != 2) {
	    mexErrMsgTxt("2 input arguments required."); 
    } else if (nlhs > 1) {
	    mexErrMsgTxt("Wrong number of output arguments."); 
    } 
    double * value = mxGetPr(prhs[0]);
    double * N = mxGetPr(prhs[1]);

    int dim = mxGetN(prhs[0]);
    int nc = mxGetM(prhs[0]);
    int nf = mxGetM(prhs[1]);
    plhs[0] = mxCreateDoubleMatrix(nf, dim, mxREAL);
    double * result = mxGetPr(plhs[0]);
    
    #pragma omp parallel for schedule(static)
    for(int i=0;i<nf;i++){
        int left = N[i] - 1;
        int right = N[i + nf] - 1;
        for(int j =0; j<dim; j++){
            result[i+nf*j] = 0.5*(value[left + nc*j] + value[right + nc * j]);
        }
    }
    return;
}


