#include <cmath>
#include <mex.h>
#include <array>
#include "mexInterpolation.cpp"

#include <iostream>

/* MEX gateway */

void mexFunction( int nlhs, mxArray *plhs[], 
                  int nrhs, const mxArray *prhs[] )
     
{ 
    // In: X (n x 1), F (n x 1), x (m x 1)
    // Out Fx (m x 1) dFdx (m x 1)
    if (nrhs == 0 && nlhs == 0){
        mexPrintf("Interpolation is ready for use.\n");
        return;
    }else if (nrhs < 3) { 
	    mexErrMsgTxt("3 or 4 input arguments required."); 
    } else if (nlhs > 2) {
	    mexErrMsgTxt("Wrong number of output arguments."); 
    }
    size_t n = mxGetNumberOfElements(prhs[0]);
    size_t n2 = mxGetNumberOfElements(prhs[0]);
    if(n != n2){
        mexErrMsgTxt("Function and sample points have mismatching dimensions.");
    }
    size_t m = mxGetNumberOfElements(prhs[2]);
    double * X = mxGetPr(prhs[0]);
    double * F = mxGetPr(prhs[1]);
    double * x = mxGetPr(prhs[2]);
    double method;
    if(nrhs == 3){
        method = 1;
    }else{
        method = mxGetScalar(prhs[3]);
    }

    plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(m, 1, mxREAL);

    double * Fx = mxGetPr(plhs[0]);
    double * dFdx = mxGetPr(plhs[1]);
    switch((int)method){
        case 1:
            interp1_binary_search(X, F, x, Fx, dFdx, n, m);
            break;
        case 2:
            interp1_binned_search(X, F, x, Fx, dFdx, n, m);
            break;
        case 3:
            double dx = X[1]-X[0];
            // mexPrintf("Bin size: %f\n", dx);
            interp1_equal_width_search(X[0], dx, F, x, Fx, dFdx, n, m);
            break;
    }
}
