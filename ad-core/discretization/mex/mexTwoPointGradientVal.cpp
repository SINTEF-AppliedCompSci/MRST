//
// include necessary system headers
//
#include <cmath>
#include <mex.h>
#include <array>

#include <iostream>

/* MEX gateway */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
    int i, j;
    // In: value (nc x m), N (nf x 2)
    // Out: nf gradient
    if (nrhs != 2) { 
	    mexErrMsgTxt("2 input arguments required."); 
    } else if (nlhs > 1) {
	    mexErrMsgTxt("Wrong number of output arguments."); 
    } 
    double * value = mxGetPr(prhs[0]);
    double * N = mxGetPr(prhs[1]);

    int nc = mxGetM(prhs[0]);
    int nf = mxGetM(prhs[1]);
    int m = mxGetN(prhs[0]);
        
    plhs[0] = mxCreateDoubleMatrix(nf, m, mxREAL);
    double * result = mxGetPr(plhs[0]);
    
    for(j=0;j<m;j++){
        #pragma omp parallel for
        for(i=0;i<2*nf;i++){
            int cell_inx = N[i]-1;
            double v = value[cell_inx +  j*nc];
            if(i<nf){
                #pragma omp atomic
                result[i + j*nf] -= v;
            }else{
                #pragma omp atomic
                result[i%nf + j*nf] += v;
            }
        }
    }
    return;
}


