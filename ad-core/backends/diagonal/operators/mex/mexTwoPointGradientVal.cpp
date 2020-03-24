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
    // In: value (nc x m), N (nf x 2)
    // Out: nf gradient
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

    int nc = mxGetM(prhs[0]);
    int nf = mxGetM(prhs[1]);
    int m = mxGetN(prhs[0]);
        
    plhs[0] = mxCreateDoubleMatrix(nf, m, mxREAL);
    double * result = mxGetPr(plhs[0]);

     #pragma omp parallel for schedule(static)
     for(int i=0;i<nf;i++){
         int left = N[i] - 1;
         int right = N[i + nf] - 1;
         for (int j = 0; j < m; j++) {
            result[i + j * nf] = value[right + j * nc] - value[left + j * nc];

        }
    }
     /*
    for(int j=0;j<m;j++){
        for(int i=0;i<2*nf;i++){
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

     */
    return;
}


