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
    // In: Cell value (nc x d), N (nf x 2), flag (nf x 1) bool
    // Out: Face value of (nf x d)
    if (nrhs != 3) { 
	    mexErrMsgTxt("3 input arguments required."); 
    } else if (nlhs > 1) {
	    mexErrMsgTxt("Wrong number of output arguments."); 
    } 
    double * value = mxGetPr(prhs[0]);
    double * N = mxGetPr(prhs[1]);
    mxLogical  * flag = mxGetLogicals(prhs[2]);

    int dim = mxGetN(prhs[0]);
    int nc = mxGetM(prhs[0]);
    int nf = mxGetM(prhs[1]);
    
    // printf("%d cells with %d faces and %d derivatives \n", nc, nf, m);
    
    plhs[0] = mxCreateDoubleMatrix(nf, dim, mxREAL);
    double * result = mxGetPr(plhs[0]);
    
    #pragma omp parallel for
    for(int i=0;i<nf;i++){
        int cell_inx;
        if(flag[i]){
            cell_inx = N[i]-1;
        }
        else{
            cell_inx = N[i + nf]-1;
        }
        for(int j =0; j<dim; j++){
            result[i+nf*j] = value[cell_inx + nc*j];
        }
    }
    return;
}


