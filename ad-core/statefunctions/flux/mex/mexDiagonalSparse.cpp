//
// include necessary system headers
//
#include <cmath>
#include <mex.h>
#include <array>
#include <omp.h>
#include <iostream>
/* MEX gateway */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
    // s = mexDiagonalSparse(D.diagonal, D.subset, D.dim);
    // In: 
    // diagonal (n x m)
    // subset (either empty or length m)
    // dim (1 x 2) with entries [l n]
    // Out: sparse
    if (nrhs != 3) { 
	    mexErrMsgTxt("3 input arguments required."); 
    } else if (nlhs > 1) {
	    mexErrMsgTxt("Wrong number of output arguments."); 
    }
    double * diagonal = mxGetPr(prhs[0]);
    double * subset = mxGetPr(prhs[1]);
    double * dims = mxGetPr(prhs[2]);

    // Number of rows (number of derivatives per element)
    int n = mxGetM(prhs[0]);
    // Number of columns (number of elements)
    int m = mxGetN(prhs[0]);

    int dim_sz = mxGetN(prhs[2]);
    if(dim_sz != 2){
        mexErrMsgTxt("Malformed dim array. Should be a 1 x 2 array.");
    }
    if(dims[1] != n){
        mexErrMsgTxt("Mismatch in Jacobian.");
    }
    int l = dims[0];
    int subset_sz = mxGetM(prhs[1]);
    
    // std::cout << "Recieved diagonal " << m << "x" << n << "\n";
    // std::cout << "l=" << l << " m=" << m << " n=" << n << " subset_sz=" << subset_sz << "\n";
    int ncols = l*n;
    int nrows = m;
    plhs[0] = mxCreateSparse(nrows, ncols, m*n, mxREAL);
    // Entries
    double * pr  = mxGetPr(plhs[0]);
    // Row indices, zero-indexed (direct entries)
    mwIndex * ir = mxGetIr(plhs[0]);
    // Column indices, zero-indexed, offset encoded of length l*n + 1
    mwIndex * jc = mxGetJc(plhs[0]);

    if(subset_sz != 0){
        mexErrMsgTxt("Routine only supports empty subsets");
    }
    #pragma omp parallel for
    for(int index = 0; index < l*n; index++){
        // One entry per column
        int row = index % m;
        int der = index / m;
        jc[index] = index;
        // Column index
        ir[index] = index % m;
        // Actual derivative
        pr[index] = diagonal[row*n + der];
    }
    // Final entry, one per column
    jc[ncols] = ncols;
}


