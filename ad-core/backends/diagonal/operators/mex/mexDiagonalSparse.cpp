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

template <bool rowMajor>
void jac_to_sparse(const int l, const int n, const int m, const double* diagonal, const double* subset,
    double* pr, mwIndex* ir, mwIndex* jc) {
#pragma omp parallel for schedule(static)
    for (int index = 0; index < l * n; index++) {
        mwIndex ji, ii;
        double d;
        if (rowMajor) {
            ji = index;
            ii = index % m;
            int der = index / m;
            d = diagonal[ii * n + der];
        } else {
            ji = index;
            ii = index % m;
            d = diagonal[index];
        }
        // One entry per column
        jc[index] = ji;
        // Column index
        ir[index] = ii;
        // Actual derivative
        pr[index] = d;
    }
    // Final entry, one per column
    jc[l*n] = l*n;
}

template <int has_subset>
void sparse_rowmajor(const int l, const int n, const int m, const int ncols, const double* diagonal, const double* subset,
    double* pr, mwIndex* ir, mwIndex* jc) {
    if (!has_subset) {
#pragma omp parallel for schedule(static)
        for (int index = 0; index < l * n; index++) {
            // One entry per column
            jc[index] = index;
            // Column index
            ir[index] = index % m;
            // Actual derivative
            pr[index] = diagonal[index];
        }
        // Final entry, one per column
        jc[ncols] = ncols;
    }
    else {
#pragma omp parallel for schedule(static)
        for (int index = 0; index < m * n; index++) {
            // One entry per column
            jc[index] = index;
            // Column index
            ir[index] = index % m;
            // Actual derivative
            pr[index] = diagonal[index];
        }
        // Final entry, one per column
        jc[ncols] = ncols;
    }
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
    // s = mexDiagonalSparse(D.diagonal, D.subset, D.dim);
    // In: 
    // diagonal (m x n)
    // subset (either empty or length m)
    // dim (1 x 2) with entries [l n]
    // rowMajor (true for row major, false otherwise)
    // Out: sparse
    if (nrhs == 0) {
        if (nlhs > 0) {
            mexErrMsgTxt("Cannot give outputs with no inputs.");
        }
        // We are being called through compilation testing. Just do nothing. 
        // If the binary was actually called, we are good to go.
        return;
    } else if (nrhs != 4) {
	    mexErrMsgTxt("4 input arguments required."); 
    } else if (nlhs > 1) {
	    mexErrMsgTxt("Wrong number of output arguments."); 
    }
    double * diagonal = mxGetPr(prhs[0]);
    double * subset = mxGetPr(prhs[1]);
    double * dims = mxGetPr(prhs[2]);
    bool rowMajor = mxGetScalar(prhs[3]);

    int n; // Number of derivatives per element
    int m; // Number of elements
    // Number of rows
    int nrows = mxGetM(prhs[0]);
    // Number of columns
    int ncols = mxGetN(prhs[0]);
    if (rowMajor) {
        n = nrows;
        m = ncols;
    } else {
        n = ncols;
        m = nrows;
    }
    int dim_sz = mxGetN(prhs[2]);
    int l = dims[0];
    // mexPrintf("n = %d, m = %d, l = %d\n", n, m, l);
    if(dim_sz != 2){
        mexErrMsgTxt("Malformed dim array. Should be a 1 x 2 array.");
    }
    if(dims[1] != n){
        mexPrintf("Dimension in D.dim(1) is %d, but I expected %d from the dimensions of D.\n", l, n);
        mexErrMsgTxt("Mismatch in Jacobian.");
    }
    
    int subset_sz = mxGetM(prhs[1]);
    int m_l;
    if (subset_sz == 0) {
        m_l = l;
    }
    else if (subset_sz != m) {
        mexErrMsgTxt("Subset was given, but does not match diagonal.");
    }
    else {
        m_l = m;
    }
    plhs[0] = mxCreateSparse(m, l * n, m * n, mxREAL);
    // Entries
    double * pr  = mxGetPr(plhs[0]);
    // Row indices, zero-indexed (direct entries)
    mwIndex * ir = mxGetIr(plhs[0]);
    // Column indices, zero-indexed, offset encoded of length l*n + 1
    mwIndex * jc = mxGetJc(plhs[0]);

    if (rowMajor) {
        jac_to_sparse<true>(m_l, n, m, diagonal, subset, pr, ir, jc);
    }
    else {
        jac_to_sparse<false>(m_l, n, m, diagonal, subset, pr, ir, jc);
    }
}


