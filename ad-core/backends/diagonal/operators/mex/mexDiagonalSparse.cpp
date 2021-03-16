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

template <bool rowMajor, typename indexType>
void jac_to_sparse(const int l, const int n, const int m, const double* diagonal, const double* subset,
    double* pr, indexType* ir, indexType* jc) {
#pragma omp parallel for
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
    } else if (nlhs != 5 && nlhs != 1 && nlhs != 0) {
	    mexErrMsgTxt("Function must either produce five dense outputs (I, J, V, n, m) or a single sparse output (sparse matrix)."); 
    }
    bool output_sparse = nlhs < 2;

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
    mwIndex* ir;
    mwIndex* jc;
    double* pr;
    if (output_sparse) {
        plhs[0] = mxCreateSparse(m, l * n, m * n, mxREAL);
        // Row indices, zero-indexed (direct entries)
        ir = mxGetIr(plhs[0]);
        // Column indices, zero-indexed, offset encoded of length l*n + 1
        jc = mxGetJc(plhs[0]);
        pr = mxGetPr(plhs[0]);
        if (rowMajor) {
            jac_to_sparse<true>(m_l, n, m, diagonal, subset, pr, ir, jc);
        }
        else {
            jac_to_sparse<false>(m_l, n, m, diagonal, subset, pr, ir, jc);
        }
    }
    else {
        int nrows_out = m * n;
        plhs[0] = mxCreateNumericMatrix(nrows_out, 1, mxUINT64_CLASS, mxREAL); // We have no way of allocating mwIndex (size_t). So we hope for the best and allocate uint64...
        plhs[1] = mxCreateNumericMatrix(l*n + 1, 1, mxUINT64_CLASS, mxREAL);
        plhs[2] = mxCreateDoubleMatrix(nrows_out, 1, mxREAL);
        // Row indices, zero-indexed (direct entries)
        ir = (mwIndex *) mxGetData(plhs[0]);
        // Column indices, zero-indexed, offset encoded of length l*n + 1
        jc = (mwIndex *) mxGetData(plhs[1]);
        // Entries
        pr = mxGetPr(plhs[2]);
        plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
        plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);
        
        double* three = mxGetPr(plhs[3]);
        three[0] = m;
        double* four = mxGetPr(plhs[4]);
        four[0] = l*n;
        if (rowMajor) {
            jac_to_sparse<true>(m_l, n, m, diagonal, subset, pr, ir, jc);
        }
        else {
            jac_to_sparse<false>(m_l, n, m, diagonal, subset, pr, ir, jc);
        }
    }

}


