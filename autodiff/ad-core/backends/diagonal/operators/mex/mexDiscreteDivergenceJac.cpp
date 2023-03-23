//
// include necessary system headers
//
#include <cmath>
#include <array>
#ifdef _OPENMP
    #include <omp.h>
#endif
#include <iostream>
#ifdef MRST_OCTEXT
    #include <octave/oct.h>
    #include <octave/dMatrix.h>
    #define octix octave_idx_type
#else
    #include <mex.h>
#endif

// In: 
// acc diagonal (nc x m) or empty
// face diagonal (nf x 2*m)
// N (nf x 2)
// facePos (nc+1 x 1)
// faces (length facePos(end)-1)
// cells
// cells_ix
// rowMajor (boolean)
// Out: either single sparse or the five sparse constructor inputs (I, J, V, n, m)
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
    } else if (nin != 8) {
        status_code = -2;
        return "8 input arguments required.";
    } else if (nout != 5 && nout != 1) {
        status_code = -3;
        return "Function must either produce five dense outputs (I, J, V, n, m) or a single sparse output (sparse matrix).";
    } else {
        // All ok.
        status_code = 0;
        return "";
    }
}

template <bool colMajor, bool lower, class index_t>
void copyFaceData(const int c, const int nf, const int m, const int diag, const int passed, 
                const int sparse_mult, const int cell_offset, const int f, const int fl, const double* diagonal, double* pr, index_t* ir, index_t* jc) {
    for (int der = 0; der < m; der++) {
        double v;
        int sparse_offset = der * sparse_mult + cell_offset;
        if (lower) { // c < 0
            // Low entry, corresponding to N(f, 1)
            if (colMajor) {
                v = -diagonal[der * nf + f];
            }
            else {
                v = -diagonal[f * 2 * m + der];
            }
            ir[sparse_offset + fl + passed] = -c;
        }
        else {
            // High entry, corresponding to N(f, 2)
            if (colMajor) {
                v = diagonal[der * nf + f + m * nf];
            }
            else {
                v = diagonal[(f * 2 + 1) * m + der];
            }
            // Set row entry
            ir[sparse_offset + fl + passed] = c;
        }
        pr[sparse_offset + fl + passed] = v;
        // Set corresponding diagonal entry
        pr[sparse_offset + diag] -= v;
    }
}

template <bool has_accumulation, bool colMajor, class index_t>
void divergenceJac(const int nf, const int nc, const int m,
    const double* N, const double* facePos, const double* faces,
    const double* cells, const double* cells_ix,
    const double* accumulation, const double* diagonal,
    double* pr, index_t* ir, index_t* jc) {
    int mv = facePos[nc];
    #pragma omp parallel for
    for (int cell = 0; cell < nc; cell++) {
        // Each cell has number of connections equal to the number of half-
        // faces for that cell plus itself multiplied by the block size
        int f_offset = facePos[cell];
        int n_local_hf = facePos[cell + 1] - f_offset;
        int diag = cells_ix[cell];
        int cell_offset = f_offset + cell;
        for (int der = 0; der < m; der++) {
            int ix = cell + der * nc;
            // Base offset taking into account how far we have come
            int base = der * (mv + nc) + cell_offset;
            jc[cell + der * nc + 1] = base + n_local_hf + 1;
            // Set diagonal entries
            int dpos = base + cells_ix[cell];
            ir[dpos] = cell;
            if (has_accumulation) {
                if (colMajor) {
                    pr[dpos] = accumulation[der * nc + cell];
                }
                else {
                    pr[dpos] = accumulation[cell * m + der];
                }
            }
            else {
                /* Not sure if this bit is needed - but the Matlab docs
                   are vague on memory initialization for sparse arrays */
                pr[dpos] = 0.0; 
            }
        }

        for (int fl = 0; fl < n_local_hf; fl++) {
            // Loop over entire column
            // Diagonal entry can be skipped - we handle this later
            // Check if we have passed diagonal entry
            int passed = (double)(fl >= diag);
            // Global face index
            int f = faces[f_offset + fl];
            // Global cell index
            int c = cells[f_offset + fl];
            // Iterate over derivatives
            int sparse_mult = mv + nc;
            // Copy the data from face jacobian stored in diagonal
            if (c < 0) {
                copyFaceData<colMajor, true>(c+1, nf, m, diag, passed, sparse_mult, cell_offset, f, fl, diagonal, pr, ir, jc);
            }
            else {
                copyFaceData<colMajor, false>(c-1, nf, m, diag, passed, sparse_mult, cell_offset, f, fl, diagonal, pr, ir, jc);
            }
        }
    }
}

template <int m, bool has_accumulation, bool colMajor, class index_t>
void divergenceJac(const int nf, const int nc,
    const double* N, const double* facePos, const double* faces,
    const double* cells, const double* cells_ix,
    const double* accumulation, const double* diagonal,
    double* pr, index_t* ir, index_t* jc) {
        divergenceJac<has_accumulation, colMajor>(nf, nc, m, N, facePos, faces, cells, cells_ix, accumulation, diagonal, pr, ir, jc);
}


template <bool has_accumulation, bool colMajor, class index_t>
void divergenceJacMain(const int nf, const int nc, const int m,
    const double* N, const double* facePos, const double* faces,
    const double* cells, const double* cells_ix,
    const double* accumulation, const double* diagonal,
    double* pr, index_t* ir, index_t* jc){
    switch (m) {
    case 1:
        divergenceJac<1, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 2:
        divergenceJac<2, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 3:
        divergenceJac<3, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 4:
        divergenceJac<4, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 5:
        divergenceJac<5, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 6:
        divergenceJac<6, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 7:
        divergenceJac<7, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 8:
        divergenceJac<8, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 9:
        divergenceJac<9, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 10:
        divergenceJac<10, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 11:
        divergenceJac<11, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 12:
        divergenceJac<12, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 13:
        divergenceJac<13, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 14:
        divergenceJac<14, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 15:
        divergenceJac<15, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 16:
        divergenceJac<16, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 17:
        divergenceJac<17, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 18:
        divergenceJac<18, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 19:
        divergenceJac<19, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 20:
        divergenceJac<20, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 21:
        divergenceJac<21, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 22:
        divergenceJac<22, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 23:
        divergenceJac<23, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 24:
        divergenceJac<24, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 25:
        divergenceJac<25, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 26:
        divergenceJac<26, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 27:
        divergenceJac<27, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 28:
        divergenceJac<28, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 29:
        divergenceJac<29, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    case 30:
        divergenceJac<30, has_accumulation, colMajor>(nf, nc, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
        break;
    default:
        divergenceJac<has_accumulation, colMajor>(nf, nc, m, N, facePos, faces, cells, cells_ix,
            accumulation, diagonal, pr, ir, jc);
}

}

#ifdef MRST_OCTEXT
    /* OCT gateway */
    DEFUN_DLD (mexDiscreteDivergenceJac, args, nargout,
               "Discrete divergence operator for MRST - diagonal Jacobian.")
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
        bool output_sparse = nargout < 2;
        // Jacobians
        const NDArray accJac   = args(0).array_value();
        const NDArray faceJac  = args(1).array_value();

        const double * accumulation = accJac.data();
        const double * diagonal = faceJac.data();
        // Grid structure
        const NDArray N_nd   = args(2).array_value();
        const double * N       = N_nd.data();
        const NDArray facePos_nd = args(3).array_value();
        const double * facePos = facePos_nd.data();
        const NDArray faces_nd = args(4).array_value();
        const double * faces   = faces_nd.data();
        const double * cells   = args(5).array_value().data();
        const double * cells_ix = args(6).array_value().data();
        bool rowMajor = args(7).scalar_value();

        octix nrows = faceJac.rows();
        octix ncols = faceJac.cols();

        octix nf = N_nd.rows();
        octix nc = facePos_nd.rows()-1;
        octix nhf = faces_nd.rows();
        octix n_acc = accJac.numel();

        bool has_accumulation = n_acc > 0;

        octix m, nf_in;
        if (rowMajor) {
            m = nrows / 2;
            nf_in = ncols;
        } else {
            m = ncols / 2;
            nf_in = nrows;
        }
        if(has_accumulation && n_acc != m*nc){
            error("Accumulation term diagonal was provided, but dimensions are incorrect.");
        }

        if(nf_in != nf){
            error("Face inputs do not match.");
        }
        // Each cell has one self-connection plus the number of half-faces, multiplied by block size
        octix nzmax = (facePos[nc] + nc)*m;
        double * pr;
        octix * ir;
        octix * jc;
        SparseMatrix jacobian (nc, nc*m, nzmax);
        pr = jacobian.data();
        ir = jacobian.ridx();
        jc = jacobian.cidx();

        if (!output_sparse) {
            error("Non-sparse output not yet supported in OCT-mode.");
        }

        if (rowMajor){
            if(has_accumulation){
                divergenceJacMain<true, false>(nf, nc, m, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
            }else{
                divergenceJacMain<false, false>(nf, nc, m, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
            }
        } else {
            if(has_accumulation){
                divergenceJacMain<true, true>(nf, nc, m, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
            }else{
                divergenceJacMain<false, true>(nf, nc, m, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
            }
        }
        return octave_value(jacobian);
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
        // Parse inputs
        bool output_sparse = nlhs < 2;
        const mxArray * accJac   = prhs[0];
        const mxArray * faceJac  = prhs[1];
        const double * N       = mxGetPr(prhs[2]);
        const double * facePos = mxGetPr(prhs[3]);
        const double * faces   = mxGetPr(prhs[4]);
        const double * cells   = mxGetPr(prhs[5]);
        const double * cells_ix = mxGetPr(prhs[6]);
        bool rowMajor = mxGetScalar(prhs[7]);

        // Dimensions of diagonals - figure out if we want row or column major solver
        int nrows = mxGetM(faceJac);
        int ncols = mxGetN(faceJac);


        /* Diagonal of accumulation term */
        const double * accumulation = mxGetPr(accJac);
        /* Diagonal of fluxes */
        const double * diagonal = mxGetPr(faceJac);

        int nf = mxGetM(prhs[2]);
        int nc = mxGetM(prhs[3])-1;
        int nhf = mxGetM(prhs[4]);
        int n_acc = mxGetNumberOfElements(accJac);
        
        bool has_accumulation = n_acc > 0;
        
        int m, nf_in;
        if (rowMajor) {
            m = nrows / 2;
            nf_in = ncols;
        } else {
            m = ncols / 2;
            nf_in = nrows;
        }
        // mexPrintf("Matrix has dimensions %d by %d. There are %d faces and %d cells (%d interfaces, %d derivatives)\n", nrows, ncols, nf, nc, nf_in, m);
        if(has_accumulation && n_acc != m*nc){
            mexErrMsgTxt("Accumulation term diagonal was provided, but dimensions are incorrect.");
        }

        if(nf_in != nf){
            mexErrMsgTxt("Face inputs do not match.");
        }
        // Each cell has one self-connection plus the number of half-faces, multiplied by block size
        mwSize nzmax = (facePos[nc] + nc)*m;
        // Row indices, zero-indexed (direct entries)
        mwIndex* ir;
        // Column indices, zero-indexed, offset encoded of length m*nc + 1
        mwIndex* jc;
        // Entries
        double* pr;
        if (output_sparse) {
            plhs[0] = mxCreateSparse(nc, nc * m, nzmax, mxREAL);

            pr = mxGetPr(plhs[0]);
            ir = mxGetIr(plhs[0]);
            jc = mxGetJc(plhs[0]);
        }
        else {
            plhs[0] = mxCreateNumericMatrix(nzmax, 1, mxUINT64_CLASS, mxREAL); // We have no way of allocating mwIndex (size_t). So we hope for the best and allocate uint64...
            plhs[1] = mxCreateNumericMatrix(m*nc+1, 1, mxUINT64_CLASS, mxREAL);
            plhs[2] = mxCreateDoubleMatrix(nzmax, 1, mxREAL);
            ir = (mwIndex*)mxGetData(plhs[0]);
            jc = (mwIndex*)mxGetData(plhs[1]);
            // Entries
            pr = mxGetPr(plhs[2]);
            plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
            plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);

            double* three = mxGetPr(plhs[3]);
            three[0] = nc;
            double* four = mxGetPr(plhs[4]);
            four[0] = nc * m;
        }
        if (rowMajor){
            if(has_accumulation){
                divergenceJacMain<true, false>(nf, nc, m, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
            }else{
                divergenceJacMain<false, false>(nf, nc, m, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
            }
        } else {
            if(has_accumulation){
                divergenceJacMain<true, true>(nf, nc, m, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
            }else{
                divergenceJacMain<false, true>(nf, nc, m, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
            }
        }
    }
#endif
