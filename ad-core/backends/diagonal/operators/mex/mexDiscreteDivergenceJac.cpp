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
#include <chrono>

template <int m, bool colMajor, bool lower>
void copyFaceData(const int c, const int nf, const int diag, const int passed, const int sparse_mult, const int cell_offset, const int f, const int fl, const double* diagonal, double* pr, mwIndex* ir, mwIndex* jc) {
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

template <int m, bool has_accumulation, bool colMajor>
void divergenceJac(const int nf, const int nc,
    const double* N, const double* facePos, const double* faces,
    const double* cells, const double* cells_ix,
    const double* accumulation, const double* diagonal,
    double* pr, mwIndex* ir, mwIndex* jc) {
    int mv = facePos[nc];
    #pragma omp parallel for schedule(dynamic, 512)
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
                copyFaceData<m, colMajor, true>(c, nf, diag, passed, sparse_mult, cell_offset, f, fl, diagonal, pr, ir, jc);
            }
            else {
                copyFaceData<m, colMajor, false>(c, nf, diag, passed, sparse_mult, cell_offset, f, fl, diagonal, pr, ir, jc);
            }
        }
    }
}

/* MEX gateway */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
    auto t0 = std::chrono::high_resolution_clock::now();
    // In: 
    // cell diagonal (nc x m) or empty
    // diagonal (nf x 2*m)
    // N (nf x 2)
    // facePos (nc+1 x 1)
    // faces (length facePos(end)-1)
    // rowMajor
    // Out: I, J, V
    if (nrhs == 0) {
        if (nlhs > 0) {
            mexErrMsgTxt("Cannot give outputs with no inputs.");
        }
        // We are being called through compilation testing. Just do nothing. 
        // If the binary was actually called, we are good to go.
        return;
    } else if (nrhs != 8) {
	    mexErrMsgTxt("8 input arguments required."); 
    } else if (nlhs != 5 && nlhs != 1 && nlhs != 0) {
        mexErrMsgTxt("Function must either produce five dense outputs (I, J, V, n, m) or a single sparse output (sparse matrix).");
    }
    bool output_sparse = nlhs < 2;
    const mxArray * accJac   = prhs[0];
    const mxArray * faceJac  = prhs[1];
    double * N       = mxGetPr(prhs[2]);
    double * facePos = mxGetPr(prhs[3]);
    double * faces   = mxGetPr(prhs[4]);
    double * cells   = mxGetPr(prhs[5]);
    double * cells_ix = mxGetPr(prhs[6]);
    bool rowMajor = mxGetScalar(prhs[7]);

    // Dimensions of diagonals - figure out if we want row or column major solver
    int nrows = mxGetM(faceJac);
    int ncols = mxGetN(faceJac);


    /* Diagonal of accumulation term */
    double * accumulation = mxGetPr(accJac);
    /* Diagonal of fluxes */
    double * diagonal = mxGetPr(faceJac);

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
    // printf("%d cells %d faces, %d half-faces and %d derivatives \n", nc, nf, nhf, m);
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
    if (rowMajor) {
        // RowMajor
        if (has_accumulation) {
            // mexPrintf("Row major, with accumulation: %d by %d.\n", nc, nf);
            switch (m) {
                case 1:
                    divergenceJac<1, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 2:
                    divergenceJac<2, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 3:
                    divergenceJac<3, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 4:
                    divergenceJac<4, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 5:
                    divergenceJac<5, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 6:
                    divergenceJac<6, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 7:
                    divergenceJac<7, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 8:
                    divergenceJac<8, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 9:
                    divergenceJac<9, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 10:
                    divergenceJac<10, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 11:
                    divergenceJac<11, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 12:
                    divergenceJac<12, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 13:
                    divergenceJac<13, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 14:
                    divergenceJac<14, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 15:
                    divergenceJac<15, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 16:
                    divergenceJac<16, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 17:
                    divergenceJac<17, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 18:
                    divergenceJac<18, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 19:
                    divergenceJac<19, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 20:
                    divergenceJac<20, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 21:
                    divergenceJac<21, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 22:
                    divergenceJac<22, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 23:
                    divergenceJac<23, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 24:
                    divergenceJac<24, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 25:
                    divergenceJac<25, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 26:
                    divergenceJac<26, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 27:
                    divergenceJac<27, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 28:
                    divergenceJac<28, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 29:
                    divergenceJac<29, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 30:
                    divergenceJac<30, true, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                default:
                    mexErrMsgTxt("%d derivatives not supported by backend.");
            }
        }
        else {
            // mexPrintf("Row major, without accumulation: %d by %d.\n", nc, nf);
            switch (m) {
                case 1:
                    divergenceJac<1, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 2:
                    divergenceJac<2, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 3:
                    divergenceJac<3, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 4:
                    divergenceJac<4, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 5:
                    divergenceJac<5, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 6:
                    divergenceJac<6, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 7:
                    divergenceJac<7, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 8:
                    divergenceJac<8, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 9:
                    divergenceJac<9, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 10:
                    divergenceJac<10, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 11:
                    divergenceJac<11, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 12:
                    divergenceJac<12, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 13:
                    divergenceJac<13, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 14:
                    divergenceJac<14, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 15:
                    divergenceJac<15, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 16:
                    divergenceJac<16, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 17:
                    divergenceJac<17, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 18:
                    divergenceJac<18, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 19:
                    divergenceJac<19, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 20:
                    divergenceJac<20, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 21:
                    divergenceJac<21, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 22:
                    divergenceJac<22, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 23:
                    divergenceJac<23, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 24:
                    divergenceJac<24, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 25:
                    divergenceJac<25, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 26:
                    divergenceJac<26, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 27:
                    divergenceJac<27, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 28:
                    divergenceJac<28, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 29:
                    divergenceJac<29, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 30:
                    divergenceJac<30, false, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                default:
                    mexErrMsgTxt("%d derivatives not supported by backend.");
            }
        }
    } else if (nrows == nf) {
        // ColMajor
        if (has_accumulation) {
            // mexPrintf("Column major, with accumulation: %d by %d.\n", nc, nf);
            switch (m) {
                case 1:
                    divergenceJac<1, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 2:
                    divergenceJac<2, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 3:
                    divergenceJac<3, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 4:
                    divergenceJac<4, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 5:
                    divergenceJac<5, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 6:
                    divergenceJac<6, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 7:
                    divergenceJac<7, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 8:
                    divergenceJac<8, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 9:
                    divergenceJac<9, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 10:
                    divergenceJac<10, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 11:
                    divergenceJac<11, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 12:
                    divergenceJac<12, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 13:
                    divergenceJac<13, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 14:
                    divergenceJac<14, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 15:
                    divergenceJac<15, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 16:
                    divergenceJac<16, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 17:
                    divergenceJac<17, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 18:
                    divergenceJac<18, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 19:
                    divergenceJac<19, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 20:
                    divergenceJac<20, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 21:
                    divergenceJac<21, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 22:
                    divergenceJac<22, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 23:
                    divergenceJac<23, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 24:
                    divergenceJac<24, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 25:
                    divergenceJac<25, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 26:
                    divergenceJac<26, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 27:
                    divergenceJac<27, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 28:
                    divergenceJac<28, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 29:
                    divergenceJac<29, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 30:
                    divergenceJac<30, true, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                default:
                    mexErrMsgTxt("%d derivatives not supported by backend.");
            }
        } else {
            // mexPrintf("Column major, without accumulation: %d by %d.\n", nc, nf);
            switch (m) {
                case 1:
                    divergenceJac<1, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 2:
                    divergenceJac<2, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 3:
                    divergenceJac<3, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 4:
                    divergenceJac<4, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 5:
                    divergenceJac<5, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 6:
                    divergenceJac<6, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 7:
                    divergenceJac<7, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 8:
                    divergenceJac<8, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 9:
                    divergenceJac<9, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 10:
                    divergenceJac<10, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 11:
                    divergenceJac<11, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 12:
                    divergenceJac<12, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 13:
                    divergenceJac<13, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 14:
                    divergenceJac<14, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 15:
                    divergenceJac<15, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 16:
                    divergenceJac<16, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 17:
                    divergenceJac<17, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 18:
                    divergenceJac<18, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 19:
                    divergenceJac<19, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 20:
                    divergenceJac<20, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 21:
                    divergenceJac<21, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 22:
                    divergenceJac<22, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 23:
                    divergenceJac<23, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 24:
                    divergenceJac<24, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 25:
                    divergenceJac<25, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 26:
                    divergenceJac<26, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 27:
                    divergenceJac<27, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 28:
                    divergenceJac<28, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 29:
                    divergenceJac<29, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                case 30:
                    divergenceJac<30, false, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                        accumulation, diagonal, pr, ir, jc);
                    break;
                default:
                    mexErrMsgTxt("%d derivatives not supported by backend.");
            }
        }
    } else {
    mexErrMsgTxt("Dimensions of face diagonal matrix does not fit either RowMajor or ColMajor.");
    }
}


/*
if (c < 0) {
    // Low entry, corresponding to N(f, 1)
    for (int der = 0; der < m; der++) {
        int sparse_offset = jc[cell + nc * der];
        double v = diagonal[f * 2 * m + der];
        // Set row entry
        ir[sparse_offset + i] = -c;
        pr[sparse_offset + i] = -v;
        // Set corresponding diagonal entry
        pr[sparse_offset + diag] += v;
    }

}
else {
    // High entry, corresponding to N(f, 2)
    for (int der = 0; der < m; der++) {
        int sparse_offset = jc[cell + nc * der];
        double v = diagonal[(f * 2 + 1)* m + der];
        // Set row entry
        ir[sparse_offset + i] = c;
        pr[sparse_offset + i] = v;
        // Set corresponding diagonal entry
        pr[sparse_offset + diag] -= v;
    }
}
*/