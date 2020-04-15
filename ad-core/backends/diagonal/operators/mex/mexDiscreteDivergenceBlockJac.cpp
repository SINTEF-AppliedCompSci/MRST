//
// include necessary system headers
//
#include <cmath>
#include "mex.h"  

#ifndef HAVE_OCTAVE
#include "matrix.h"
#endif

#include <array>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <iostream>
#include <chrono>
#include <vector>

#define DIVBLOCK(bz) divergenceJacBlock<bz>(nf, nc, n_jacs, facePos, faces, cells, cells_ix, cellDiagonals, faceDiagonals, pr, ir, jc)



// template <int m, bool has_accumulation, bool colMajor>
template <int nder>
void divergenceJacBlock(const int nf, const int nc, const int njac,
    const double* facePos, const double* faces,
    const double* cells, const double* cells_ix,
    std::vector<double*> &cellDiagonals, std::vector<double*> &faceDiagonals,
    double* pr, size_t * ir, size_t * jc) {
    int mv = facePos[nc];
    int jac_val_width = nder * njac;
#pragma omp parallel for schedule(dynamic, 512)
    for (int cell = 0; cell < nc; cell++) {
        // Each cell has number of connections equal to the number of half-
        // faces for that cell plus itself multiplied by the block size
        int f_offset = facePos[cell];
        int n_local_hf = facePos[cell + 1] - f_offset;
        int row_width = n_local_hf + 1; // diagonal + faces
        int diag = cells_ix[cell];
        
        int row_start = f_offset + cell;
        int row_end = row_start + row_width;
        // Base offset taking into account how far we have come
        jc[cell + 1] = row_end;
        // Set diagonal entries
        int global_diag_pos = row_start + diag;
        ir[global_diag_pos] = cell;

        //mexPrintf("Row start %d -> global_diag_pos = %d\n", row_start, global_diag_pos);
        for (int jacNo = 0; jacNo < njac; jacNo++) {
            // cell_offset

            int start = global_diag_pos * nder * njac  // row at diagonal * block size * number of block-sized Jacobians
                                        + nder * jacNo;  // Offset 
            // mexPrintf("Cell %d Jac %d Writing indices start %d -> %d\n", cell, jacNo, start, start + nder - 1);
            double* accumulation = cellDiagonals[jacNo];
            for (int der = 0; der < nder; der++) {
                pr[start + der] = accumulation[cell*nder + der];
                
            }
            
        }

        for (int fl = 0; fl < n_local_hf; fl++) {
            // Loop over entire column
            // Diagonal entry can be skipped - we handle this later
            // Check if we have passed diagonal entry
            int passed = (double)(fl >= diag);
            int f_i = fl + passed;
            // Global face index
            int f = faces[f_offset + fl];
            // Global cell index
            int c = cells[f_offset + fl];
            int s;
            if (c > 0) {
                s = 1;
            } else {
                s = -1;
            }
            // Iterate over derivatives
            int sparse_mult = mv + nc;
            ir[row_start + f_i] = abs(c);

            int face_start = (row_start + f_i) * nder * njac;
            int diag_start = global_diag_pos * nder * njac;
            for (int jacNo = 0; jacNo < njac; jacNo++) {
                int diff = nder * jacNo;
                int start = face_start + diff;  // Offset 
                int diag  = diag_start + diff;  // Offset 

                double* facejac = faceDiagonals[jacNo];
                if (c < 0) {
                    // Low entry, corresponding to N(f, 1)
                    for (int der = 0; der < nder; der++) {
                        pr[start + der] = facejac[(2 * f + 1) * nder + der];
                        // Set corresponding diagonal entry
                        pr[diag + der] += facejac[(2 * f) * nder + der];
                    }
                }
                else {
                    for (int der = 0; der < nder; der++) {
                        pr[start + der] = -facejac[(2 * f) * nder + der];
                        // Set corresponding diagonal entry
                        pr[diag + der] -= facejac[(2 * f + 1) * nder + der];
                    }
                }
            }
        }
    }
}

/* MEX gateway */

void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[])

{
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
    }
    else if (nrhs != 9) {
        mexErrMsgTxt("9 input arguments required.");
    }
    else if (nlhs != 5 && nlhs != 0) {
        mexErrMsgTxt("Function must either produce five dense outputs (I, J, V, n, m) or no output.");
    }
    bool output_sparse = nlhs < 2;
    const mxArray* accJac = prhs[0];
    const mxArray* faceJac = prhs[1]; // Cell arrays
    double* N = mxGetPr(prhs[2]);
    double* facePos = mxGetPr(prhs[3]);
    double* faces = mxGetPr(prhs[4]);
    double* cells = mxGetPr(prhs[5]);
    double* cells_ix = mxGetPr(prhs[6]);
    // accjac, facejac, N, facePos, faces, cells, cells_ix, nder, is_row_major
    int nder = mxGetScalar(prhs[7]);
    bool rowMajor = mxGetScalar(prhs[8]);
    size_t n_jacs = mxGetNumberOfElements(faceJac);


    
    // mexPrintf("Cell %d: %d -> %d\n", cell, start + der, cell * nder + der);
    // std::vector<bool> hasAccumulation;
    std::vector<double*> cellDiagonals;
    std::vector<double*> faceDiagonals;


    for (int jacNo = 0; jacNo < n_jacs; jacNo++) {
        /*

        mxArray* tmp2 = mxGetCell(&faceJacs[jacNo], jacNo);
        cellDiagonals[jacNo] = (double*)mxGetPr(tmp1);
        faceDiagonals[jacNo] = (double*)mxGetPr(tmp2);
        */
        cellDiagonals.push_back(mxGetPr(mxGetCell(prhs[0], jacNo)));
        faceDiagonals.push_back(mxGetPr(mxGetCell(prhs[1], jacNo)));
    }
    /*
    double* accumulation = mxGetPr(accJac);
    double* diagonal = mxGetPr(faceJac);
    */

    int nf = mxGetM(prhs[2]);
    int nc = mxGetM(prhs[3]) - 1;
    int nhf = mxGetM(prhs[4]);

    // mexPrintf("%d Jacobians, %d derivatives, %d faces and %d cells\n", n_jacs, nder, nf, nc);

    // int n_acc = mxGetNumberOfElements(accJac);

    // bool has_accumulation = n_acc > 0;

    int m = nder;
    // mexPrintf("Matrix has dimensions %d by %d. There are %d faces and %d cells (%d interfaces, %d derivatives)\n", nrows, ncols, nf, nc, nf_in, m);

    // Each cell has one self-connection plus the number of half-faces
    mwSize nzmax = (facePos[nc] + nc);
    // printf("%d cells %d faces, %d half-faces and %d derivatives \n", nc, nf, nhf, m);
    // Row indices, zero-indexed (direct entries)
    size_t* ir;
    // Column indices, zero-indexed, offset encoded of length m*nc + 1
    size_t* jc;
    // Entries
    double* pr;

    plhs[0] = mxCreateNumericMatrix(nzmax, 1, mxUINT64_CLASS, mxREAL); // We have no way of allocating mwIndex (size_t). So we hope for the best and allocate uint64...
    plhs[1] = mxCreateNumericMatrix(nc + 1, 1, mxUINT64_CLASS, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(n_jacs * nder, nzmax, mxREAL);
    ir = (mwIndex*)mxGetData(plhs[0]);
    jc = (mwIndex*)mxGetData(plhs[1]);
    // Entries
    pr = mxGetPr(plhs[2]);
    plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);

    double* three = mxGetPr(plhs[3]);
    three[0] = nc;
    double* four = mxGetPr(plhs[4]);
    four[0] = nc;

    switch (nder) {
        case 1:
            DIVBLOCK(1);
            break;
        case 2:
            DIVBLOCK(2);
            break;
        case 3:
            DIVBLOCK(3);
            break;
        case 4:
            DIVBLOCK(4);
            break;
        case 5:
            DIVBLOCK(5);
            break;
        case 6:
            DIVBLOCK(6);
            break;
        case 7:
            DIVBLOCK(7);
            break;
        case 8:
            DIVBLOCK(8);
            break;
        case 9:
            DIVBLOCK(9);
            break;
        case 10:
            DIVBLOCK(10);
            break;
        case 11:
            DIVBLOCK(11);
            break;
        case 12:
            DIVBLOCK(12);
            break;
        default:
            mexErrMsgTxt("Block size not supported!");
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