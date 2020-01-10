//
// include necessary system headers
//
#include <cmath>
#include <mex.h>
#include <array>
#include <omp.h>
#include <iostream>
#include <chrono>

template <int m, bool has_accumulation>
void divergenceJac(int nf, int nc,
                   const double * N, const double * facePos, const double * faces, 
                   const double * cells, const double * cells_ix,
                   const double * accumulation, const double * diagonal,
                   double * pr, mwIndex * ir, mwIndex * jc){
    int mv = facePos[nc];
    #pragma omp parallel for
    for(int index = 0; index < m*nc; index++){
        int col = index % nc;
        int der = index / nc;
        // Each cell has number of connections equal to the number of half-
        // faces for that cell plus itself multiplied by the block size
        int prev = facePos[col];
        int n_local_hf = facePos[col+1] - prev;
        int ix = col + der*nc;
        // Base offset taking into account how far we have come
        int base = der*(mv + nc) + (prev + col);
        jc[ix+1] = base + n_local_hf+1;
        // Set diagonal entries
        int dpos = base + cells_ix[col];
        ir[dpos] = col;
        if(has_accumulation){
            pr[dpos] = accumulation[index];
        }else{
            pr[dpos] = 0.0;
        }
    }
    // Loop over cells and assemble Jacobian
    #pragma omp parallel for
    for(int outer_ix = 0; outer_ix < nc*m; outer_ix++){
        int cell = outer_ix % nc;
        int der = outer_ix / nc;
        int f_offset = facePos[cell];
        // Number of local faces
        int nlf = facePos[cell+1] - f_offset;
        int diag = cells_ix[cell];

        int start = cell + nc*der;
        int sparse_offset = jc[start];
        int width = jc[start+1]-jc[start];
        for(int i=0; i<width; i++){
            // Loop over entire column
            // Local face index
            int fl = i % (nlf+1);
            // Diagonal entry can be skipped - we handle this later
            if(fl == diag){
                continue;
            }
            // Check if we have passed diagonal entry
            int passed = (double)(fl > diag);
            // Global face index
            int f = faces[f_offset + fl - passed];
            // Global cell index
            int c = cells[f_offset + fl- passed];
            double other_val, v;
            int diag_index, sgn;
            if(c >= 0){
                // High entry, corresponding to N(f, 2)
                v = diagonal[der*2*nf + f + nf];
            }else{
                // Low entry, corresponding to N(f, 1)
                v = -diagonal[der*2*nf + f];
            }
            // Set row entry
            ir[sparse_offset + i] = abs(c);
            pr[sparse_offset + i] = v;
            // Set corresponding diagonal entry
            #pragma omp atomic
            pr[sparse_offset + diag] -= v;
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
    // diagonal (nf x m)
    // N (nf x 2)
    // facePos (nc+1 x 1)
    // faces (length facePos(end)-1)
    // Out: I, J, V
    if (nrhs != 7) { 
	    mexErrMsgTxt("7 input arguments required."); 
    } else if (nlhs > 1) {
	    mexErrMsgTxt("Wrong number of output arguments."); 
    }
    const mxArray * accJac = prhs[0];
    const mxArray * diagJac = prhs[1];
    double * N = mxGetPr(prhs[2]);
    double * facePos = mxGetPr(prhs[3]);
    double * faces = mxGetPr(prhs[4]);
    double * cells = mxGetPr(prhs[5]);
    double * cells_ix = mxGetPr(prhs[6]);

    /* Diagonal of accumulation term */
    double * accumulation = mxGetPr(accJac);
    /* Diagonal of fluxes */
    double * diagonal = mxGetPr(prhs[1]);

    int m = mxGetN(diagJac);
    int nf = mxGetM(prhs[2]);
    int nc = mxGetM(prhs[3])-1;
    int nhf = mxGetM(prhs[4]);
    int n_acc = mxGetNumberOfElements(accJac);
    
    bool has_accumulation = n_acc > 0;
    
    if(has_accumulation && n_acc != m*nc){
        mexErrMsgTxt("Accumulation term diagonal was provided, but dimensions are incorrect.");
    }
    {
        int nf_in = mxGetM(diagJac);
        if(nf_in != 2*nf){
            mexErrMsgTxt("Face inputs do not match.");
        }
    }
    // Each cell has one self-connection plus the number of half-faces, multiplied by block size
    mwSize nzmax = (facePos[nc] + nc)*m;
    // printf("%d cells %d faces, %d half-faces and %d derivatives \n", nc, nf, nhf, m);
    plhs[0] = mxCreateSparse(nc, nc*m, nzmax, mxREAL);
    
    // Entries
    double * pr  = mxGetPr(plhs[0]);
    // Row indices, zero-indexed (direct entries)
    mwIndex * ir = mxGetIr(plhs[0]);
    // Column indices, zero-indexed, offset encoded of length m*nc + 1
    mwIndex * jc = mxGetJc(plhs[0]);
    if (has_accumulation){
        switch (m){
            case 1:
                divergenceJac<1, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 2:
                divergenceJac<2, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 3:
                divergenceJac<3, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 4:
                divergenceJac<4, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 5:
                divergenceJac<5, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 6:
                divergenceJac<6, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 7:
                divergenceJac<7, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 8:
                divergenceJac<8, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 9:
                divergenceJac<9, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 10:
                divergenceJac<10, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 11:
                divergenceJac<11, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 12:
                divergenceJac<12, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 13:
                divergenceJac<13, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 14:
                divergenceJac<14, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 15:
                divergenceJac<15, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 16:
                divergenceJac<16, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 17:
                divergenceJac<17, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 18:
                divergenceJac<18, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 19:
                divergenceJac<19, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 20:
                divergenceJac<20, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 21:
                divergenceJac<21, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 22:
                divergenceJac<22, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 23:
                divergenceJac<23, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 24:
                divergenceJac<24, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 25:
                divergenceJac<25, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 26:
                divergenceJac<26, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 27:
                divergenceJac<27, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 28:
                divergenceJac<28, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 29:
                divergenceJac<29, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 30:
                divergenceJac<30, true>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            default:
                mexErrMsgTxt("%d derivatives not supported by backend."); 
        }
    }else{
        switch (m){
            case 1:
                divergenceJac<1, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 2:
                divergenceJac<2, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 3:
                divergenceJac<3, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 4:
                divergenceJac<4, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 5:
                divergenceJac<5, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 6:
                divergenceJac<6, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 7:
                divergenceJac<7, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 8:
                divergenceJac<8, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 9:
                divergenceJac<9, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 10:
                divergenceJac<10, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 11:
                divergenceJac<11, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 12:
                divergenceJac<12, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 13:
                divergenceJac<13, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 14:
                divergenceJac<14, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 15:
                divergenceJac<15, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 16:
                divergenceJac<16, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 17:
                divergenceJac<17, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 18:
                divergenceJac<18, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 19:
                divergenceJac<19, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 20:
                divergenceJac<20, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 21:
                divergenceJac<21, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 22:
                divergenceJac<22, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 23:
                divergenceJac<23, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 24:
                divergenceJac<24, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 25:
                divergenceJac<25, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 26:
                divergenceJac<26, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 27:
                divergenceJac<27, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 28:
                divergenceJac<28, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 29:
                divergenceJac<29, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            case 30:
                divergenceJac<30, false>(nf, nc, N, facePos, faces, cells, cells_ix,
                            accumulation, diagonal, pr, ir, jc);
                break;
            default:
                mexErrMsgTxt("%d derivatives not supported by backend."); 
        }
    }
}


