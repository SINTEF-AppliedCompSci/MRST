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
#define TIME_NOW std::chrono::high_resolution_clock::now

template <int m>
void copyElements(const int offset_face, const int offset_cell, const double * celldata, double * facedata) {
    for (int der = 0; der < m; der++) {
        facedata[offset_face + der] = celldata[offset_cell + der];
    }
}

template <int m>
void upwindJacColMajor(const int nf, const int nc, const mxLogical * flag, const double * diagonal, const double * N, double * result){
    #pragma omp parallel for schedule(static)
    for(int i=0;i<nf;i++){
        int cell_inx, offset;
        if (flag[i]) {
            cell_inx = N[i] - 1;
            offset = 0;
        }
        else {
            cell_inx = N[i + nf] - 1;
            offset = m * nf;
        }
        for (int j = 0; j < m; j++) {
            result[j * nf + offset + i] = diagonal[nc * j + cell_inx];
        }
    }
    return;
}

template <int m>
void upwindJacRowMajor(const int nf, const int nc, const mxLogical* flag, const double* diagonal, const double* N, double* result) {
    /*
#pragma omp parallel for schedule(static)
    for (int face = 0; face < nf; face++) {
        int cell_inx, offset;
        if (flag[face]) {
            cell_inx = N[face] - 1;
            offset = 0;
        }
        else {
            cell_inx = N[face + nf] - 1;
            offset = m;
        }
        for (int der = 0; der < m; der++) {
            result[face * 2 * m + der + offset] = diagonal[m * cell_inx + der];
        }
    }
    */
    /**/

    /*
#pragma omp parallel
    {
    #pragma omp for schedule(static) nowait
            for (int face = 0; face < nf; face++) {
                if (!flag[face]) {
                    continue;
                }
                int cell = N[face] - 1;
                for (int der = 0; der < m; der++) {
                    result[face * 2 * m + der] = diagonal[m * cell + der];
                }
                // memcpy(&result[face * 2 * m + m * do_up], &diagonal[m * cell_inx], m * sizeof(double));
            }
    #pragma omp for schedule(static)
            for (int face = 0; face < nf; face++) {
                if (flag[face]) {
                    continue;
                }
                int cell = N[face + nf] - 1;
                for (int der = 0; der < m; der++) {
                    result[(face * 2 + 1 )* m + der] = diagonal[m * cell + der];
                }
                // memcpy(&result[face * 2 * m + m * do_up], &diagonal[m * cell_inx], m * sizeof(double));
            }

    }

    */
    
#pragma omp parallel for schedule(static)
    for (int face = 0; face < nf; face++) {
        bool do_up = !flag[face];
        int cell_inx = N[face + do_up * nf] - 1;
        copyElements<m>(face * 2 * m + m * do_up, m * cell_inx, diagonal, result);
        /*
        for (int der = 0; der < m; der++) {
            result[face * 2 * m + der + m * do_up] = diagonal[m * cell_inx + der];
        }
        */
        // memcpy(&result[face * 2 * m + m * do_up], &diagonal[m * cell_inx], m * sizeof(double));
    }
}

/* MEX gateway */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
    // INPUTS:
    //  - cell_diagonal<double> [nc x m] if column major or [m x nc] if row major)
    //  - N<double>             [nf x 2]
    //  - flag<bool>            [nf x 1]
    //  - nc<double>            [scalar]
    // OUTPUT:
    //  - face_diagonal<double> [nf x 2*m] if column major or [2*m x nf] if row major
    // auto t0 = TIME_NOW();
    if (nrhs == 0) {
        if (nlhs > 0) {
            mexErrMsgTxt("Cannot give outputs with no inputs.");
        }
        // We are being called through compilation testing. Just do nothing. 
        // If the binary was actually called, we are good to go.
        return;
    } else if (nrhs != 5) {
        mexErrMsgTxt("5 input arguments required: Diagonal, N, flag, number of cells, rowMajor indicator");
    }
    else if (nlhs > 1) {
        mexErrMsgTxt("Too many outputs requested. Function has a single output argument.");
    }
    double* diagonal = mxGetPr(prhs[0]);
    double* N = mxGetPr(prhs[1]);
    mxLogical* flag = mxGetLogicals(prhs[2]);
    int nc = mxGetScalar(prhs[3]);
    int nf = mxGetM(prhs[1]);

    bool rowMajor = mxGetScalar(prhs[4]);

    // Dimensions of diagonals - figure out if we want row or column major solver
    int nrows = mxGetM(prhs[0]);
    int ncols = mxGetN(prhs[0]);
    // auto t1 = TIME_NOW();
    // mexPrintf("Matrix has dimensions %d by %d. There are %d faces and %d cells\n", nrows, ncols, nf, nc);
    if (!rowMajor) {
        // ColMajor
        int m = ncols;
        plhs[0] = mxCreateDoubleMatrix(nf, 2 * m, mxREAL);
        double* result = mxGetPr(plhs[0]);
        switch (m) {
        case 1:
            upwindJacColMajor<1>(nf, nc, flag, diagonal, N, result);
            break;
        case 2:
            upwindJacColMajor<2>(nf, nc, flag, diagonal, N, result);
            break;
        case 3:
            upwindJacColMajor<3>(nf, nc, flag, diagonal, N, result);
            break;
        case 4:
            upwindJacColMajor<4>(nf, nc, flag, diagonal, N, result);
            break;
        case 5:
            upwindJacColMajor<5>(nf, nc, flag, diagonal, N, result);
            break;
        case 6:
            upwindJacColMajor<6>(nf, nc, flag, diagonal, N, result);
            break;
        case 7:
            upwindJacColMajor<7>(nf, nc, flag, diagonal, N, result);
            break;
        case 8:
            upwindJacColMajor<8>(nf, nc, flag, diagonal, N, result);
            break;
        case 9:
            upwindJacColMajor<9>(nf, nc, flag, diagonal, N, result);
            break;
        case 10:
            upwindJacColMajor<10>(nf, nc, flag, diagonal, N, result);
            break;
        case 11:
            upwindJacColMajor<11>(nf, nc, flag, diagonal, N, result);
            break;
        case 12:
            upwindJacColMajor<12>(nf, nc, flag, diagonal, N, result);
            break;
        case 13:
            upwindJacColMajor<13>(nf, nc, flag, diagonal, N, result);
            break;
        case 14:
            upwindJacColMajor<14>(nf, nc, flag, diagonal, N, result);
            break;
        case 15:
            upwindJacColMajor<15>(nf, nc, flag, diagonal, N, result);
            break;
        case 16:
            upwindJacColMajor<16>(nf, nc, flag, diagonal, N, result);
            break;
        case 17:
            upwindJacColMajor<17>(nf, nc, flag, diagonal, N, result);
            break;
        case 18:
            upwindJacColMajor<18>(nf, nc, flag, diagonal, N, result);
            break;
        case 19:
            upwindJacColMajor<19>(nf, nc, flag, diagonal, N, result);
            break;
        case 20:
            upwindJacColMajor<20>(nf, nc, flag, diagonal, N, result);
            break;
        case 21:
            upwindJacColMajor<21>(nf, nc, flag, diagonal, N, result);
            break;
        case 22:
            upwindJacColMajor<22>(nf, nc, flag, diagonal, N, result);
            break;
        case 23:
            upwindJacColMajor<23>(nf, nc, flag, diagonal, N, result);
            break;
        case 24:
            upwindJacColMajor<24>(nf, nc, flag, diagonal, N, result);
            break;
        case 25:
            upwindJacColMajor<25>(nf, nc, flag, diagonal, N, result);
            break;
        case 26:
            upwindJacColMajor<26>(nf, nc, flag, diagonal, N, result);
            break;
        case 27:
            upwindJacColMajor<27>(nf, nc, flag, diagonal, N, result);
            break;
        case 28:
            upwindJacColMajor<28>(nf, nc, flag, diagonal, N, result);
            break;
        case 29:
            upwindJacColMajor<29>(nf, nc, flag, diagonal, N, result);
            break;
        case 30:
            upwindJacColMajor<30>(nf, nc, flag, diagonal, N, result);
            break;
        default:
            mexErrMsgTxt("%d derivatives not supported by backend.");
        }
    }
    else {
        // RowMajor
        int m = nrows;
        plhs[0] = mxCreateDoubleMatrix(2 * m, nf, mxREAL);
        double* result = mxGetPr(plhs[0]);
        auto t2 = TIME_NOW();
        switch (m) {
        case 1:
            upwindJacRowMajor<1>(nf, nc, flag, diagonal, N, result);
            break;
        case 2:
            upwindJacRowMajor<2>(nf, nc, flag, diagonal, N, result);
            break;
        case 3:
            upwindJacRowMajor<3>(nf, nc, flag, diagonal, N, result);
            break;
        case 4:
            upwindJacRowMajor<4>(nf, nc, flag, diagonal, N, result);
            break;
        case 5:
            upwindJacRowMajor<5>(nf, nc, flag, diagonal, N, result);
            break;
        case 6:
            upwindJacRowMajor<6>(nf, nc, flag, diagonal, N, result);
            break;
        case 7:
            upwindJacRowMajor<7>(nf, nc, flag, diagonal, N, result);
            break;
        case 8:
            upwindJacRowMajor<8>(nf, nc, flag, diagonal, N, result);
            break;
        case 9:
            upwindJacRowMajor<9>(nf, nc, flag, diagonal, N, result);
            break;
        case 10:
            upwindJacRowMajor<10>(nf, nc, flag, diagonal, N, result);
            break;
        case 11:
            upwindJacRowMajor<11>(nf, nc, flag, diagonal, N, result);
            break;
        case 12:
            upwindJacRowMajor<12>(nf, nc, flag, diagonal, N, result);
            break;
        case 13:
            upwindJacRowMajor<13>(nf, nc, flag, diagonal, N, result);
            break;
        case 14:
            upwindJacRowMajor<14>(nf, nc, flag, diagonal, N, result);
            break;
        case 15:
            upwindJacRowMajor<15>(nf, nc, flag, diagonal, N, result);
            break;
        case 16:
            upwindJacRowMajor<16>(nf, nc, flag, diagonal, N, result);
            break;
        case 17:
            upwindJacRowMajor<17>(nf, nc, flag, diagonal, N, result);
            break;
        case 18:
            upwindJacRowMajor<18>(nf, nc, flag, diagonal, N, result);
            break;
        case 19:
            upwindJacRowMajor<19>(nf, nc, flag, diagonal, N, result);
            break;
        case 20:
            upwindJacRowMajor<20>(nf, nc, flag, diagonal, N, result);
            break;
        case 21:
            upwindJacRowMajor<21>(nf, nc, flag, diagonal, N, result);
            break;
        case 22:
            upwindJacRowMajor<22>(nf, nc, flag, diagonal, N, result);
            break;
        case 23:
            upwindJacRowMajor<23>(nf, nc, flag, diagonal, N, result);
            break;
        case 24:
            upwindJacRowMajor<24>(nf, nc, flag, diagonal, N, result);
            break;
        case 25:
            upwindJacRowMajor<25>(nf, nc, flag, diagonal, N, result);
            break;
        case 26:
            upwindJacRowMajor<26>(nf, nc, flag, diagonal, N, result);
            break;
        case 27:
            upwindJacRowMajor<27>(nf, nc, flag, diagonal, N, result);
            break;
        case 28:
            upwindJacRowMajor<28>(nf, nc, flag, diagonal, N, result);
            break;
        case 29:
            upwindJacRowMajor<29>(nf, nc, flag, diagonal, N, result);
            break;
        case 30:
            upwindJacRowMajor<30>(nf, nc, flag, diagonal, N, result);
            break;
        default:
            mexErrMsgTxt("%d derivatives not supported by backend.");
        }
        /*
        if (1) {
            auto t3 = TIME_NOW();
            auto p_input = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
            auto p_alloc = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
            auto p_compute = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count();
            mexPrintf("Input %d us, alloc %d us, compute %d us\n", p_input, p_alloc, p_compute);
        }
        */
    }
    return;
}


