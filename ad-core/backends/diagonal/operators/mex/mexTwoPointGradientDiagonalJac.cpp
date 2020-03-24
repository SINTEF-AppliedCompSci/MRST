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
/* Templated function for main operation */

template <int m>
void gradientJacColMajor(const int nf, const int nc, const double * diagonal, const double * N, double * result){
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < nf; i++) {
        int left = N[i] - 1;
        int right = N[i + nf] - 1;
        for (int j = 0; j < m; j++) {
            result[j * nf + i]      = -diagonal[nc * j + left];
            result[j * nf + i + m*nf] =  diagonal[nc * j + right];
        }
    }
    return;
}

template <int m>
void gradientJacRowMajor(const int nf, const int nc, const double* diagonal, const double* N, double* result) {
    /*
#pragma omp parallel
    {
#pragma omp for nowait schedule(static)
        for (int i = 0; i < nf; i++) {
            int cell = N[i] - 1;
            for (int j = 0; j < m; j++) {
                result[i * 2 * m + j] = -diagonal[m * cell + j];
            }
        }
#pragma omp for nowait schedule(static)
        for (int i = 0; i < nf; i++) {
            int cell = N[i + nf] - 1;
            for (int j = 0; j < m; j++) {
                result[(i * 2 + 1) * m + j] = diagonal[m * cell + j];
            }
        }
    }
    return;
    */
#pragma omp parallel for schedule(static)
    for (int i = 0; i < nf; i++) {
        int left = N[i] - 1;
        int right = N[i + nf] - 1;
        for (int j = 0; j < m; j++) {
            result[i * 2 * m + j] = -diagonal[m * left + j];
            result[i * 2 * m + j + m] = diagonal[m * right + j];
        }
    }
    return;
}
/* MEX gateway */
void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
    // INPUTS:
    //  - cell_diagonal<double> [nc x m] if column major or [m x nc] if row major)
    //  - N<double>             [nf x 2]
    //  - nc<double>            [scalar]
    // OUTPUT:
    //  - face_diagonal<double> [nf x 2*m] if column major or [2*m x nf] if row major
    if (nrhs == 0) {
        if (nlhs > 0) {
            mexErrMsgTxt("Cannot give outputs with no inputs.");
        }
        // We are being called through compilation testing. Just do nothing. 
        // If the binary was actually called, we are good to go.
        return;
    } else if (nrhs != 4) {
        mexErrMsgTxt("4 input arguments required: Diagonal, N, number of cells and rowMajor bool");
    } else if (nlhs > 1) {
        mexErrMsgTxt("Too many outputs requested. Function has a single output argument.");
    }
    double* diagonal = mxGetPr(prhs[0]);
    double* N = mxGetPr(prhs[1]);
    int nc = mxGetScalar(prhs[2]);
    int nf = mxGetM(prhs[1]);
    bool rowMajor = mxGetScalar(prhs[3]);

    // Dimensions of diagonals - figure out if we want row or column major solver
    int nrows = mxGetM(prhs[0]);
    int ncols = mxGetN(prhs[0]);

    // mexPrintf("Matrix has dimensions %d by %d. There are %d faces and %d cells\n", nrows, ncols, nf, nc);
    if (rowMajor) {
        // RowMajor
        int m = nrows;
        plhs[0] = mxCreateDoubleMatrix(2 * m, nf, mxREAL);
        double* result = mxGetPr(plhs[0]);
        switch (m) {
        case 1:
            gradientJacRowMajor<1>(nf, nc, diagonal, N, result);
            break;
        case 2:
            gradientJacRowMajor<2>(nf, nc, diagonal, N, result);
            break;
        case 3:
            gradientJacRowMajor<3>(nf, nc, diagonal, N, result);
            break;
        case 4:
            gradientJacRowMajor<4>(nf, nc, diagonal, N, result);
            break;
        case 5:
            gradientJacRowMajor<5>(nf, nc, diagonal, N, result);
            break;
        case 6:
            gradientJacRowMajor<6>(nf, nc, diagonal, N, result);
            break;
        case 7:
            gradientJacRowMajor<7>(nf, nc, diagonal, N, result);
            break;
        case 8:
            gradientJacRowMajor<8>(nf, nc, diagonal, N, result);
            break;
        case 9:
            gradientJacRowMajor<9>(nf, nc, diagonal, N, result);
            break;
        case 10:
            gradientJacRowMajor<10>(nf, nc, diagonal, N, result);
            break;
        case 11:
            gradientJacRowMajor<11>(nf, nc, diagonal, N, result);
            break;
        case 12:
            gradientJacRowMajor<12>(nf, nc, diagonal, N, result);
            break;
        case 13:
            gradientJacRowMajor<13>(nf, nc, diagonal, N, result);
            break;
        case 14:
            gradientJacRowMajor<14>(nf, nc, diagonal, N, result);
            break;
        case 15:
            gradientJacRowMajor<15>(nf, nc, diagonal, N, result);
            break;
        case 16:
            gradientJacRowMajor<16>(nf, nc, diagonal, N, result);
            break;
        case 17:
            gradientJacRowMajor<17>(nf, nc, diagonal, N, result);
            break;
        case 18:
            gradientJacRowMajor<18>(nf, nc, diagonal, N, result);
            break;
        case 19:
            gradientJacRowMajor<19>(nf, nc, diagonal, N, result);
            break;
        case 20:
            gradientJacRowMajor<20>(nf, nc, diagonal, N, result);
            break;
        case 21:
            gradientJacRowMajor<21>(nf, nc, diagonal, N, result);
            break;
        case 22:
            gradientJacRowMajor<22>(nf, nc, diagonal, N, result);
            break;
        case 23:
            gradientJacRowMajor<23>(nf, nc, diagonal, N, result);
            break;
        case 24:
            gradientJacRowMajor<24>(nf, nc, diagonal, N, result);
            break;
        case 25:
            gradientJacRowMajor<25>(nf, nc, diagonal, N, result);
            break;
        case 26:
            gradientJacRowMajor<26>(nf, nc, diagonal, N, result);
            break;
        case 27:
            gradientJacRowMajor<27>(nf, nc, diagonal, N, result);
            break;
        case 28:
            gradientJacRowMajor<28>(nf, nc, diagonal, N, result);
            break;
        case 29:
            gradientJacRowMajor<29>(nf, nc, diagonal, N, result);
            break;
        case 30:
            gradientJacRowMajor<30>(nf, nc, diagonal, N, result);
            break;
        default:
            mexErrMsgTxt("%d derivatives not supported by backend.");
        }
    } else {
        // ColMajor
        int m = ncols;
        plhs[0] = mxCreateDoubleMatrix(nf, 2 * m, mxREAL);
        double* result = mxGetPr(plhs[0]);
        switch (m) {
        case 1:
            gradientJacColMajor<1>(nf, nc, diagonal, N, result);
            break;
        case 2:
            gradientJacColMajor<2>(nf, nc, diagonal, N, result);
            break;
        case 3:
            gradientJacColMajor<3>(nf, nc, diagonal, N, result);
            break;
        case 4:
            gradientJacColMajor<4>(nf, nc, diagonal, N, result);
            break;
        case 5:
            gradientJacColMajor<5>(nf, nc, diagonal, N, result);
            break;
        case 6:
            gradientJacColMajor<6>(nf, nc, diagonal, N, result);
            break;
        case 7:
            gradientJacColMajor<7>(nf, nc, diagonal, N, result);
            break;
        case 8:
            gradientJacColMajor<8>(nf, nc, diagonal, N, result);
            break;
        case 9:
            gradientJacColMajor<9>(nf, nc, diagonal, N, result);
            break;
        case 10:
            gradientJacColMajor<10>(nf, nc, diagonal, N, result);
            break;
        case 11:
            gradientJacColMajor<11>(nf, nc, diagonal, N, result);
            break;
        case 12:
            gradientJacColMajor<12>(nf, nc, diagonal, N, result);
            break;
        case 13:
            gradientJacColMajor<13>(nf, nc, diagonal, N, result);
            break;
        case 14:
            gradientJacColMajor<14>(nf, nc, diagonal, N, result);
            break;
        case 15:
            gradientJacColMajor<15>(nf, nc, diagonal, N, result);
            break;
        case 16:
            gradientJacColMajor<16>(nf, nc, diagonal, N, result);
            break;
        case 17:
            gradientJacColMajor<17>(nf, nc, diagonal, N, result);
            break;
        case 18:
            gradientJacColMajor<18>(nf, nc, diagonal, N, result);
            break;
        case 19:
            gradientJacColMajor<19>(nf, nc, diagonal, N, result);
            break;
        case 20:
            gradientJacColMajor<20>(nf, nc, diagonal, N, result);
            break;
        case 21:
            gradientJacColMajor<21>(nf, nc, diagonal, N, result);
            break;
        case 22:
            gradientJacColMajor<22>(nf, nc, diagonal, N, result);
            break;
        case 23:
            gradientJacColMajor<23>(nf, nc, diagonal, N, result);
            break;
        case 24:
            gradientJacColMajor<24>(nf, nc, diagonal, N, result);
            break;
        case 25:
            gradientJacColMajor<25>(nf, nc, diagonal, N, result);
            break;
        case 26:
            gradientJacColMajor<26>(nf, nc, diagonal, N, result);
            break;
        case 27:
            gradientJacColMajor<27>(nf, nc, diagonal, N, result);
            break;
        case 28:
            gradientJacColMajor<28>(nf, nc, diagonal, N, result);
            break;
        case 29:
            gradientJacColMajor<29>(nf, nc, diagonal, N, result);
            break;
        case 30:
            gradientJacColMajor<30>(nf, nc, diagonal, N, result);
            break;
        default:
            mexErrMsgTxt("%d derivatives not supported by backend.");
        }
    }
    return;
}


