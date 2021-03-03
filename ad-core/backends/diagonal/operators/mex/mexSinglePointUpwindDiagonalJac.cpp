//
// include necessary system headers
//
#include <cmath>
#include <array>
#ifdef _OPENMP
    #include <omp.h>
#endif
#include <iostream>
#include <chrono>
#ifdef MRST_OCTEXT
    #include <octave/oct.h>
    #include <octave/dMatrix.h>
#else
    #include <mex.h>
#endif

#define TIME_NOW std::chrono::high_resolution_clock::now

// INPUTS:
//  - cell_diagonal<double> [nc x m] if column major or [m x nc] if row major)
//  - N<double>             [nf x 2]
//  - flag<bool>            [nf x 1]
//  - nc<double>            [scalar]
//  - rowMajor<bool>        [scalar]
// OUTPUT:
//  - face_diagonal<double> [nf x 2*m] if column major or [2*m x nf] if row major
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
    } else if (nin != 5) {
        status_code = -2;
        return "5 input arguments required: Diagonal, N, flag, number of cells and rowMajor bool";
    } else if (nout > 1) {
        status_code = -3;
        return "Too many outputs requested. Function has a single output argument.";
    } else {
        // All ok.
        status_code = 0;
        return "";
    }
}


void copyElements(const int offset_face, const int offset_cell, const int m, const double * celldata, double * facedata) {
    for (int der = 0; der < m; der++) {
        facedata[offset_face + der] = celldata[offset_cell + der];
    }
}

void zeroElements(const int offset_face, const int m, double * facedata) {
    for (int der = 0; der < m; der++) {
        facedata[offset_face + der] = 0;
    }
}

/* Templated function for main operation */
template <bool rowMajor, class logic_type>
void upwindJac(const int nf, const int nc, const int m, const logic_type* flag, const double* diagonal, const double* N, double* result) {
    if (rowMajor){
        #pragma omp parallel for schedule(static)
        for (int face = 0; face < nf; face++) {
            // Copy / zero out logic. We are working with uninitialized arrays so we need to set both zero and value.
            int copy_offset, zero_offset, copy_cell;
            if(flag[face]){
                copy_offset = 0;
                zero_offset = m;
                copy_cell = N[face]-1;
            }else{
                copy_offset = m;
                zero_offset = 0;
                copy_cell = N[face + nf]-1;
            }
            copyElements(face*2*m + copy_offset, m * copy_cell, m, diagonal, result);
            zeroElements(face*2*m + zero_offset, m, result);
        }
    }else{
        #pragma omp parallel for schedule(static)
        for(int i=0;i<nf;i++){
            int left = N[i] - 1;
            int right = N[i + nf] - 1;

            if (flag[i]) {
                for (int j = 0; j < m; j++) {
                    result[j * nf  + i] = diagonal[nc * j + left];
                    result[j * nf + m*nf + i] = 0;
                }
            }
            else {
                for (int j = 0; j < m; j++) {
                    result[j * nf  + m*nf + i] = diagonal[nc * j + right];
                    result[j * nf + i] = 0;
                }
            }
        }
        return;
    }
}

template <int m, bool rowMajor, class logic_type>
void upwindJac(const int nf, const int nc, const logic_type * flag, const double* diagonal, const double* N, double* result) {
    upwindJac<rowMajor, logic_type>(nf, nc, m, flag, diagonal, N, result);
}

template <bool rowMajor, class logic_type>
void upwindJacMain(const int m, const int nf, const int nc, const logic_type * flag, const double * diagonal, const double * N, double * result){
    switch (m) {
        case 1:
            upwindJac<1, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 2:
            upwindJac<2, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 3:
            upwindJac<3, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 4:
            upwindJac<4, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 5:
            upwindJac<5, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 6:
            upwindJac<6, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 7:
            upwindJac<7, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 8:
            upwindJac<8, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 9:
            upwindJac<9, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 10:
            upwindJac<10, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 11:
            upwindJac<11, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 12:
            upwindJac<12, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 13:
            upwindJac<13, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 14:
            upwindJac<14, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 15:
            upwindJac<15, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 16:
            upwindJac<16, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 17:
            upwindJac<17, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 18:
            upwindJac<18, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 19:
            upwindJac<19, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 20:
            upwindJac<20, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 21:
            upwindJac<21, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 22:
            upwindJac<22, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 23:
            upwindJac<23, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 24:
            upwindJac<24, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 25:
            upwindJac<25, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 26:
            upwindJac<26, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 27:
            upwindJac<27, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 28:
            upwindJac<28, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 29:
            upwindJac<29, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        case 30:
            upwindJac<30, rowMajor>(nf, nc, flag, diagonal, N, result);
            break;
        default:
            upwindJac<rowMajor>(nf, nc, m, flag, diagonal, N, result);
    }
}

#ifdef MRST_OCTEXT
    /* OCT gateway */
    DEFUN_DLD (mexSinglePointUpwindDiagonalJac, args, nargout,
               "Single point upwind operator for MRST - diagonal Jacobian.")
    {
        //auto start = high_resolution_clock::now(); 

        const int nrhs = args.length();
        const int nlhs = nargout;

        int status_code = 0;
        auto msg = inputCheck(nrhs, nlhs, status_code);
        
        if(status_code < 0){
            // Some kind of error
            error(msg);
        } else if (status_code == 1){
            // Early return
            return octave_value_list();
        }
        
        const NDArray diagonal_nd = args(0).array_value();
        const NDArray N_nd = args(1).array_value();

        const double * diagonal = diagonal_nd.data();
        const double * N = N_nd.data();

        bool rowMajor = args(4).scalar_value();
        int nc = args(3).scalar_value();
        int nf = N_nd.rows();

        int nrows = diagonal_nd.rows();
        int ncols = diagonal_nd.cols();

        const NDArray flag_nd = args(2).array_value();
        const double * flag = flag_nd.data();

        int outrow, outcol, m;
        if(rowMajor){
            m = nrows;
            outrow = 2*m;
            outcol = nf;
        }else{
            m = ncols;
            outrow = nf;
            outcol = 2*m;
        }
        NDArray output({outrow, outcol});
        double * result = output.fortran_vec();

        if (rowMajor){
            upwindJacMain<true, double>(m, nf, nc, flag, diagonal, N, result);
        } else {
            upwindJacMain<false, double>(m, nf, nc, flag, diagonal, N, result);
        }
        return octave_value (output);
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
        double* diagonal = mxGetPr(prhs[0]);
        double* N = mxGetPr(prhs[1]);
        mxLogical* flag = mxGetLogicals(prhs[2]);
        int nc = mxGetScalar(prhs[3]);
        int nf = mxGetM(prhs[1]);

        bool rowMajor = mxGetScalar(prhs[4]);

        // Dimensions of diagonals - figure out if we want row or column major solver
        int nrows = mxGetM(prhs[0]);
        int ncols = mxGetN(prhs[0]);

        int outrow, outcol, m;
        if(rowMajor){
            m = nrows;
            outrow = 2*m;
            outcol = nf;
        } else {
            m = ncols;
            outrow = nf;
            outcol = 2*m;
        }
        plhs[0] = mxCreateUninitNumericMatrix(outrow, outcol, mxDOUBLE_CLASS, mxREAL);
        // plhs[0] = mxCreateDoubleMatrix(outrow, outcol, mxREAL);
        double* result = mxGetPr(plhs[0]);
        if (rowMajor){
            upwindJacMain<true, mxLogical>(m, nf, nc, flag, diagonal, N, result);
        } else {
            upwindJacMain<false, mxLogical>(m, nf, nc, flag, diagonal, N, result);
        }
        return;
    }
#endif


