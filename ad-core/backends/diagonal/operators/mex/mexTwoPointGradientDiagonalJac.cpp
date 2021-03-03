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
#else
    #include <mex.h>
#endif

// INPUTS:
//  - cell_diagonal<double> [nc x m] if column major or [m x nc] if row major)
//  - N<double>             [nf x 2]
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
    } else if (nin != 4) {
        status_code = -2;
        return "4 input arguments required: Diagonal, N, number of cells and rowMajor bool";
    } else if (nout > 1) {
        status_code = -3;
        return "Too many outputs requested. Function has a single output argument.";
    } else {
        // All ok.
        status_code = 0;
        return "";
    }
}

const char* dimensionCheck(const int nc, const int nrows, const int ncols, int & status_code){
    if(nrows != nc && ncols != nc){
        status_code = -5;
        return "Malformed input. No dimension of input matrix matches number of cells: Dimensions of diagonal matrix does not fit either RowMajor or ColMajor";
    } else {
        return "";
    }
}


/* Templated function for main operation */
template <bool rowMajor>
void gradientJac(const int nf, const int nc, const int m, const double * diagonal, const double * N, double * result){
    #pragma omp parallel for
    for (int i = 0; i < nf; i++) {
        int left = N[i] - 1;
        int right = N[i + nf] - 1;
        for (int j = 0; j < m; j++) {
            if (rowMajor) {
                result[i * 2 * m + j] = -diagonal[m * left + j];
                result[i * 2 * m + j + m] = diagonal[m * right + j];

            } else {
                result[j * nf + i]      = -diagonal[nc * j + left];
                result[j * nf + i + m*nf] =  diagonal[nc * j + right];
            }
        }
    }
    return;
}

template <int m, bool rowMajor>
void gradientJac(const int nf, const int nc, const double* diagonal, const double* N, double* result) {
    gradientJac<rowMajor>(nf, nc, m, diagonal, N, result);
}

template <bool rowMajor>
void gradientJacMain(const int m, const int nf, const int nc, const double * diagonal, const double * N, double * result){
    switch (m) {
        case 1:
            gradientJac<1, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 2:
            gradientJac<2, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 3:
            gradientJac<3, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 4:
            gradientJac<4, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 5:
            gradientJac<5, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 6:
            gradientJac<6, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 7:
            gradientJac<7, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 8:
            gradientJac<8, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 9:
            gradientJac<9, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 10:
            gradientJac<10, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 11:
            gradientJac<11, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 12:
            gradientJac<12, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 13:
            gradientJac<13, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 14:
            gradientJac<14, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 15:
            gradientJac<15, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 16:
            gradientJac<16, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 17:
            gradientJac<17, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 18:
            gradientJac<18, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 19:
            gradientJac<19, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 20:
            gradientJac<20, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 21:
            gradientJac<21, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 22:
            gradientJac<22, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 23:
            gradientJac<23, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 24:
            gradientJac<24, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 25:
            gradientJac<25, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 26:
            gradientJac<26, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 27:
            gradientJac<27, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 28:
            gradientJac<28, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 29:
            gradientJac<29, rowMajor>(nf, nc, diagonal, N, result);
            break;
        case 30:
            gradientJac<30, rowMajor>(nf, nc, diagonal, N, result);
            break;
        default:
            gradientJac<rowMajor>(nf, nc, m, diagonal, N, result);
    }
}



#ifdef MRST_OCTEXT
    /* OCT gateway */
    DEFUN_DLD (mexTwoPointGradientDiagonalJac, args, nargout,
               "Two point gradient average operator for MRST - diagonal Jacobian.")
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

        bool rowMajor = args(3).scalar_value();
        int nc = args(2).scalar_value();
        int nf = N_nd.rows();

        int nrows = diagonal_nd.rows();
        int ncols = diagonal_nd.cols();

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
            gradientJacMain<true>(m, nf, nc, diagonal, N, result);
        }else{
            gradientJacMain<false>(m, nf, nc, diagonal, N, result);
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

        double * diagonal = mxGetPr(prhs[0]);
        double * N = mxGetPr(prhs[1]);
        bool rowMajor = mxGetScalar(prhs[3]);

        int nc = mxGetScalar(prhs[2]);
        int nf = mxGetM(prhs[1]);

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
        // plhs[0] = mxCreateDoubleMatrix(outrow, outcol, mxREAL);
        plhs[0] = mxCreateUninitNumericMatrix(outrow, outcol, mxDOUBLE_CLASS, mxREAL);
        double* result = mxGetPr(plhs[0]);
        if (rowMajor){
            gradientJacMain<true>(m, nf, nc, diagonal, N, result);
        } else {
            gradientJacMain<false>(m, nf, nc, diagonal, N, result);
        }
        return;
    }
#endif


