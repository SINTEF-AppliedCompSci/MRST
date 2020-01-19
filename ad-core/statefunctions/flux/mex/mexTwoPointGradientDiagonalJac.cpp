//
// include necessary system headers
//
#include <cmath>
#include <mex.h>
#include <array>

#include <iostream>
/* Templated function for main operation */

template <int m>
void gradientJac(const int nf, const int nc, const double * diagonal, const double * N, double * result){
    #pragma omp parallel for
    // Loop over each face
    for(int i=0;i<nf;i++){
        int cell_left = N[i]-1;
        int cell_right = N[i + nf]-1;
        // Loop over derivatives
        for(int j=0;j<m;j++){
            result[i*m + j] = -diagonal[m*cell_left + j];
            result[(i + nf)*m + j] = diagonal[m*cell_right + j];
        }
    }
    return;
}
/* MEX gateway */
void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
    // In: diagonal (m x nc), N (nf x 2)
    // Out: Face diagonal of (m x 2*nf)
    if (nrhs != 2) { 
	    mexErrMsgTxt("2 input arguments required."); 
    } else if (nlhs > 1) {
	    mexErrMsgTxt("Wrong number of output arguments."); 
    } 
    double * diagonal = mxGetPr(prhs[0]);
    double * N = mxGetPr(prhs[1]);

    int nc = mxGetN(prhs[0]);
    int m = mxGetM(prhs[0]);
    int nf = mxGetM(prhs[1]);
    
    // mexPrintf("%d cells with %d faces and %d derivatives \n", nc, nf, m);
    
    plhs[0] = mxCreateDoubleMatrix(m, 2*nf, mxREAL);
    double * result = mxGetPr(plhs[0]);
    switch (m){
            case 1:
            gradientJac<1>(nf, nc, diagonal, N, result);
            break;
            case 2:
            gradientJac<2>(nf, nc, diagonal, N, result);
            break;
            case 3:
            gradientJac<3>(nf, nc, diagonal, N, result);
            break;
            case 4:
            gradientJac<4>(nf, nc, diagonal, N, result);
            break;
            case 5:
            gradientJac<5>(nf, nc, diagonal, N, result);
            break;
            case 6:
            gradientJac<6>(nf, nc, diagonal, N, result);
            break;
            case 7:
            gradientJac<7>(nf, nc, diagonal, N, result);
            break;
            case 8:
            gradientJac<8>(nf, nc, diagonal, N, result);
            break;
            case 9:
            gradientJac<9>(nf, nc, diagonal, N, result);
            break;
            case 10:
            gradientJac<10>(nf, nc, diagonal, N, result);
            break;
            case 11:
            gradientJac<11>(nf, nc, diagonal, N, result);
            break;
            case 12:
            gradientJac<12>(nf, nc, diagonal, N, result);
            break;
            case 13:
            gradientJac<13>(nf, nc, diagonal, N, result);
            break;
            case 14:
            gradientJac<14>(nf, nc, diagonal, N, result);
            break;
            case 15:
            gradientJac<15>(nf, nc, diagonal, N, result);
            break;
            case 16:
            gradientJac<16>(nf, nc, diagonal, N, result);
            break;
            case 17:
            gradientJac<17>(nf, nc, diagonal, N, result);
            break;
            case 18:
            gradientJac<18>(nf, nc, diagonal, N, result);
            break;
            case 19:
            gradientJac<19>(nf, nc, diagonal, N, result);
            break;
            case 20:
            gradientJac<20>(nf, nc, diagonal, N, result);
            break;
            case 21:
            gradientJac<21>(nf, nc, diagonal, N, result);
            break;
            case 22:
            gradientJac<22>(nf, nc, diagonal, N, result);
            break;
            case 23:
            gradientJac<23>(nf, nc, diagonal, N, result);
            break;
            case 24:
            gradientJac<24>(nf, nc, diagonal, N, result);
            break;
            case 25:
            gradientJac<25>(nf, nc, diagonal, N, result);
            break;
            case 26:
            gradientJac<26>(nf, nc, diagonal, N, result);
            break;
            case 27:
            gradientJac<27>(nf, nc, diagonal, N, result);
            break;
            case 28:
            gradientJac<28>(nf, nc, diagonal, N, result);
            break;
            case 29:
            gradientJac<29>(nf, nc, diagonal, N, result);
            break;
            case 30:
            gradientJac<30>(nf, nc, diagonal, N, result);
            break;
        default:
            mexErrMsgTxt("%d derivatives not supported by backend."); 
    }
    
    return;
}


