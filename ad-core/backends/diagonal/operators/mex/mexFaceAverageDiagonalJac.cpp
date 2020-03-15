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


template <int m>
void faceAverageJac(const int nf, const int nc, const double * diagonal, const double * N, double * result){
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < nf; i++) {
        int left = N[i] - 1;
        int right = N[i + nf] - 1;
        for(int j=0;j<m;j++){
            result[j * nf + i] = 0.5 * diagonal[nc * j + left];
            result[j * nf + i + m*nf] = 0.5 * diagonal[nc * j + right];
        }
    }
}

/* MEX gateway */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
    // In: diagonal (nc x m), N (nf x 2)
    // Out: Face diagonal of (nf x 2*m)
    if (nrhs != 2) { 
	    mexErrMsgTxt("2 input arguments required."); 
    } else if (nlhs > 1) {
	    mexErrMsgTxt("Wrong number of output arguments."); 
    } 
    double * diagonal = mxGetPr(prhs[0]);
    double * N = mxGetPr(prhs[1]);

    int nc = mxGetM(prhs[0]);
    int m = mxGetN(prhs[0]);
    int nf = mxGetM(prhs[1]);
        
    plhs[0] = mxCreateDoubleMatrix(nf, 2*m, mxREAL);
    double * result = mxGetPr(plhs[0]);
    switch (m){
            case 1:
            faceAverageJac<1>(nf, nc, diagonal, N, result);
            break;
            case 2:
            faceAverageJac<2>(nf, nc, diagonal, N, result);
            break;
            case 3:
            faceAverageJac<3>(nf, nc, diagonal, N, result);
            break;
            case 4:
            faceAverageJac<4>(nf, nc, diagonal, N, result);
            break;
            case 5:
            faceAverageJac<5>(nf, nc, diagonal, N, result);
            break;
            case 6:
            faceAverageJac<6>(nf, nc, diagonal, N, result);
            break;
            case 7:
            faceAverageJac<7>(nf, nc, diagonal, N, result);
            break;
            case 8:
            faceAverageJac<8>(nf, nc, diagonal, N, result);
            break;
            case 9:
            faceAverageJac<9>(nf, nc, diagonal, N, result);
            break;
            case 10:
            faceAverageJac<10>(nf, nc, diagonal, N, result);
            break;
            case 11:
            faceAverageJac<11>(nf, nc, diagonal, N, result);
            break;
            case 12:
            faceAverageJac<12>(nf, nc, diagonal, N, result);
            break;
            case 13:
            faceAverageJac<13>(nf, nc, diagonal, N, result);
            break;
            case 14:
            faceAverageJac<14>(nf, nc, diagonal, N, result);
            break;
            case 15:
            faceAverageJac<15>(nf, nc, diagonal, N, result);
            break;
            case 16:
            faceAverageJac<16>(nf, nc, diagonal, N, result);
            break;
            case 17:
            faceAverageJac<17>(nf, nc, diagonal, N, result);
            break;
            case 18:
            faceAverageJac<18>(nf, nc, diagonal, N, result);
            break;
            case 19:
            faceAverageJac<19>(nf, nc, diagonal, N, result);
            break;
            case 20:
            faceAverageJac<20>(nf, nc, diagonal, N, result);
            break;
            case 21:
            faceAverageJac<21>(nf, nc, diagonal, N, result);
            break;
            case 22:
            faceAverageJac<22>(nf, nc, diagonal, N, result);
            break;
            case 23:
            faceAverageJac<23>(nf, nc, diagonal, N, result);
            break;
            case 24:
            faceAverageJac<24>(nf, nc, diagonal, N, result);
            break;
            case 25:
            faceAverageJac<25>(nf, nc, diagonal, N, result);
            break;
            case 26:
            faceAverageJac<26>(nf, nc, diagonal, N, result);
            break;
            case 27:
            faceAverageJac<27>(nf, nc, diagonal, N, result);
            break;
            case 28:
            faceAverageJac<28>(nf, nc, diagonal, N, result);
            break;
            case 29:
            faceAverageJac<29>(nf, nc, diagonal, N, result);
            break;
            case 30:
            faceAverageJac<30>(nf, nc, diagonal, N, result);
            break;
        default:
            mexErrMsgTxt("%d derivatives not supported by backend."); 
    }

    return;
}


