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
void upwindJac(const int nf, const int nc, const mxLogical * flag, const double * diagonal, const double * N, double * result){
    #ifdef _MSC_VER
        #pragma omp parallel for
    #else
        #pragma omp parallel for collapse(2)
    #endif
    for(int j=0;j<m;j++){
        for(int i=0;i<2*nf;i++){
            int cell_inx = N[i]-1;
            if(flag[i % nf] == i < nf){
                result[j*2*nf + i] = diagonal[nc*j + cell_inx];
            }
        }
    }
    return;
}

/* MEX gateway */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
    // In: diagonal (nc x m), N (nf x 2), flag (nf x 1) bool
    // Out: Face diagonal of (2*nf x m)
    if (nrhs != 3) { 
	    mexErrMsgTxt("3 input arguments required."); 
    } else if (nlhs > 1) {
	    mexErrMsgTxt("Wrong number of output arguments."); 
    } 
    double * diagonal = mxGetPr(prhs[0]);
    double * N = mxGetPr(prhs[1]);
    mxLogical  * flag = mxGetLogicals(prhs[2]);

    int nc = mxGetM(prhs[0]);
    int m = mxGetN(prhs[0]);
    int nf = mxGetM(prhs[1]);
    
    // printf("%d cells with %d faces and %d derivatives \n", nc, nf, m);
    
    plhs[0] = mxCreateDoubleMatrix(2*nf, m, mxREAL);
    double * result = mxGetPr(plhs[0]);
    switch (m){
            case 1:
            upwindJac<1>(nf, nc, flag, diagonal, N, result);
            break;
            case 2:
            upwindJac<2>(nf, nc, flag, diagonal, N, result);
            break;
            case 3:
            upwindJac<3>(nf, nc, flag, diagonal, N, result);
            break;
            case 4:
            upwindJac<4>(nf, nc, flag, diagonal, N, result);
            break;
            case 5:
            upwindJac<5>(nf, nc, flag, diagonal, N, result);
            break;
            case 6:
            upwindJac<6>(nf, nc, flag, diagonal, N, result);
            break;
            case 7:
            upwindJac<7>(nf, nc, flag, diagonal, N, result);
            break;
            case 8:
            upwindJac<8>(nf, nc, flag, diagonal, N, result);
            break;
            case 9:
            upwindJac<9>(nf, nc, flag, diagonal, N, result);
            break;
            case 10:
            upwindJac<10>(nf, nc, flag, diagonal, N, result);
            break;
            case 11:
            upwindJac<11>(nf, nc, flag, diagonal, N, result);
            break;
            case 12:
            upwindJac<12>(nf, nc, flag, diagonal, N, result);
            break;
            case 13:
            upwindJac<13>(nf, nc, flag, diagonal, N, result);
            break;
            case 14:
            upwindJac<14>(nf, nc, flag, diagonal, N, result);
            break;
            case 15:
            upwindJac<15>(nf, nc, flag, diagonal, N, result);
            break;
            case 16:
            upwindJac<16>(nf, nc, flag, diagonal, N, result);
            break;
            case 17:
            upwindJac<17>(nf, nc, flag, diagonal, N, result);
            break;
            case 18:
            upwindJac<18>(nf, nc, flag, diagonal, N, result);
            break;
            case 19:
            upwindJac<19>(nf, nc, flag, diagonal, N, result);
            break;
            case 20:
            upwindJac<20>(nf, nc, flag, diagonal, N, result);
            break;
            case 21:
            upwindJac<21>(nf, nc, flag, diagonal, N, result);
            break;
            case 22:
            upwindJac<22>(nf, nc, flag, diagonal, N, result);
            break;
            case 23:
            upwindJac<23>(nf, nc, flag, diagonal, N, result);
            break;
            case 24:
            upwindJac<24>(nf, nc, flag, diagonal, N, result);
            break;
            case 25:
            upwindJac<25>(nf, nc, flag, diagonal, N, result);
            break;
            case 26:
            upwindJac<26>(nf, nc, flag, diagonal, N, result);
            break;
            case 27:
            upwindJac<27>(nf, nc, flag, diagonal, N, result);
            break;
            case 28:
            upwindJac<28>(nf, nc, flag, diagonal, N, result);
            break;
            case 29:
            upwindJac<29>(nf, nc, flag, diagonal, N, result);
            break;
            case 30:
            upwindJac<30>(nf, nc, flag, diagonal, N, result);
            break;
        default:
            mexErrMsgTxt("%d derivatives not supported by backend."); 
    }
    return;
}


