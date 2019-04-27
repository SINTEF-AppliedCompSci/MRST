//
// include necessary system headers
//
#define _USE_MATH_DEFINES
#include <cmath>
#include <mex.h>
#include <array>
#include "matrix.h"


#include <iostream>
#include <fstream>
/* MEX interfaces */
/* Block system support */
#include <iomanip>
#include <cmath>
#include <limits>
#include "mrst_duneistl.hpp"
//#include <boost/program_options.hpp>
// #include <dune/common/parallel/mpihelper.hh>
// #include <dune/istl/bcrsmatrix.hh>
// #include <dune/istl/matrixmarket.hh>
// #include <dune/common/fmatrix.hh>
// #include <dune/istl/solvers.hh>
// #include <dune/istl/preconditioners.hh>
// #include <dune/istl/umfpack.hh>
// #include <dune/istl/solvers.hh>


void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[] )

{
  
    // double *result;
    // double *rhs;
    // double *err;
    // double *it_count;
    // mwSize m,n,nnz;
    // mwIndex * cols;
    // mwIndex * rows;
    // const mxArray * pa;

    // double * entries;
    // std::string relaxParam;
    // std::string coarsenParam;

    // if (nrhs != 6) {
    // 	    mexErrMsgTxt("6 input arguments required.");
    // } else if (nlhs > 3) {
    // 	    mexErrMsgTxt("Wrong number of output arguments.");
    // }

    // m = mxGetM(prhs[0]);
    // n = mxGetN(prhs[0]);
    // int M = (int)m;

    // if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||  !mxIsSparse(prhs[0]) ) {
    // 	    mexErrMsgTxt("Matrix should be a real sparse matrix.");
    //     return;
    // }
    // if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) {
    // 	    mexErrMsgTxt("Right hand side must be real double column vector.");
    //     return;
    // }
    // // main();
    // plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    // plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);

    //input_buf = mxArrayToString(prhs[0]);
    // result = mxGetPr(plhs[0]);
    // err = mxGetPr(plhs[1]);
    // it_count = mxGetPr(plhs[2]);

    // cols    = mxGetJc(prhs[0]);
    // rows    = mxGetIr(prhs[0]);
    // entries = mxGetPr(prhs[0]);
    // nnz  = mxGetNzmax(prhs[0]);
    // rhs     = mxGetPr(prhs[1]);
    // pa = prhs[2];
    // double tolerance = mxGetScalar(prhs[3]);
    // int maxiter = mxGetScalar(prhs[4]);
    // int solver_strategy_id = mxGetScalar(prhs[5]);
    
    char *matrixfilename = mxArrayToString(prhs[0]);
    char *rhsfilename = mxArrayToString(prhs[1]);
    //int m  = mxArrayToInt(prhs[2]);
    double mm = mxGetScalar(prhs[2]);
    int m(mm);
    double bzm = mxGetScalar(prhs[3]);
    int bz(bzm);
    /* Assign a pointer to the output */
    //y = mxGetDoubles(plhs[0]);
    std::string rhsfile(rhsfilename);
    std::string matrixfile(matrixfilename);
    plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);
    double* result = mxGetPr(plhs[0]);

    if(bz==1){    
      mrst::BlockIlu0Solver<1> solver;
      solver.solve(result,matrixfile,rhsfile);
    }else if(bz == 2){
      mrst::BlockIlu0Solver<2> solver;
      solver.solve(result,matrixfile,rhsfile); 
    }else if(bz ==3){
      mrst::BlockIlu0Solver<3> solver;
      solver.solve(result,matrixfile,rhsfile); 
    }else{
      std::cout<< "BlockIlu0 solver not implemented for blocksize " << bz << std::endl;
    }
    return;
}
