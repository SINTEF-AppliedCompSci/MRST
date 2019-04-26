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
//#include <boost/program_options.hpp>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixmarket.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/umfpack.hh>
#include <dune/istl/solvers.hh>
/* MEX gateway */

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
    constexpr int bz = 3;
    char *matrixfilename = mxArrayToString(prhs[0]);
    std::string matrixfile(matrixfilename);
    char *rhsfilename = mxArrayToString(prhs[1]);
    std::string rhsfile(rhsfilename);
    std::cout << matrixfile << std::endl;
    typedef Dune::BCRSMatrix< Dune::FieldMatrix< double, bz, bz > > MatrixType;
    typedef Dune::BlockVector< Dune::FieldVector< double, bz > > VectorType;
    MatrixType matrix;
    VectorType rhs;
    {
      std::ifstream infile(rhsfile);
      if(!infile){
	throw std::runtime_error("Rhs file not read");
      }
      Dune::readMatrixMarket(rhs,infile);
    }
    {
      std::ifstream infile(matrixfile);
      if(!infile){
	throw std::runtime_error("Matrix file not read");
      }
      Dune::readMatrixMarket(matrix,infile);
    }
    Dune::Timer perfTimer;
    perfTimer.start();
    double tol = 1e-4;
    int maxiter = 200;
    int verbosity = 10;
    Dune::SeqILU0<MatrixType, VectorType, VectorType> preconditioner(matrix, 1.0);
    Dune::MatrixAdapter<MatrixType, VectorType, VectorType> linearOperator(matrix);
    Dune::BiCGSTABSolver<VectorType> linsolver(linearOperator,
					       preconditioner,
					       tol, // desired residual reduction factor
					       maxiter, // maximum number of iterations
					       verbosity); 
    

    int m = bz*rhs.size();
    plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);
    double* result = mxGetPr(plhs[0]);
    VectorType x(rhs.size());
    Dune::InverseOperatorResult res;
    linsolver.apply(x, rhs, res);
    double time = perfTimer.stop();
//x[0]=5;
    int i = 0;
    for(size_t ic = 0; ic < rhs.size(); ic++){
        for(size_t ib = 0; ib < bz; ib++){
         result[i] = x[ic][ib];
          i++;
         }
    }
    // int    iters;
    // double error;
    // std::vector<double> x(M, 0.0);
    
   
    
    // for(int ix=0; ix < M; ix++){
    //     result[ix] = x[ix];
    // }
    // x.clear();
    // b.clear();
    // err[0] = error;
    // it_count[0] = iters;
    return;
}
