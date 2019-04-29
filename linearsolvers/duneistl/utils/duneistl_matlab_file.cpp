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
