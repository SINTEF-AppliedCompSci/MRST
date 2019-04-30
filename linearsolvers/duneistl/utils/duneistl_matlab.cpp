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
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[] )

{
  
  
  //mxDouble *rhs;
  double *rhs;
  double *err;
  double *it_count;
  mwSize nnz;
  //mxDouble * id;
  //mxDouble * jd;
  //mxDouble * vald;
  double * id;
  double * jd;
  double * vald;
  // const mxArray * pa;

  if (nrhs != 8) {
    mexErrMsgTxt("7 input arguments required.");
  } else if (nlhs > 3) {
    mexErrMsgTxt("Wrong number of output arguments.");
  }

  //id = mxGetDoubles(prhs[0]);
  //jd = mxGetDoubles(prhs[1]);
  //vald = mxGetDoubles(prhs[2]);
  id = mxGetPr(prhs[0]);
  jd = mxGetPr(prhs[1]);
  vald = mxGetPr(prhs[2]);
  nnz = mxGetM(prhs[0]);
  //rhs = mxGetDoubles(prhs[3]);
  rhs = mxGetPr(prhs[3]);
  size_t rows = mxGetM(prhs[3]);
  //double mm = mxGetScalar(prhs[4]);
  //int m(mm);
  double bzm = mxGetScalar(prhs[4]);
  size_t bz(bzm);
  double tol = mxGetScalar(prhs[5]);
  double dmaxiter = mxGetScalar(prhs[6]);
  int maxiter(dmaxiter);
  //std::cout << "All values set" << std::endl;
  /* Assign a pointer to the output */
  //y = mxGetDoubles(plhs[0]);
  namespace pt = boost::property_tree;
  pt::ptree prm;
  {
    char *chopt = mxArrayToString(prhs[7]);  
    std::istringstream str(chopt);
    pt::read_json(str, prm);
    std::cout << "options in prm" << std::endl;
    std::ofstream file("options.json");
    pt::write_json(std::cout, prm);
    pt::write_json(file, prm);
  }
  std::vector<int> i(nnz);
  std::vector<int> j(nnz);
  std::vector<double> val(nnz);
  for(size_t kk=0; kk<val.size(); ++kk){
    i[kk] = id[kk];
    j[kk] = jd[kk];
    val[kk] = vald[kk];
  }
  std::vector<double> orhs(rows);
  for(size_t kk=0; kk<orhs.size(); ++kk){
    orhs[kk] = rhs[kk];
  }
  // std::cout << "Start solving " << std::endl;
  // std::cout << "Eq size " << rows << std::endl;
  // std::cout << "nnz org " << val.size() << std::endl;
  // std::cout << "tol " << tol << std::endl;
  // std::cout << "maxiter " << maxiter << std::endl;
  // std::cout << "******rhs*********" << std::endl;
  // for(auto x : orhs){
  //   std::cout << x << std::endl;
  // }
  // std::cout << "*******matrix**********" << std::endl;
  // for(int kk=0; kk < val.size(); ++kk){
  //   std::cout << i[kk] << " " << j[kk] << " " << val[kk] << std::endl;
  // }
  // std::cout << "**********************" << std::endl;

  plhs[0] = mxCreateDoubleMatrix(orhs.size(), 1, mxREAL);
  double* result = mxGetPr(plhs[0]);
  if(bz==1){    
    mrst::BlockIlu0Solver<1> solver;
    solver.solve(result, i, j, val, rows, orhs, tol, maxiter);
  }else if(bz == 2){
    mrst::BlockIlu0Solver<2> solver;
    solver.solve(result, i, j, val, rows, orhs, tol, maxiter); 
  }else if(bz ==3){
    mrst::BlockIlu0Solver<3> solver;
    solver.solve(result, i, j, val, rows, orhs, tol, maxiter); 
  }else{
    std::cout<< "BlockIlu0 solver not implemented for blocksize " << bz << std::endl;
  }
  return;
}
