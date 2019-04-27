#ifndef MRST_DUNEISTL_HPP
#define MRST_DUNEISTL_HPP

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

namespace mrst{
  template <int  bz>
  class BlockIlu0Solver{
  public:
    BlockIlu0Solver(){};
    void solve(double* result,std::string matrixfile,std::string rhsfile){
      std::cout << matrixfile << std::endl;
      std::cout << rhsfile << std::endl;
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
	//Dune::writeMatrixMarket(rhs,std::cout);
      }
      {
	std::ifstream infile(matrixfile);
	if(!infile){
	  throw std::runtime_error("Matrix file not read");
	}
	Dune::readMatrixMarket(matrix,infile);
	std::string tmpmatrix("matrixtmp.txt");
	std::ofstream outfile(tmpmatrix);
	if(!outfile){
	  throw std::runtime_error("Could not write");
	}
	Dune::writeMatrixMarket(matrix,outfile);
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
      //plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);
      //double* result = mxGetPr(plhs[0]);
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
    }
  };
}
#endif /* MRST_DUNEISTL_HPP */
