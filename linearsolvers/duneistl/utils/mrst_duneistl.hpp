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


namespace Dune
{


  
  namespace MatrixMarketImpl{
    template<typename T, typename A, int brows, int bcols, typename D>
    void makeSparseEntries(Dune::BCRSMatrix<Dune::FieldMatrix<T,brows,bcols>,A>& matrix,
			   std::vector<int>& i, std::vector<int>& j, std::vector<double>& val,
			   std::size_t entries,
			   const D&)
    {
      // addapted from readSparseEntries
      typedef Dune::BCRSMatrix<Dune::FieldMatrix<T,brows,bcols>,A> Matrix;
      std::vector<std::set<IndexData<D> > > rows(matrix.N()*brows);
      for(int kk=0; kk<i.size(); ++kk) {
	std::size_t row;
	IndexData<D> data;
	row= i[kk];
	assert(row/bcols<matrix.N());
	data.number = val[kk];
	data.index = j[kk];
	assert(data.index/bcols<matrix.M());
	rows[row].insert(data);
      }
      
      // Setup the matrix sparsity pattern
      int nnz=0;
      for(typename Matrix::CreateIterator iter=matrix.createbegin();
	  iter!= matrix.createend(); ++iter)
	{
	  for(std::size_t brow=iter.index()*brows, browend=iter.index()*brows+brows;
	      brow<browend; ++brow)
	    {
	      typedef typename std::set<IndexData<D> >::const_iterator Siter;
	      for(Siter siter=rows[brow].begin(), send=rows[brow].end();
		  siter != send; ++siter, ++nnz)
		iter.insert(siter->index/bcols);
	    }
	}
      
      //Set the matrix values
      matrix=0;
      
      MatrixValuesSetter<D,brows,bcols> Setter;
      
      Setter(rows, matrix);
    }
  } // end namespace MatrixMarketImpl

  template<typename T, typename A, int brows, int bcols>
  void makeMatrixMarket(Dune::BCRSMatrix<Dune::FieldMatrix<T,brows,bcols>,A>& matrix,
			std::vector<int> i, std::vector<int> j, std::vector<T> val,
			size_t rows,
			size_t cols,
			size_t entries)
  {
    // addapted from readMatrixMarket
    using namespace MatrixMarketImpl;
    //std::size_t rows, cols, entries;
    //std::size_t nnz, blockrows, blockcols;
    //std::tie(blockrows, blockcols, nnz) = calculateNNZ<brows, bcols>(rows, cols, entries, header);
    std::size_t blockrows=rows/brows;
    std::size_t blockcols=cols/bcols;
    std::size_t blocksize=brows*bcols;
    std::size_t blockentries=0;
    blockentries = entries/blocksize;//nnz
    matrix.setSize(blockrows, blockcols);
    matrix.setBuildMode(Dune::BCRSMatrix<Dune::FieldMatrix<T,brows,bcols>,A>::row_wise);
    makeSparseEntries(matrix, i, j, val, entries, NumericWrapper<T>());
  }
}



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
    
    void solve(double* result,
	       std::vector<int>& i,
	       std::vector<int>& j,
	       std::vector<double>& val,
	       size_t rows,
	       std::vector<double>& orhs,
	       double tol,
	       int maxiter){
      typedef Dune::BCRSMatrix< Dune::FieldMatrix< double, bz, bz > > MatrixType;
      typedef Dune::BlockVector< Dune::FieldVector< double, bz > > VectorType;
      
      // copy rhs to block vector
      VectorType rhs(rows/bz);
      {
	int lind=0;
	for(size_t ic = 0; ic < rhs.size(); ic++){
	  for(size_t ib = 0; ib < bz; ib++){
	    rhs[ic][ib] = orhs[lind];
	    lind++;
	  }
	}
      }
      
      // make block matrix
      MatrixType matrix;
      makeMatrixMarket(matrix,
		       i,j,val,
		       rows,
		       rows,
		       val.size());
  
      Dune::Timer perfTimer;
      perfTimer.start();
      //double tol = 1e-4;
      //int maxiter = 200;
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
      {
	int lind = 0;
	for(size_t ic = 0; ic < rhs.size(); ic++){
	  for(size_t ib = 0; ib < bz; ib++){
	    result[lind] = x[ic][ib];
	    lind++;
	  }
	}
      }
    }
  };
}
#endif /* MRST_DUNEISTL_HPP */
