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

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp> 

#include <dune/istl/paamg/fastamg.hh>
#include <dune/istl/paamg/amg.hh>
//#include "PressureSolverPolicy.hpp"
//#include "PressureTransferPolicy.hpp"

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
  template<class MatrixType,class VectorType>
  std::shared_ptr< Dune::Preconditioner<VectorType,VectorType> > getSeqPreconditioner(Dune::MatrixAdapter<MatrixType, VectorType, VectorType> & linearoperator,
										      boost::property_tree::ptree& prm)
  {
    auto& matrix = linearoperator.getmat();
    std::shared_ptr< Dune::Preconditioner<VectorType,VectorType> > preconditioner;
    int verbosity = 10;
    double w=prm.get<int>("w");;
    int n=prm.get<int>("n");
    std::string precond(prm.get<std::string>("preconditioner"));
    if(precond == "ILU0"){
      preconditioner.reset(new Dune::SeqILU0<MatrixType, VectorType, VectorType>(matrix,w));
    }else if(precond == "Jac"){
      preconditioner.reset(new Dune::SeqJac<MatrixType, VectorType, VectorType>(matrix,n, w));
    }else if(precond == "GS"){
      preconditioner.reset(new Dune::SeqGS<MatrixType, VectorType, VectorType>(matrix,n, w));
    }else if(precond == "SOR"){
      preconditioner.reset(new Dune::SeqSOR<MatrixType, VectorType, VectorType>(matrix,n, w));
    }else if(precond == "SSOR"){
      preconditioner.reset(new Dune::SeqSSOR<MatrixType, VectorType, VectorType>(matrix,n, w));
    }else if(precond == "ILUn"){
      preconditioner.reset(new Dune::SeqILUn<MatrixType, VectorType, VectorType>(matrix,n, w));
    }else{
      std::string msg("No such seq preconditioner : ");
      msg += precond;
	throw std::runtime_error(msg);;
    }
    return preconditioner;
  }

  template<class Smoother,class MatrixType,class VectorType>
  std::shared_ptr< Dune::Preconditioner<VectorType,VectorType> >
  makeAmgPreconditioner(Dune::MatrixAdapter<MatrixType, VectorType, VectorType> & linearoperator,
			boost::property_tree::ptree& global_prm){
      boost::property_tree::ptree prm = global_prm.get_child("amg");
      typedef Dune::MatrixAdapter<MatrixType, VectorType, VectorType> OperatorType;
      std::shared_ptr< Dune::Preconditioner<VectorType,VectorType> > preconditioner;
      typedef Dune::Amg::AggregationCriterion<
	Dune::Amg::SymmetricMatrixDependency<MatrixType,Dune::Amg::FirstDiagonal> > CriterionBase;
      typedef Dune::Amg::CoarsenCriterion<CriterionBase> Criterion;
      int coarsenTarget = prm.get<int>("coarsenTarget");
      int ml = prm.get<int>("maxlevel");
      Criterion criterion(15,coarsenTarget);
      criterion.setDefaultValuesIsotropic(2);
      criterion.setAlpha(.67);
      criterion.setBeta(1.0e-4);
      criterion.setMaxLevel(ml);
      criterion.setSkipIsolated(false);
      Dune::Amg::Parameters parms;
      if(global_prm.get<std::string>("preconditioner") == "famg"){ 
	preconditioner.reset(new Dune::Amg::FastAMG<OperatorType,VectorType>(linearoperator, criterion, parms));
      }else{
	typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments SmootherArgs;
	SmootherArgs smootherArgs;
	smootherArgs.iterations = 1;
	//smootherArgs.overlap=SmootherArgs::vertex;
	//smootherArgs.overlap=SmootherArgs::none;
	//smootherArgs.overlap=SmootherArgs::aggregate;
	smootherArgs.relaxationFactor = 1;
	//Smoother smoother(linearoperator.getmat(), 1);// 1);
	preconditioner.reset(
			     new
			     Dune::Amg::AMG<OperatorType,VectorType,Smoother>(linearoperator,
									      criterion,
									      smootherArgs));
	
      }
      return preconditioner;
  }
  



  template<class MatrixType,class VectorType>
  std::shared_ptr< Dune::Preconditioner<VectorType,VectorType> >
  makeAmgPreconditioners(Dune::MatrixAdapter<MatrixType, VectorType, VectorType> & linearoperator,
			 boost::property_tree::ptree& prm){
    std::shared_ptr< Dune::Preconditioner<VectorType,VectorType> > preconditioner;
    int verbosity = 10;
    std::string precond = prm.get<std::string>("amg.smoother");
    if(precond == "ILU0"){
      preconditioner = makeAmgPreconditioner<
	Dune::SeqILU0<MatrixType, VectorType, VectorType>,
	MatrixType, VectorType>(linearoperator,prm);
    }else if(precond == "Jac"){
      preconditioner = makeAmgPreconditioner<
	Dune::SeqJac<MatrixType, VectorType, VectorType>,
	MatrixType, VectorType>(linearoperator,prm);
    // }else if(precond == "GS"){
    //   preconditioner = makeAmgPreconditioner<
    // 	Dune::SeqGS<MatrixType, VectorType, VectorType>,
    // 	MatrixType, VectorType>(linearoperator,prm);
    }else if(precond == "SOR"){
      preconditioner = makeAmgPreconditioner<
	Dune::SeqSOR<MatrixType, VectorType, VectorType>,
	MatrixType, VectorType>(linearoperator,prm);
    }else if(precond == "SSOR"){
      preconditioner = makeAmgPreconditioner<
	Dune::SeqSSOR<MatrixType, VectorType, VectorType>,
	MatrixType, VectorType>(linearoperator,prm);
    }else if(precond == "ILUn"){
      preconditioner = makeAmgPreconditioner<
	Dune::SeqILUn<MatrixType, VectorType, VectorType>,
	MatrixType, VectorType>(linearoperator,prm);
    }else{
      std::string msg("No such seq preconditioner : ");
      msg += precond;
      throw std::runtime_error(msg);;
    }
    return preconditioner;
  }
  


    template <class int bz>
    class OwningTwoLevelPreconditioner : public Dune::Preconditioner<Dune::BlockVector< Dune::FieldVector< double, bz > >,
                                                                     Dune::BlockVector< Dune::FieldVector< double, bz > >>
    {
    public:
        typedef boost::property_tree::ptree pt;
        typedef Dune::BCRSMatrix< Dune::FieldMatrix< double, bz, bz > > MatrixType;
        typedef Dune::BlockVector< Dune::FieldVector< double, bz > > VectorType;
        typedef Dune::MatrixAdapter<MatrixType, VectorType, VectorType> OperatorType;

        OwningTwoLevelPreconditioner(const OperatorType& linearoperator, const pt& prm)
            : finesmoother_(makePreconditioner(prm.get_child("finesmoother")))
            , weights_
            , levelTransferPolicy_(Communication(), wei ghts_)
            , coarseSolverPolicy_(prm.get_child("coarsesolver"))
            , twolevel_method_(linearoperator,
                               *finesmoother_,
                               levelTransferPolicy_,
                               coarseSolverPolicy_,
                               0,
                               1);
        {
        }

        virtual void pre(VectorType& x, VectorType& b)
        {
            twolevel_method_.pre(x, b);
        }
        virtual void apply(VectorType& v, const VectorType& d)
        {
            twolevel_method_.apply(v, d);
        }
        virtual void post(VectorType& x)
        {
            twolevel_method_.post(x);
        }
    private:
        // for cpr
        typedef Dune::BCRSMatrix< Dune::FieldMatrix< double, 1, 1 > > PressureMatrixType;
        typedef Dune::BlockVector< Dune::FieldVector< double, 1 > > PressureVectorType;
        typedef Dune::MatrixAdapter<PressureMatrixType, PressureVectorType, PressureVectorType> CoarseOperatorType;
        using Communication =  Dune::Amg::SequentialInformation;
        constexpr int pressureVarIndex = 0;
        using LevelTransferPolicy = Opm::PressureTransferPolicy<OperatorType,
                                                                CoarseOperatorType,
                                                                Communication,
                                                                pressureVarIndex>;
        using CoarseSolverPolicy   =
          Dune::Amg::PressureSolverPolicy<CoarseOperatorType,
                                          LevelTransferPolicy>;
        using TwoLevelMethod =
            Dune::Amg::TwoLevelMethod<FineOperatorType,
                                      CoarseSolverPolicy,
                                      FineSmootherType>;
        std::shared_ptr< Dune::Preconditioner<VectorType,VectorType> > finesmoother_;
        VectorType weights_;
        std::unique_ptr<LevelTransferPolicy> levelTransferPolicy_;
        std::unique_ptr<CoarseSolverPolicy> coarseSolverPolicy_;
        TwoLevelMethod twolevel_method_;
    };



    template<class MatrixType,class VectorType>
    std::shared_ptr< Dune::Preconditioner<VectorType,VectorType> >
    makeTwoLevelPreconditioner(Dune::MatrixAdapter<MatrixType, VectorType, VectorType>& linearoperator,
                               boost::property_tree::ptree& prm)
    {
         int pressureVarIndex = 0;
         boost::property_tree::ptree prm_cpr = prm.get_child("cpr");
         boost::property_tree::ptree prm_finesmoother = prm_cpr.get_child("finesmoother");
         finesmoother_ = getSeqPreconditioner<MatrixType,VectorType>(matrix_,prm_finesmoother);
         Opm::Amg::getQuasiImpesWeights(matrix_, pressureVarIndex, weights_);
         levelTransferPolicy_.reset(new LevelTransferPolicy(comm_, weights_));
         boost::property_tree::ptree prm_coarsesolver = prm_cpr.get_child("coarsesolver");
         coarseSolverPolicy_.reset(new CoarseSolverPolicy(prm_coarsesolver));
         std::shared_ptr< Dune::Preconditioner<VectorType,VectorType> > preconditioner;
         preconditioner.reset(new TwoLevelMethod(linearoperator_,
                                                 finesmoother_,
                                                 *levelTransferPolicy_,
                                                 *coarseSolverPolicy_,
                                                 0,
                                                 1));
         return preconditioner;
    }

    template<class MatrixType,class VectorType>
    std::shared_ptr< Dune::Preconditioner<VectorType,VectorType> >
    makePreconditioner(Dune::MatrixAdapter<MatrixType, VectorType, VectorType>& linearoperator,
                       boost::property_tree::ptree& prm)
    {
        if ((prm.get<std::string>("preconditioner") == "famg") or
            (prm.get<std::string>("preconditioner") == "amg")) {
            return makeAmgPreconditioners<MatrixType, VectorType>(linearoperator, prm);
        } else if( prm.get<std::string>("preconditioner") == "cpr") {
            return makeTwoLevelPreconditioner(linearoperator, prm);
        } else {
            return getSeqPreconditioner<MatrixType,VectorType>(linearoperator, prm);
        }
    }



  template <int  bz>
  class BlockIlu0Solver{
  public:
    typedef boost::property_tree::ptree pt;
    typedef Dune::BCRSMatrix< Dune::FieldMatrix< double, bz, bz > > MatrixType;
    typedef Dune::BlockVector< Dune::FieldVector< double, bz > > VectorType;
    typedef Dune::MatrixAdapter<MatrixType, VectorType, VectorType> OperatorType;
    // for cpr
    typedef Dune::BCRSMatrix< Dune::FieldMatrix< double, 1, 1 > > PressureMatrixType;
    typedef Dune::BlockVector< Dune::FieldVector< double, 1 > > PressureVectorType;
    typedef Dune::MatrixAdapter<PressureMatrixType, PressureVectorType, PressureVectorType> CoarseOperatorType;
    // using Communication =  Dune::Amg::SequentialInformation;
    // using Criterion  =
    //   Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<PressureMatrixType,
    //     							Dune::Amg::FirstDiagonal> >;
    
    //constexpr int pressureVarIndex = 0;
    // using LevelTransferPolicy = Opm::PressureTransferPolicy<OperatorType,
    //     						    CoarseOperatorType,
    //     						    Communication,
    //     						    0>;
    // using CoarseSolverPolicy   =
    //   Dune::Amg::PressureSolverPolicy<CoarseOperatorType,
    //     			      LevelTransferPolicy>;
        
    
    BlockIlu0Solver(boost::property_tree::ptree prm): prm_(prm){};

    
    void solve(double* result,std::string matrixfile,std::string rhsfile){
      std::cout << matrixfile << std::endl;
      std::cout << rhsfile << std::endl;
      makeSystem(matrixfile, rhsfile);
           
      Dune::Timer perfTimer;
      perfTimer.start();
      double tol = 1e-4;
      int maxiter = 200;
      this->makeSolver(tol, maxiter);      
      int m = bz*rhs_.size();
      VectorType x(rhs_.size());
      Dune::InverseOperatorResult res;
      linsolver_->apply(x, rhs_, res);
      double time = perfTimer.stop();
      this->makeResult(result,x);
    }
    

    void solve(double* result,
	       std::vector<int>& i,
	       std::vector<int>& j,
	       std::vector<double>& val,
	       size_t rows,
	       std::vector<double>& orhs,
	       double tol,
	       int maxiter){
      this->makeSystem(i,
		       j,
		       val,
		       rows,
		       orhs);
      Dune::Timer perfTimer;
      perfTimer.start();
      //double tol = 1e-4;
      //int maxiter = 200;
      this->makeSolver(tol, maxiter);    
      int m = bz*rhs_.size();
      VectorType x(rhs_.size());
      Dune::InverseOperatorResult res;
      linsolver_->apply(x, rhs_, res);
      double time = perfTimer.stop();
      this->makeResult(result,x);
      // result is returned
    }

    void solve(VectorType& x, VectorType& rhs){
      Dune::InverseOperatorResult res;
      linsolver_->apply(x, rhs, res);
    }
    
    void makeSolver(double tol, int maxiter,MatrixType& matrix){
      matrix_ = matrix;
      makeSolver(tol, maxiter);
    }
    
    private:

    
    void makeSolver(double tol, int maxiter){
      linearoperator_.reset(new Dune::MatrixAdapter<MatrixType, VectorType, VectorType>(matrix_));
      preconditioner_ = makePreconditioner<MatrixType, VectorType>(*linearoperator_,prm_);
      int verbosity = 10;
      linsolver_.reset(new Dune::BiCGSTABSolver<VectorType>(*linearoperator_,
							 *preconditioner_,
							 tol, // desired residual reduction factor
							 maxiter, // maximum number of iterations
							 verbosity));   
    }
    
    
    
    void makeResult(double* result,VectorType& x){
      int i = 0;
      for(size_t ic = 0; ic < rhs_.size(); ic++){
        for(size_t ib = 0; ib < bz; ib++){
	  result[i] = x[ic][ib];
          i++;
	}
      }
    }
    
    void makeSystem(std::vector<int>& i,
		    std::vector<int>& j,
		    std::vector<double>& val,
		    size_t rows,
		    std::vector<double>& orhs){
      // copy rhs to block vector
      rhs_.resize(rows/bz);
      {
	int lind=0;
	for(size_t ic = 0; ic < rhs_.size(); ic++){
	  for(size_t ib = 0; ib < bz; ib++){
	    rhs_[ic][ib] = orhs[lind];
	    lind++;
	  }
	}
      }    
      
      // make block matrix
      makeMatrixMarket(matrix_,
		       i,j,val,
		       rows,
		       rows,
		       val.size());
    
    }
    
    void makeSystem(std::string matrixfile,std::string rhsfile){
      {
	std::ifstream infile(rhsfile);
	if(!infile){
	  throw std::runtime_error("Rhs file not read");
	}
	Dune::readMatrixMarket(rhs_,infile);
	//Dune::writeMatrixMarket(rhs,std::cout);
      }
      {
	std::ifstream infile(matrixfile);
	if(!infile){
	  throw std::runtime_error("Matrix file not read");
	}
	Dune::readMatrixMarket(matrix_,infile);
	std::string tmpmatrix("matrixtmp.txt");
	std::ofstream outfile(tmpmatrix);
	if(!outfile){
	  throw std::runtime_error("Could not write");
	}
	Dune::writeMatrixMarket(matrix_,outfile);
      }
    }
    

    boost::property_tree::ptree prm_;
    MatrixType matrix_;
    VectorType rhs_;
    std::shared_ptr< Dune::Preconditioner<VectorType,VectorType> > preconditioner_;
    std::shared_ptr< Dune::MatrixAdapter<MatrixType, VectorType, VectorType> > linearoperator_;
    std::shared_ptr< Dune::IterativeSolver<VectorType,VectorType> > linsolver_;
    // or cpr
    // std::shared_ptr< Dune::Preconditioner<VectorType,VectorType> > finesmoother_;
    // std::unique_ptr<CoarseSolverPolicy> coarseSolverPolicy_;
    // std::unique_ptr<LevelTransferPolicy> levelTransferPolicy_;
    // VectorType weights_;
    
    };
}
#endif /* MRST_DUNEISTL_HPP */
