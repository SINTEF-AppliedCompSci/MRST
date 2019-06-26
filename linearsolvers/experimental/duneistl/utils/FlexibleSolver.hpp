#ifndef FLEXIBLESOLVER_HPP
#define FLEXIBLESOLVER_HPP



#include <fstream>
#include <iostream>
/* MEX interfaces */
/* Block system support */
#include <cmath>
#include <iomanip>
#include <limits>
//#include <boost/program_options.hpp>
#include <dune/common/fmatrix.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixmarket.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/umfpack.hh>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/fastamg.hh>
#include "makePreconditioner.hpp"
namespace Dune
{

    
  template <class Preconditioner>
  class Preconditoner2InverseOperator :
    public Dune::InverseOperator<typename Preconditioner::domain_type,
				 typename Preconditioner::range_type>
  {
  public:
    typedef typename Preconditioner::domain_type domain_type;
    //! \brief The range type of the preconditioner.
    typedef typename Preconditioner::range_type range_type;
    typedef typename Preconditioner::field_type field_type;
    using X = domain_type;
    using Y = range_type;
    Preconditoner2InverseOperator(Preconditioner& preconditioner):
      preconditioner_(preconditioner)
    {
    }
    
    
    virtual void apply (X& v, Y& d, InverseOperatorResult& res){
      preconditioner_.pre(v,d);
      preconditioner_.apply(v,d);
      preconditioner_.post(v);
    }
    
    virtual void apply (X& v, Y& d,double reduction, InverseOperatorResult& res){
      this->apply(v,d, res);
    }

    virtual SolverCategory::Category category() const{
      return preconditioner_.category();
    }
  private:
    Preconditioner& preconditioner_;  

  };


  


namespace MatrixMarketImpl
{
    template <typename T, typename A, int brows, int bcols, typename D>
    void makeSparseEntries(Dune::BCRSMatrix<Dune::FieldMatrix<T, brows, bcols>, A>& matrix,
                           std::vector<int>& i,
                           std::vector<int>& j,
                           std::vector<double>& val,
                           std::size_t entries,
                           const D&)
    {
        // addapted from readSparseEntries
        typedef Dune::BCRSMatrix<Dune::FieldMatrix<T, brows, bcols>, A> Matrix;
        std::vector<std::set<IndexData<D>>> rows(matrix.N() * brows);
        for (size_t kk = 0; kk < i.size(); ++kk) {
            std::size_t row;
            IndexData<D> data;
            row = i[kk];
            assert(row / bcols < matrix.N());
            data.number = val[kk];
            data.index = j[kk];
            assert(data.index / bcols < matrix.M());
            rows[row].insert(data);
        }

        // Setup the matrix sparsity pattern
        int nnz = 0;
        for (typename Matrix::CreateIterator iter = matrix.createbegin(); iter != matrix.createend(); ++iter) {
            for (std::size_t brow = iter.index() * brows, browend = iter.index() * brows + brows; brow < browend;
                 ++brow) {
                typedef typename std::set<IndexData<D>>::const_iterator Siter;
                for (Siter siter = rows[brow].begin(), send = rows[brow].end(); siter != send; ++siter, ++nnz)
                    iter.insert(siter->index / bcols);
            }
        }

        // Set the matrix values
        matrix = 0;

        MatrixValuesSetter<D, brows, bcols> Setter;

        Setter(rows, matrix);
    }
} // end namespace MatrixMarketImpl

template <typename T, typename A, int brows, int bcols>
void
makeMatrixMarket(Dune::BCRSMatrix<Dune::FieldMatrix<T, brows, bcols>, A>& matrix,
                 std::vector<int> i,
                 std::vector<int> j,
                 std::vector<T> val,
                 size_t rows,
                 size_t cols,
                 size_t entries)
{
    // addapted from readMatrixMarket
    using namespace MatrixMarketImpl;
    // std::size_t rows, cols, entries;
    // std::size_t nnz, blockrows, blockcols;
    // std::tie(blockrows, blockcols, nnz) = calculateNNZ<brows, bcols>(rows, cols, entries, header);
    std::size_t blockrows = rows / brows;
    std::size_t blockcols = cols / bcols;
    //std::size_t blocksize = brows * bcols;
    //std::size_t blockentries = 0;
    //blockentries = entries / blocksize; // nnz
    matrix.setSize(blockrows, blockcols);
    matrix.setBuildMode(Dune::BCRSMatrix<Dune::FieldMatrix<T, brows, bcols>, A>::row_wise);
    makeSparseEntries(matrix, i, j, val, entries, NumericWrapper<T>());
}



template <int bz>
class FlexibleSolver
{
public:
    typedef boost::property_tree::ptree pt;
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double, bz, bz>> MatrixType;
    typedef Dune::BlockVector<Dune::FieldVector<double, bz>> VectorType;
    typedef Dune::MatrixAdapter<MatrixType, VectorType, VectorType> OperatorType;
    // for cpr
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>> PressureMatrixType;
    typedef Dune::BlockVector<Dune::FieldVector<double, 1>> PressureVectorType;
    typedef Dune::MatrixAdapter<PressureMatrixType, PressureVectorType, PressureVectorType> CoarseOperatorType;

    explicit FlexibleSolver(boost::property_tree::ptree prm)
        : prm_(prm) {};


    void solve(double* result, std::string matrixfile, std::string rhsfile)
    {
        std::cout << matrixfile << std::endl;
        std::cout << rhsfile << std::endl;
        makeSystem(matrixfile, rhsfile);

        Dune::Timer perfTimer;
        perfTimer.start();
        double tol = 1e-4;
        int maxiter = 200;
        this->makeSolver(tol, maxiter);
        int m = bz * rhs_.size();
        VectorType x(rhs_.size());
        Dune::InverseOperatorResult res;
        linsolver_->apply(x, rhs_, res);
        double time = perfTimer.stop();
        this->makeResult(result, x);
    }


    void solve(double* result,
               std::vector<int>& i,
               std::vector<int>& j,
               std::vector<double>& val,
               size_t rows,
               std::vector<double>& orhs,
               double tol,
               int maxiter,
               boost::property_tree::ptree& out)
    {
        Dune::Timer buildTimer;
        this->makeSystem(i, j, val, rows, orhs);
        double btime = buildTimer.stop();
        out.put("time.matrixtime", btime);
        std::cout << "Dune build matrix time " << btime << std::endl;
        Dune::Timer perfTimer;
        perfTimer.start();
        // double tol = 1e-4;
        // int maxiter = 200;
        Dune::Timer sperfTimer;

        this->makeSolver(tol, maxiter);
        double stime = sperfTimer.stop();
        out.put("time.setup", stime);
        std::cout << "Dune setup time " << stime << std::endl;
        //int m = bz * rhs_.size();
        VectorType x(rhs_.size());
        Dune::InverseOperatorResult res;
        linsolver_->apply(x, rhs_, res);
        double time = perfTimer.stop();
        this->makeResult(result, x);
        std::cout << "Dune Solver time " << time << std::endl;
        out.put("time.solvetime", time - stime);
        out.put("time.solvetotal", time);
        out.put("res.iterations", res.iterations);
        out.put("res.reduction", res.reduction);
        out.put("res.conv_rate", res.conv_rate);
        out.put("res.elapsed", res.elapsed);

        // result is returned
    }

    void solve(VectorType& x, VectorType& rhs)
    {
        Dune::InverseOperatorResult res;
        linsolver_->apply(x, rhs, res);
    }

    void apply(VectorType& x, VectorType& rhs, Dune::InverseOperatorResult& res)
    {
        linsolver_->apply(x, rhs, res);
    }

    void makeSolver(double tol, int maxiter, const MatrixType& matrix)
    {
        matrix_ = matrix;
        makeSolver(tol, maxiter);
    }

private:
    void makeSolver(double tol, int maxiter)
    {
        linearoperator_.reset(new Dune::MatrixAdapter<MatrixType, VectorType, VectorType>(matrix_));
        preconditioner_ = Dune::makePreconditioner<MatrixType, VectorType, bz>(*linearoperator_, prm_);
        int verbosity = prm_.get<int>("verbosity");
        std::string solver_type = prm_.get<std::string>("solver");
        if (solver_type == "bicgstab") {
            linsolver_.reset(new Dune::BiCGSTABSolver<VectorType>(*linearoperator_,
                                                                  *preconditioner_,
                                                                  tol, // desired residual reduction factor
                                                                  maxiter, // maximum number of iterations
                                                                  verbosity));
        } else if (solver_type == "loopsolver") {
            linsolver_.reset(new Dune::LoopSolver<VectorType>(*linearoperator_,
                                                              *preconditioner_,
                                                              tol, // desired residual reduction factor
                                                              maxiter, // maximum number of iterations
                                                              verbosity));
        } else if (solver_type == "gmres") {
            int restart = prm_.get<int>("restart");
            linsolver_.reset(new Dune::RestartedGMResSolver<VectorType>(*linearoperator_,
                                                                        *preconditioner_,
                                                                        tol,
                                                                        restart, // desired residual reduction factor
                                                                        maxiter, // maximum number of iterations
                                                                        verbosity));
#if HAVE_SUITESPARSE_UMFPACK
 	} else if (solver_type == "umfpack"){
 	  bool dummy = false;
 	  linsolver_.reset(new Dune::UMFPack<MatrixType>(linearoperator_->getmat(),verbosity, dummy) );
#endif
	} else if (solver_type == "precond"){
	  
	  linsolver_.reset(new
			   Dune::Preconditoner2InverseOperator<
			   Dune::Preconditioner<VectorType, VectorType>
			   >(*preconditioner_));
	  
        } else {
            std::string msg("Solver not known ");
            msg += solver_type;
            throw std::runtime_error(msg);
        }
    }
 


    void makeResult(double* result, VectorType& x)
    {
        size_t i = 0;
        for (size_t ic = 0; ic < rhs_.size(); ic++) {
            for (size_t ib = 0; ib < bz; ib++) {
                result[i] = x[ic][ib];
                i++;
            }
        }
    }

    void makeSystem(
        std::vector<int>& i, std::vector<int>& j, std::vector<double>& val, size_t rows, std::vector<double>& orhs)
    {
        // copy rhs to block vector
        rhs_.resize(rows / bz);
        {
            size_t lind = 0;
            for (size_t ic = 0; ic < rhs_.size(); ic++) {
                for (size_t ib = 0; ib < bz; ib++) {
                    rhs_[ic][ib] = orhs[lind];
                    lind++;
                }
            }
        }

        // make block matrix
        makeMatrixMarket(matrix_, i, j, val, rows, rows, val.size());
    }

    void makeSystem(std::string matrixfile, std::string rhsfile)
    {
        {
            std::ifstream infile(rhsfile);
            if (!infile) {
                throw std::runtime_error("Rhs file not read");
            }
            Dune::readMatrixMarket(rhs_, infile);
            // Dune::writeMatrixMarket(rhs,std::cout);
        }
        {
            std::ifstream infile(matrixfile);
            if (!infile) {
                throw std::runtime_error("Matrix file not read");
            }
            Dune::readMatrixMarket(matrix_, infile);
            std::string tmpmatrix("matrixtmp.txt");
            std::ofstream outfile(tmpmatrix);
            if (!outfile) {
                throw std::runtime_error("Could not write");
            }
            Dune::writeMatrixMarket(matrix_, outfile);
        }
    }


    boost::property_tree::ptree prm_;
    MatrixType matrix_;
    VectorType rhs_;
    std::shared_ptr<Dune::Preconditioner<VectorType, VectorType>> preconditioner_;
    std::shared_ptr<Dune::MatrixAdapter<MatrixType, VectorType, VectorType>> linearoperator_;
  //std::shared_ptr<Dune::IterativeSolver<VectorType, VectorType>> linsolver_;
    std::shared_ptr<Dune::InverseOperator<VectorType, VectorType>> linsolver_;
};

} // namespace Dune



#endif /* FLEXIBLESOLVER_HPP */
