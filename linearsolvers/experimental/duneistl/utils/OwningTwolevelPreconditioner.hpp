/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.
*/

#ifndef OWNINGTWOLEVELPRECONDITIONER_HEADER_INCLUDED
#define OWNINGTWOLEVELPRECONDITIONER_HEADER_INCLUDED

#include "GetQuasiImpesWeights.hpp"
#include "PressureSolverPolicy.hpp"
#include "PressureTransferPolicy.hpp"
#include "PressureTransferPolicyTranspose.hpp"

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/paamg/amg.hh>

#include <boost/property_tree/ptree.hpp>

#include <cmath>


namespace Dune
{

// Circular dependency between makePreconditioner() [which can make an OwningTwolevelPreconditioner]
// and OwningTwoLevelPreconditioner [which uses makePreconditioner() to choose the fine-level smoother]
// must be broken, accomplished by forward-declaration here.
template <class MatrixType, class VectorType, int bz>
std::shared_ptr<Dune::Preconditioner<VectorType, VectorType>>
makePreconditioner(Dune::MatrixAdapter<MatrixType, VectorType, VectorType>& linearoperator,
                   boost::property_tree::ptree& prm);

// Must forward-declare FlexibleSolver as we want to use it as solver for the pressure system.
template <int bz>
class FlexibleSolver;

template <int bz>
class OwningTwoLevelPreconditioner : public Dune::Preconditioner<Dune::BlockVector<Dune::FieldVector<double, bz>>,
                                                                 Dune::BlockVector<Dune::FieldVector<double, bz>>>
{
public:
    typedef boost::property_tree::ptree pt;
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double, bz, bz>> MatrixType;
    typedef Dune::BlockVector<Dune::FieldVector<double, bz>> VectorType;
    typedef Dune::MatrixAdapter<MatrixType, VectorType, VectorType> OperatorType;

    OwningTwoLevelPreconditioner(OperatorType& linearoperator, pt& prm)
        : finesmoother_(makePreconditioner<MatrixType, VectorType, bz>(linearoperator, prm.get_child("finesmoother")))
        , comm_()
        , weights_(
		   Opm::Amg::getQuasiImpesWeights<MatrixType, VectorType>(linearoperator.getmat(), prm.get<int>("pressure_var_index"))) // TODO
        , levelTransferPolicy_(comm_, weights_, prm.get<int>("pressure_var_index"))
        , coarseSolverPolicy_(prm.get_child("coarsesolver"))
        , twolevel_method_(linearoperator,
			   finesmoother_,
			   levelTransferPolicy_,
			   coarseSolverPolicy_,
			   prm.get<int>("pre_smooth"),
			   prm.get<int>("post_smooth"))
    {
      if(prm.get<int>("verbosity")>10){
	std::ofstream outfile("weight_cpr.txt");
	if(!outfile){
	  throw std::runtime_error("Could not write weights");
	}
	Dune::writeMatrixMarket(weights_,outfile);
      }
    }

    virtual void pre(VectorType& x, VectorType& b) override
    {
        twolevel_method_.pre(x, b);
    }
    virtual void apply(VectorType& v, const VectorType& d) override
    {
        twolevel_method_.apply(v, d);
    }
    virtual void post(VectorType& x) override
    {
        twolevel_method_.post(x);
    }
    virtual Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::sequential;
    }

private:
    // for cpr
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>> PressureMatrixType;
    typedef Dune::BlockVector<Dune::FieldVector<double, 1>> PressureVectorType;
    typedef Dune::MatrixAdapter<PressureMatrixType, PressureVectorType, PressureVectorType> CoarseOperatorType;
    using Communication = Dune::Amg::SequentialInformation;
    constexpr static int pressureVarIndex = 1;
    using LevelTransferPolicy
    = Opm::PressureTransferPolicy<OperatorType, CoarseOperatorType, Communication>;// pressureVarIndex>;
    using CoarseSolverPolicy
        = Dune::Amg::PressureSolverPolicy<CoarseOperatorType, LevelTransferPolicy, FlexibleSolver<1>>;
    using TwoLevelMethod
        = Dune::Amg::TwoLevelMethod<OperatorType, CoarseSolverPolicy, Dune::Preconditioner<VectorType, VectorType>>;

    std::shared_ptr<Dune::Preconditioner<VectorType, VectorType>> finesmoother_;
    Communication comm_;
    VectorType weights_;
    LevelTransferPolicy levelTransferPolicy_;
    CoarseSolverPolicy coarseSolverPolicy_;
    TwoLevelMethod twolevel_method_;
};

template <int bz>
class OwningTwoLevelPreconditionerTranspose
    : public Dune::Preconditioner<Dune::BlockVector<Dune::FieldVector<double, bz>>,
                                  Dune::BlockVector<Dune::FieldVector<double, bz>>>
{
public:
    typedef boost::property_tree::ptree pt;
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double, bz, bz>> MatrixType;
    typedef Dune::BlockVector<Dune::FieldVector<double, bz>> VectorType;
    typedef Dune::MatrixAdapter<MatrixType, VectorType, VectorType> OperatorType;

    OwningTwoLevelPreconditionerTranspose(OperatorType& linearoperator, pt& prm)
        : finesmoother_(makePreconditioner<MatrixType, VectorType, bz>(linearoperator, prm.get_child("finesmoother")))
        , comm_()
        , weights_(
              Opm::Amg::getQuasiImpesWeights<MatrixType, VectorType>(linearoperator.getmat(), prm.get<int>("pressure_var_index"), true))
        , levelTransferPolicy_(comm_, weights_, prm.get<int>("pressure_var_index"))
        , coarseSolverPolicy_(prm.get_child("coarsesolver"))
        , twolevel_method_(linearoperator, finesmoother_, levelTransferPolicy_, coarseSolverPolicy_, 1, 0)
    {
      if(prm.get<int>("verbosity")>10){
	std::ofstream outfile("weight_cprt.txt");
	if(!outfile){
	  throw std::runtime_error("Could not write weights");
	}
	Dune::writeMatrixMarket(weights_,outfile);
      }
    }

    virtual void pre(VectorType& x, VectorType& b) override
    {
        twolevel_method_.pre(x, b);
    }
    virtual void apply(VectorType& v, const VectorType& d) override
    {
        twolevel_method_.apply(v, d);
    }
    virtual void post(VectorType& x) override
    {
        twolevel_method_.post(x);
    }
    virtual Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::sequential;
    }

private:
    // for cpr
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>> PressureMatrixType;
    typedef Dune::BlockVector<Dune::FieldVector<double, 1>> PressureVectorType;
    typedef Dune::MatrixAdapter<PressureMatrixType, PressureVectorType, PressureVectorType> CoarseOperatorType;
    using Communication = Dune::Amg::SequentialInformation;
  //constexpr static int pressureVarIndex = 1;
    using LevelTransferPolicy
    = Opm::PressureTransferPolicyTranspose<OperatorType, CoarseOperatorType, Communication>;// pressureVarIndex>;
    using CoarseSolverPolicy
        = Dune::Amg::PressureSolverPolicy<CoarseOperatorType, LevelTransferPolicy, FlexibleSolver<1>>;
    using TwoLevelMethod = Dune::Amg::
        TwoLevelMethod<OperatorType, CoarseSolverPolicy, Dune::Preconditioner<VectorType, VectorType>>;

    std::shared_ptr<Dune::Preconditioner<VectorType, VectorType>> finesmoother_;
    Communication comm_;
    VectorType weights_;
    LevelTransferPolicy levelTransferPolicy_;
    CoarseSolverPolicy coarseSolverPolicy_;
    TwoLevelMethod twolevel_method_;
};



} // namespace Dune




#endif // OWNINGTWOLEVELPRECONDITIONER_HEADER_INCLUDED
