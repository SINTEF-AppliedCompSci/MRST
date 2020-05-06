/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.
*/

#ifndef MAKEPRECONDITIONER_HEADER_INCLUDED
#define MAKEPRECONDITIONER_HEADER_INCLUDED

#include "OwningTwolevelPreconditioner.hpp"

#include <dune/common/fmatrix.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixmarket.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/fastamg.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/umfpack.hh>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <cmath>
#include <iomanip>
#include <limits>
#ifdef EXTRA_PRECONDITIONERS
#include "MakeXXXPreconditioner.hpp"
#endif

namespace Dune
{

template <class MatrixType, class VectorType>
std::shared_ptr<Dune::Preconditioner<VectorType, VectorType>>
makeSeqPreconditioner(Dune::MatrixAdapter<MatrixType, VectorType, VectorType>& linearoperator,
                      boost::property_tree::ptree& prm)
{
    auto& matrix = linearoperator.getmat();
    std::shared_ptr<Dune::Preconditioner<VectorType, VectorType>> preconditioner;
    double w = prm.get<double>("w");
    ;
    int n = prm.get<int>("n");
    std::string precond(prm.get<std::string>("preconditioner"));
    if (precond == "ILU0") {
        preconditioner.reset(new Dune::SeqILU0<MatrixType, VectorType, VectorType>(matrix, w));
    } else if (precond == "Jac") {
        preconditioner.reset(new Dune::SeqJac<MatrixType, VectorType, VectorType>(matrix, n, w));
    } else if (precond == "GS") {
        preconditioner.reset(new Dune::SeqGS<MatrixType, VectorType, VectorType>(matrix, n, w));
    } else if (precond == "SOR") {
        preconditioner.reset(new Dune::SeqSOR<MatrixType, VectorType, VectorType>(matrix, n, w));
    } else if (precond == "SSOR") {
        preconditioner.reset(new Dune::SeqSSOR<MatrixType, VectorType, VectorType>(matrix, n, w));
    } else if (precond == "ILUn") {
        preconditioner.reset(new Dune::SeqILUn<MatrixType, VectorType, VectorType>(matrix, n, w));
    } else {
        std::string msg("No such seq preconditioner : ");
        msg += precond;
        throw std::runtime_error(msg);
        ;
    }
    return preconditioner;
}

template <class Smoother, class MatrixType, class VectorType>
std::shared_ptr<Dune::Preconditioner<VectorType, VectorType>>
makeAmgPreconditioner(Dune::MatrixAdapter<MatrixType, VectorType, VectorType>& linearoperator,
                      boost::property_tree::ptree& global_prm)
{
    boost::property_tree::ptree prm = global_prm.get_child("amg");
    typedef Dune::MatrixAdapter<MatrixType, VectorType, VectorType> OperatorType;
    std::shared_ptr<Dune::Preconditioner<VectorType, VectorType>> preconditioner;
    typedef Dune::Amg::AggregationCriterion<Dune::Amg::SymmetricMatrixDependency<MatrixType, Dune::Amg::FirstDiagonal>>
        CriterionBase;
    typedef Dune::Amg::CoarsenCriterion<CriterionBase> Criterion;
    int coarsenTarget = prm.get<int>("coarsenTarget");
    int ml = prm.get<int>("maxlevel");
    Criterion criterion(15, coarsenTarget);
    criterion.setDefaultValuesIsotropic(2);
    criterion.setAlpha(prm.get<double>("alpha"));
    criterion.setBeta(prm.get<double>("beta"));
    criterion.setMaxLevel(ml);
    criterion.setSkipIsolated(false);
    criterion.setDebugLevel(prm.get<int>("verbosity"));
    criterion.setNoPreSmoothSteps(prm.get<int>("pre_smooth"));
    criterion.setNoPostSmoothSteps(prm.get<int>("post_smooth"));
    Dune::Amg::Parameters parms;
    if (global_prm.get<std::string>("preconditioner") == "famg") {
        preconditioner.reset(new Dune::Amg::FastAMG<OperatorType, VectorType>(linearoperator, criterion, parms));
    } else {
        typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments SmootherArgs;
        SmootherArgs smootherArgs;
        smootherArgs.iterations = prm.get<int>("n");
        // smootherArgs.overlap=SmootherArgs::vertex;
        // smootherArgs.overlap=SmootherArgs::none;
        // smootherArgs.overlap=SmootherArgs::aggregate;
        smootherArgs.relaxationFactor = prm.get<double>("w");
        ;
        // Smoother smoother(linearoperator.getmat(), 1);// 1);
        preconditioner.reset(
            new Dune::Amg::AMG<OperatorType, VectorType, Smoother>(linearoperator, criterion, smootherArgs));
    }
    return preconditioner;
}




template <class MatrixType, class VectorType>
std::shared_ptr<Dune::Preconditioner<VectorType, VectorType>>
makeAmgPreconditioners(Dune::MatrixAdapter<MatrixType, VectorType, VectorType>& linearoperator,
                       boost::property_tree::ptree& prm)
{
    std::shared_ptr<Dune::Preconditioner<VectorType, VectorType>> preconditioner;
    if (prm.get<std::string>("preconditioner") == "famg") {
        // smoother type should not be used
        preconditioner
            = makeAmgPreconditioner<Dune::SeqILU0<MatrixType, VectorType, VectorType>, MatrixType, VectorType>(
                linearoperator, prm);
        return preconditioner;
    }

    std::string precond = prm.get<std::string>("amg.smoother");
    if (precond == "ILU0") {
        preconditioner
            = makeAmgPreconditioner<Dune::SeqILU0<MatrixType, VectorType, VectorType>, MatrixType, VectorType>(
                linearoperator, prm);
    } else if (precond == "Jac") {
        preconditioner
            = makeAmgPreconditioner<Dune::SeqJac<MatrixType, VectorType, VectorType>, MatrixType, VectorType>(
                linearoperator, prm);
        // }else if(precond == "GS"){
        //   preconditioner = makeAmgPreconditioner<
        // 	Dune::SeqGS<MatrixType, VectorType, VectorType>,
        // 	MatrixType, VectorType>(linearoperator,prm);
    } else if (precond == "SOR") {
        preconditioner
            = makeAmgPreconditioner<Dune::SeqSOR<MatrixType, VectorType, VectorType>, MatrixType, VectorType>(
                linearoperator, prm);
    } else if (precond == "SSOR") {
        preconditioner
            = makeAmgPreconditioner<Dune::SeqSSOR<MatrixType, VectorType, VectorType>, MatrixType, VectorType>(
                linearoperator, prm);
    } else if (precond == "ILUn") {
        preconditioner
            = makeAmgPreconditioner<Dune::SeqILUn<MatrixType, VectorType, VectorType>, MatrixType, VectorType>(
                linearoperator, prm);
    } else {
        std::string msg("No such seq preconditioner : ");
        msg += precond;
        throw std::runtime_error(msg);
        ;
    }
    return preconditioner;
}

template <class MatrixType, class VectorType, int bz>
std::shared_ptr<Dune::Preconditioner<VectorType, VectorType>>
makeTwoLevelPreconditioner(Dune::MatrixAdapter<MatrixType, VectorType, VectorType>& linearoperator,
                           boost::property_tree::ptree& global_prm)
{
    boost::property_tree::ptree prm = global_prm.get_child("cpr");
    std::shared_ptr<Dune::Preconditioner<VectorType, VectorType>> preconditioner;
    if (global_prm.get<std::string>("preconditioner") == "cpr") {
        preconditioner.reset(new OwningTwoLevelPreconditioner<bz>(linearoperator, prm));
    } else if (global_prm.get<std::string>("preconditioner") == "cprt") {
        preconditioner.reset(new OwningTwoLevelPreconditionerTranspose<bz>(linearoperator, prm));
    } else {
        std::string msg("Wrong cpr Should not happen");
        throw std::runtime_error(msg);
    }
    return preconditioner;
}

template <class MatrixType, class VectorType, int bz>
std::shared_ptr<Dune::Preconditioner<VectorType, VectorType>>
makePreconditioner(Dune::MatrixAdapter<MatrixType, VectorType, VectorType>& linearoperator,
                   boost::property_tree::ptree& prm)
{
    if ((prm.get<std::string>("preconditioner") == "famg") or (prm.get<std::string>("preconditioner") == "amg")) {
        return makeAmgPreconditioners<MatrixType, VectorType>(linearoperator, prm);
    } else if ((prm.get<std::string>("preconditioner") == "cpr")
               or (prm.get<std::string>("preconditioner") == "cprt")) {
        return makeTwoLevelPreconditioner<MatrixType, VectorType, bz>(linearoperator, prm);
#ifdef EXTRA_PRECONDITIONERS	
    } else if (	prm.get<std::string>("preconditioner") == "xxx"){
      return XXX::makePreconditioner<bz>(linearoperator, prm);
#endif
    } else {
        return makeSeqPreconditioner<MatrixType, VectorType>(linearoperator, prm);
    }
}


} // namespace Dune


#endif // MAKEPRECONDITIONER_HEADER_INCLUDED
