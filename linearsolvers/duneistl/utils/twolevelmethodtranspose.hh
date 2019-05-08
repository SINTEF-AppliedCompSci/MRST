// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ISTL_TWOLEVELMETHOD_TRANSPOSE_HH
#define DUNE_ISTL_TWOLEVELMETHOD_TRANSPOSE_HH


#include<dune/istl/operators.hh>
#include<dune/istl/solver.hh>
#include<dune/common/unused.hh>


/**
 * @addtogroup ISTL_PAAMG
 * @{
 * @file
 * @author Markus Blatt
 * @brief Algebraic twolevel methods.
 */
namespace Dune
{
  namespace Amg
  {


    /**
     * @tparam FO The type of the fine level linear operator.
     * @tparam CSP The type of the coarse level solver policy.
     * @tparam S The type of the fine level smoother used.
     */
    template<class FO, class CSP, class S>
    class TwoLevelMethodTranspose :
      public Preconditioner<typename FO::domain_type, typename FO::range_type>
    {
    public:
      /** @brief The type of the policy for constructing the coarse level solver. */
      typedef CSP CoarseLevelSolverPolicy;
      /** @brief The type of the coarse level solver. */
      typedef typename CoarseLevelSolverPolicy::CoarseLevelSolver CoarseLevelSolver;
      /**
       * @brief The linear operator of the finel level system. Has to be
       * derived from AssembledLinearOperator.
       */
      typedef FO FineOperatorType;
      /**
       * @brief The type of the range of the fine level operator.
       */
      typedef typename FineOperatorType::range_type FineRangeType;
      /**
       * @brief The type of the domain of the fine level operator.
       */
      typedef typename FineOperatorType::domain_type FineDomainType;
      /**
       * @brief The linear operator of the finel level system. Has to be
       * derived from AssembledLinearOperator.
       */
      typedef typename CSP::Operator CoarseOperatorType;
      /**
       * @brief The type of the range of the coarse level operator.
       */
      typedef typename CoarseOperatorType::range_type CoarseRangeType;
      /**
       * @brief The type of the domain of the coarse level operator.
       */
      typedef typename CoarseOperatorType::domain_type CoarseDomainType;
      /**
       * @brief The type of the fine level smoother.
       */
      typedef S SmootherType;

      /**
       * @brief Constructs a two level method.
       *
       * @tparam CoarseSolverPolicy The policy for constructing the coarse
       * solver, e.g. OneStepAMGCoarseSolverPolicy
       * @param op The fine level operator.
       * @param smoother The fine level smoother.
       * @param policy The level transfer policy.
       * @param coarsePolicy The policy for constructing the coarse level solver.
       * @param preSteps The number of smoothing steps to apply before the coarse
       * level correction.
       * @param preSteps The number of smoothing steps to apply after the coarse
       * level correction.
       */
      TwoLevelMethodTranspose(const FineOperatorType& op,
                              std::shared_ptr<SmootherType> smoother,
                              const LevelTransferPolicy<FineOperatorType,
                              CoarseOperatorType>& policy,
                              CoarseLevelSolverPolicy& coarsePolicy,
                              std::size_t preSteps=1, std::size_t postSteps=1)
        : operator_(&op), smoother_(smoother),
          preSteps_(preSteps), postSteps_(postSteps)
      {
        policy_ = policy.clone();
        policy_->createCoarseLevelSystem(*operator_);
        coarseSolver_=coarsePolicy.createCoarseLevelSolver(*policy_);
      }

      TwoLevelMethodTranspose(const TwoLevelMethodTranspose& other)
        : operator_(other.operator_), coarseSolver_(new CoarseLevelSolver(*other.coarseSolver_)),
          smoother_(other.smoother_), policy_(other.policy_->clone()),
          preSteps_(other.preSteps_), postSteps_(other.postSteps_)
      {}

      ~TwoLevelMethodTranspose()
      {
        // Each instance has its own policy.
        delete policy_;
        delete coarseSolver_;
      }

      void pre(FineDomainType& x, FineRangeType& b)
      {
        smoother_->pre(x,b);
      }

      void post(FineDomainType& x)
      {
        DUNE_UNUSED_PARAMETER(x);    
      }

      void apply(FineDomainType& v, const FineRangeType& d)
      {
        FineDomainType u(v);
        FineRangeType rhs(d);
        LevelContext context;
        SequentialInformation info;
        context.pinfo=&info;
        context.lhs=&u;
        context.update=&v;
        context.smoother=smoother_;
        context.rhs=&rhs;
        context.matrix=operator_;
        // Presmoothing
        //presmooth(context, preSteps_);
        //Coarse grid correction
        policy_->moveToCoarseLevel(*context.rhs);
        InverseOperatorResult res;
        coarseSolver_->apply(policy_->getCoarseLevelLhs(), policy_->getCoarseLevelRhs(), res);
        *context.lhs=0;
        policy_->moveToFineLevel(*context.lhs);
        //presmooth(context, preSteps_);
        *context.update += *context.lhs;
        // Postsmoothing
        {
          //                                  matrix*rhs
          auto& levelContext=context;
          // update defect
          SmootherApplier<typename LevelContext::SmootherType>
            ::postSmooth(*levelContext.smoother, *levelContext.lhs, *levelContext.rhs);

          levelContext.matrix->applyscaleadd(-1, *levelContext.lhs,
                                             *levelContext.rhs);
          *levelContext.lhs=0;
          levelContext.pinfo->project(*levelContext.rhs);
          
          // Accumulate update
          *levelContext.update += *levelContext.lhs;
        }      
        //postsmooth(context, postSteps_);

      }

      //! Category of the preconditioner (see SolverCategory::Category)
      virtual SolverCategory::Category category() const
      {
        return SolverCategory::sequential;
      }

    private:
      /**
       * @brief Struct containing the level information.
       */
      struct LevelContext
      {
        /** @brief The type of the smoother used. */
        typedef S SmootherType;
        /** @brief A pointer to the smoother. */
        std::shared_ptr<SmootherType> smoother;
        /** @brief The left hand side passed to the and returned by the smoother. */
        FineDomainType* lhs;
        /*
         * @brief The right hand side holding the current residual.
         *
         * This is passed to the smoother as the right hand side.
         */
        FineRangeType* rhs;
        /**
         * @brief The total update calculated by the preconditioner.
         *
         * I.e. all update from smoothing and coarse grid correction summed up.
         */
        FineDomainType* update;
        /** @parallel information */
        SequentialInformation* pinfo;
        /**
         * @brief The matrix that we are solving.
         *
         * Needed to update the residual.
         */
        const FineOperatorType* matrix;
      };
      const FineOperatorType* operator_;
      /** @brief The coarse level solver. */
      CoarseLevelSolver* coarseSolver_;
      /** @brief The fine level smoother. */
      std::shared_ptr<S> smoother_;
      /** @brief Policy for prolongation, restriction, and coarse level system creation. */
      LevelTransferPolicy<FO,typename CSP::Operator>* policy_;
      /** @brief The number of presmoothing steps to apply. */
      std::size_t preSteps_;
      /** @brief The number of postsmoothing steps to apply. */
      std::size_t postSteps_;
    };
  }// end namespace Amg
}// end namespace Dune

/** @} */
#endif
