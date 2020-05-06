#ifndef OPM_PRESSURE_SOLVER_POLICY_HEADER_INCLUDED
#define OPM_PRESSURE_SOLVER_POLICY_HEADER_INCLUDED

#include <boost/property_tree/ptree.hpp>

#include <dune/istl/solver.hh>

namespace Dune
{
namespace Amg
{
    namespace pt = boost::property_tree;

    template <class O, class P, class S>
    class PressureSolverPolicy
    {
    public:
        typedef P LevelTransferPolicy;
        /** @brief The type of the linear operator used. */
        typedef O Operator;
        /** @brief The type of the range and domain of the operator. */
        typedef typename O::range_type X;
        /** @brief The type of the solver for the pressure system. */
        typedef S Solver;
        /**
         * @brief Constructs the coarse solver policy.
         * @param args The arguments used for constructing the smoother.
         * @param c The crition used for the aggregation within AMG.
         */
        PressureSolverPolicy(const pt::ptree prm)
            : prm_(prm)
        {
        }
        /** @brief Copy constructor. */
        PressureSolverPolicy(const PressureSolverPolicy& other)
            : coarseOperator_(other.coarseOperator_)
            , prm_(prm_)
        {
        }

    private:
        /**
         * @brief A wrapper that makes an inverse operator out of AMG.
         *
         * The operator will use one step of AMG to approximately solve
         * the coarse level system.
         */
        struct PressureInverseOperator : public Dune::InverseOperator<X, X> {
            PressureInverseOperator(Operator& op, const boost::property_tree::ptree& prm)
                : linsolver_()
                , op_(op) //, prm_(prm)
            {
                linsolver_.reset(new Solver(prm));
                int maxiter = prm.get<int>("maxiter");
                double tol = prm.get<double>("tol");
                linsolver_->makeSolver(tol, maxiter, op_.getmat());
            }

            Dune::SolverCategory::Category category() const override
            {
                return Dune::SolverCategory::sequential;
            }

            void apply(X& x, X& b, double reduction, Dune::InverseOperatorResult& res) override
            {
                linsolver_->apply(x, b, res);
            }

            void apply(X& x, X& b, Dune::InverseOperatorResult& res) override
            {
                linsolver_->apply(x, b, res);
            }

        private:
            std::shared_ptr<Solver> linsolver_;
            Operator& op_;
        };

    public:
        /** @brief The type of solver constructed for the coarse level. */
        typedef PressureInverseOperator CoarseLevelSolver;

        /**
         * @brief Constructs a coarse level solver.
         *
         * @param transferPolicy The policy describing the transfer between levels.
         * @return A pointer to the constructed coarse level solver.
         */
        template <class LTP>
        void setCoarseOperator(LTP& transferPolicy)
        {
            coarseOperator_ = transferPolicy.getCoarseLevelOperator();
        }
        template <class LTP>
        CoarseLevelSolver* createCoarseLevelSolver(LTP& transferPolicy)
        {
            coarseOperator_ = transferPolicy.getCoarseLevelOperator();
            //const LevelTransferPolicy& transfer = reinterpret_cast<const LevelTransferPolicy&>(transferPolicy);
            PressureInverseOperator* inv = new PressureInverseOperator(*coarseOperator_, prm_);

            return inv; // std::shared_ptr<InverseOperator<X,X> >(inv);
        }

    private:
        /** @brief The coarse level operator. */
        std::shared_ptr<Operator> coarseOperator_;
        pt::ptree prm_;
    };
} // namespace Amg
} // namespace Dune


#endif
