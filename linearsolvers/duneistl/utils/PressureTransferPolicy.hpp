#ifndef OPM_PRESSURE_TRANSFER_POLICY_HEADER_INCLUDED
#define OPM_PRESSURE_TRANSFER_POLICY_HEADER_INCLUDED

//#include <opm/autodiff/ParallelOverlappingILU0.hpp>
//#include <opm/autodiff/FlowLinearSolverParameters.hpp>
//#include <opm/autodiff/CPRPreconditioner.hpp>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/paamg/aggregates.hh>
#include <dune/istl/paamg/twolevelmethod.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/schwarz.hh>



namespace Opm
{
template <class FineOperator, class CoarseOperator, class Communication, std::size_t VARIABLE_INDEX>
class PressureTransferPolicy : public Dune::Amg::LevelTransferPolicy<FineOperator, CoarseOperator>
{
public:
    typedef Dune::Amg::LevelTransferPolicy<FineOperator, CoarseOperator> FatherType;
    typedef Communication ParallelInformation;
    typedef typename FineOperator::domain_type FineVectorType;

public:
    PressureTransferPolicy(const Communication& comm, const FineVectorType& weights)
        : communication_(&const_cast<Communication&>(comm))
        , weights_(weights)
    {
    }

    void createCoarseLevelSystem(const FineOperator& fineOperator)
    {
        using CoarseMatrix = typename CoarseOperator::matrix_type;
        const auto& fineLevelMatrix = fineOperator.getmat();
        coarseLevelMatrix_.reset(new CoarseMatrix(fineLevelMatrix.N(), fineLevelMatrix.M(), CoarseMatrix::row_wise));
        auto createIter = coarseLevelMatrix_->createbegin();

        for (const auto& row : fineLevelMatrix) {
            for (auto col = row.begin(), cend = row.end(); col != cend; ++col) {
                createIter.insert(col.index());
            }
            ++createIter;
        }

        auto coarseRow = coarseLevelMatrix_->begin();
        for (auto row = fineLevelMatrix.begin(), rowEnd = fineLevelMatrix.end(); row != rowEnd; ++row) {
            auto coarseCol = coarseRow->begin();
            // auto& row = *rowit;
            for (auto col = row->begin(), cend = row->end(); col != cend; ++col, ++coarseCol) {
                assert(col.index() == coarseCol.index());
                double matrix_el = 0;
                auto bw = weights_[row.index()];
                for (int i = 0; i < bw.size(); ++i) {
                    matrix_el += (*col)[i][VARIABLE_INDEX] * bw[i];
                }
                *coarseCol = matrix_el;
            }
            ++coarseRow;
        }
        coarseLevelCommunication_.reset(communication_, [](Communication*) {});


        this->lhs_.resize(this->coarseLevelMatrix_->M());
        this->rhs_.resize(this->coarseLevelMatrix_->N());
        using OperatorArgs = typename Dune::Amg::ConstructionTraits<CoarseOperator>::Arguments;
        OperatorArgs oargs(*coarseLevelMatrix_, *coarseLevelCommunication_);
        this->operator_.reset(Dune::Amg::ConstructionTraits<CoarseOperator>::construct(oargs));
    }

    // compleately unsafe!!!!!!
    void calculateCoarseEntries(const FineOperator& fineOperator) // const M& fineMatrix)
    {
        const auto& fineMatrix = fineOperator.getmat();
        *coarseLevelMatrix_ = 0;
        for (auto row = fineMatrix.begin(), rowEnd = fineMatrix.end(); row != rowEnd; ++row) {
            const auto& i = row.index();
            for (auto entry = row->begin(), entryEnd = row->end(); entry != entryEnd; ++entry) {
                double matrix_el = 0;
                auto bw = weights_[i];
                for (int ii = 0; ii < bw.size(); ++ii) {
                    matrix_el += (*entry)[ii][VARIABLE_INDEX] * bw[ii];
                }
                const auto& j = entry.index();
                (*coarseLevelMatrix_)[i][j] = matrix_el;
            }
        }
    }

    void moveToCoarseLevel(const typename FatherType::FineRangeType& fine)
    {
        // Set coarse vector to zero
        this->rhs_ = 0;

        auto end = fine.end(), begin = fine.begin();

        for (auto block = begin; block != end; ++block) {
            auto bw = weights_[block.index()];
            double rhs_el = 0.0;
            for (int i = 0; i < block->size(); ++i) {
                rhs_el += (*block)[i] * bw[i];
            }
            this->rhs_[block - begin] = rhs_el;
        }


        this->lhs_ = 0;
    }

    void moveToFineLevel(typename FatherType::FineDomainType& fine)
    {

        auto end = fine.end(), begin = fine.begin();

        for (auto block = begin; block != end; ++block) {
            (*block)[VARIABLE_INDEX] = this->lhs_[block - begin];
        }
    }

    PressureTransferPolicy* clone() const { return new PressureTransferPolicy(*this); }

    const Communication& getCoarseLevelCommunication() const { return *coarseLevelCommunication_; }

private:
    Communication* communication_;
    const FineVectorType& weights_;
    std::shared_ptr<Communication> coarseLevelCommunication_;
    std::shared_ptr<typename CoarseOperator::matrix_type> coarseLevelMatrix_;
};

} // namespace Opm
#endif
