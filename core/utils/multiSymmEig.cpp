#include <mex.h>
#include <lapack.h>

#include <algorithm>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP

#include <omp.h>

#else  // ! defined(_OPENMP)

#define omp_get_thread_num()  0
#define omp_get_num_threads() 1

#endif  // _OPENMP

class EVProblem
{
public:
    bool solve(const mwSignedIndex n, const double* A)
    {
        resize(n);

        A_.assign(A, n, n);

        return solve();
    }

    void eigenValues(double* pr) const
    {
        const auto& R = E_.Lambda_.Data();
        std::copy(R.begin(), R.end(), pr);
    }

    void eigenVectors(double* ev) const
    {
        const auto& V = E_.Vectors_.Data();
        std::copy(V.begin(), V.end(), ev);
    }

private:
    template <typename T = double>
    class Array {
    public:
        void resize(const mwSignedIndex nRows,
                    const mwSignedIndex nCols = 1)
        {
            x_.resize(nRows * nCols, T(0));
            ld_ = nRows;
        }

        void assign(const T*            x,
                    const mwSignedIndex nRows,
                    const mwSignedIndex nCols = 1)
        {
            std::copy(x, x + (nRows * nCols), x_.begin());
            ld_ = nRows;
        }

        const std::vector<T>& Data() const { return x_; }

        T*            Data()       { return x_.data(); }
        mwSignedIndex LD  () const { return ld_; }

    private:
        std::vector<T> x_;
        mwSignedIndex  ld_;
    };

    struct EigenResults {
        Array<> Lambda_;
        Array<> Vectors_;

        void resize(const mwSignedIndex n)
        {
            Lambda_ .resize(n);
            Vectors_.resize(n, n);
        }
    };

    struct ProblemCharacteristics {
        const char job_z_ { 'V' }; // Do compute eigenvectors
        const char range_ { 'A' }; // Compute all eigenvalues
        const char uplo_  { 'L' }; // Reference lower triangle

        const double vl_ { - std::numeric_limits<double>::max() };
        const double vu_ {   std::numeric_limits<double>::max() };

        // Work arrays
        Array<>              Work_;
        Array<mwSignedIndex> IWork_;
        Array<mwSignedIndex> ISuppZ_;

        void resize(const mwSignedIndex n,
                    const mwSignedIndex lwork,
                    const mwSignedIndex liwork)
        {
            Work_  .resize(lwork);
            IWork_ .resize(liwork);
            ISuppZ_.resize(2 * n);
        }
    };

    // System matrix.
    Array<> A_;

    // Numerical output
    EigenResults E_;

    //
    ProblemCharacteristics C_;

    bool resize(const mwSignedIndex n)
    {
        A_.resize(n, n);
        E_.resize(n);

        mwSignedIndex LWork, LIWork;
        if (optimalWorkSize(n, LWork, LIWork)) {
            C_.resize(n, LWork, LIWork);

            return true;
        }

        return false;
    }

    bool optimalWorkSize(const mwSignedIndex n,
                         mwSignedIndex&      lwork,
                         mwSignedIndex&      liwork) const
    {
        mwSignedIndex N = n;

        auto JobZ  = C_.job_z_;
        auto Range = C_.range_;
        auto UpLo  = C_.uplo_;

        double* A = nullptr;  mwSignedIndex lda = n;

        auto vl = C_.vl_;
        auto vu = C_.vu_;

        mwSignedIndex il = 1;
        mwSignedIndex iu = n;

        double abstol;
        {
            char param[] = "Safe minimum";
            abstol = dlamch(param);
        }

        mwSignedIndex m;

        double* w = nullptr;
        double* z = nullptr;
        mwSignedIndex ldz = n;

        mwSignedIndex* isuppz = nullptr;

        double        work [1];
        mwSignedIndex iwork[1];

        mwSignedIndex info = 0;

        lwork  = -1;
        liwork = -1;

        dsyevr(&JobZ, &Range, &UpLo,
               &N, A, &lda, &vl, &vu, &il, &iu,
               &abstol, &m, w, z, &ldz, isuppz,
               work, &lwork, iwork, &liwork,
               &info);

        if (info == 0) {
            lwork  = static_cast<mwSignedIndex>(work[0]);
            liwork = iwork[0];
        }

        return info == 0;
    }

    bool solve()
    {
        auto JobZ  = C_.job_z_;
        auto Range = C_.range_;
        auto UpLo  = C_.uplo_;

        auto N    = A_.LD();

        auto* A   = A_.Data();
        auto  lda = A_.LD();

        auto vl = C_.vl_;
        auto vu = C_.vu_;

        mwSignedIndex il = 1, iu = N;

        double abstol;
        {
            char param[] = "Safe minimum";
            abstol = dlamch(param);
        }

        auto* w = E_.Lambda_.Data();
        auto* z = E_.Vectors_.Data();

        auto ldz = E_.Vectors_.LD();

        auto* isuppz = C_.ISuppZ_.Data();

        auto* work  = C_.Work_.Data();
        auto  lwork = C_.Work_.LD();

        auto* iwork  = C_.IWork_.Data();
        auto  liwork = C_.IWork_.LD();

        mwSignedIndex m, info = 0;

        dsyevr(&JobZ, &Range, &UpLo,
               &N, A, &lda, &vl, &vu, &il, &iu,
               &abstol, &m, w, z, &ldz, isuppz,
               work, &lwork, iwork, &liwork,
               &info);

        return info == 0;
    }
};


class BlockBoundaries
{
public:
    explicit BlockBoundaries(const mxArray* BSZ)
    {
        const auto nblk = mxGetNumberOfElements(BSZ);

        if (mxIsInt32(BSZ))       { construct<int>   (BSZ, nblk); }
        else if (mxIsDouble(BSZ)) { construct<double>(BSZ, nblk); }
        else {
            mexErrMsgTxt("SZ must be DOUBLE or INT32");
        }
    }

    using BlockID  = std::size_t;
    using SizeType = std::size_t;

    BlockID numBlocks() const { return p1_.size() - 1; }

    SizeType n1() const { return p1_.back(); }
    SizeType n2() const { return p2_.back(); }

    SizeType p1(const BlockID blk) const
    {
        mxAssert (blk < numBlocks(), "Internal Error");

        return p1_[blk];
    }

    SizeType p2(const BlockID blk) const
    {
        mxAssert (blk < numBlocks(), "Internal Error");

        return p2_[blk];
    }

    SizeType size(const BlockID blk) const
    {
        mxAssert (blk < numBlocks(), "Internal Error");

        return p1_[blk + 1] - p1_[blk + 0];
    }

private:
    using SizeVector = std::vector<SizeType>;

    SizeVector p1_;
    SizeVector p2_;

    template <typename BlockSizeType>
    void construct(const mxArray*    BSZ,
                   const std::size_t nblk)
    {
        p1_.reserve(nblk + 1);  p1_.push_back(0);
        p2_.reserve(nblk + 1);  p2_.push_back(0);

        auto bsz = static_cast<const BlockSizeType*>(mxGetData(BSZ));

        for (auto end = bsz + nblk; bsz != end; ++bsz) {
            const auto n = static_cast<SizeType>(*bsz);

            p1_.push_back(p1_.back() + n);
            p2_.push_back(p2_.back() + (n * n));
        }
    }
};


class MXDoubleVector
{
public:
    explicit MXDoubleVector(const mwSize n)
        : x_(mxCreateDoubleMatrix(n, 1, mxREAL))
    {}

    double* Data(const mwSize i = 0)
    {
        mxAssert (i < Size(), "Internal Error");

        return mxGetPr(Array()) + i;
    }

    mxArray* ReleaseMXArray()
    {
        return x_.release();
    }

    std::size_t Size() const
    {
        return mxGetNumberOfElements(Array());
    }

private:
    struct Delete {
        void operator()(mxArray* x)
        {
            if (x != nullptr) {
                mxDestroyArray(x);
            }
        }
    };

    using MXArray = std::unique_ptr<mxArray, Delete>;

    MXArray x_;

    mxArray* Array() const
    {
        return x_.get();
    }
};


class MEXResult
{
public:
    MEXResult(const int              nlhs,
              const BlockBoundaries& blocks)
        : blocks_(blocks)
    {
        mxAssert ((nlhs == 1) || (nlhs == 2),
                  "Must be one or two return values.");

        result_.emplace_back(blocks_.n1());

        if (nlhs > 1) {
            result_.emplace_back(blocks_.n2());
        }
    }

    double* EigenValues(const BlockBoundaries::BlockID blockID)
    {
        mxAssert (blockID < blocks_.numBlocks(), "Internal Error");

        return result_[Lambda].Data(blocks_.p1(blockID));
    }

    double* EigenVectors(const BlockBoundaries::BlockID blockID)
    {
        mxAssert (blockID < blocks_.numBlocks(), "Internal Error");

        if (result_.size() < (InvSubspace + 1)) {
            return nullptr;
        }

        return result_[InvSubspace].Data(blocks_.p2(blockID));
    }

    void ExtractResultArrays(mxArray* plhs[])
    {
        for (auto& x : result_) {
            *plhs++ = x.ReleaseMXArray();
        }
    }

private:
    enum { Lambda      = 0 ,
           InvSubspace = 1 };

    const BlockBoundaries&      blocks_;
    std::vector<MXDoubleVector> result_;
};


namespace {
    //  d     = multiSymmEig(A, sz)
    // [d, v] = multiSymmEig(A, sz)
    bool args_ok(const int nlhs, const int nrhs, const mxArray* prhs[])
    {
        auto ok = nrhs == 2;

        ok = ok && ((nlhs == 1) || (nlhs == 2));
        ok = ok && (!mxIsEmpty(prhs[0]) && mxIsDouble(prhs[0]));
        ok = ok && (!mxIsEmpty(prhs[1]) &&
                    (mxIsDouble(prhs[1]) || mxIsInt32(prhs[1])));

        return ok;
    }

    void solveEigenProblem(const BlockBoundaries::BlockID blockID,
                           const double* const            Ai,
                           const mwSignedIndex            n,
                           MEXResult&                     result)
    {
#define VERBOSE_PRINT 0

#if VERBOSE_PRINT
        mexPrintf("Subproblem %llu, Size %lld, Thread %d/%d\n",
                  static_cast<unsigned long long>(blockID),
                  static_cast<long long>(n),
                  omp_get_thread_num() + 1, omp_get_num_threads());
#endif  // VERBOSE_PRINT

        auto p = EVProblem();

        if (p.solve(n, Ai)) {
            p.eigenValues(result.EigenValues(blockID));

            if (auto v = result.EigenVectors(blockID)) {
                p.eigenVectors(v);
            }
        }
    }
} // Anonymous


//  d     = multiSymmEig(A, sz)
// [d, v] = multiSymmEig(A, sz)
void mexFunction(int nlhs, mxArray*       plhs[],
                 int nrhs, const mxArray* prhs[])
{
    if (args_ok(nlhs, nrhs, prhs)) {
        const auto blocks = BlockBoundaries(prhs[1]);
        auto       result = MEXResult(nlhs, blocks);

        const double* const A = mxGetPr(prhs[0]);

#pragma omp parallel                                            \
    if ((blocks.numBlocks() > (50 * omp_get_num_threads())))

#pragma omp single
        {
            for (decltype(blocks.numBlocks())
                     b = 0, nb = blocks.numBlocks();
                 b < nb; ++b)
            {
                const double* const Ai = A + blocks.p2(b);
                const auto          n  = blocks.size(b);

#pragma omp task
                solveEigenProblem(b, Ai, n, result);
            }
        }

        result.ExtractResultArrays(plhs);
    }
    else {
        const std::string fn = mexFunctionName();
        std::ostringstream msg;

        msg << "Syntax:\n\t"
            << " d     = " << fn << "(A, sz) % or\n\t"
            << "[d, v] = " << fn << "(A, sz)";

        mexErrMsgTxt(msg.str().c_str());
    }
}
