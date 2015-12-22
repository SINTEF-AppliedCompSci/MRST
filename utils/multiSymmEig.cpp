#include <mex.h>
#include <lapack.h>

#include <algorithm>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

class EVProblem
{
public:
    EVProblem() {}

    explicit EVProblem(const mwSignedIndex n)
    {
        resize(n);
    }

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


struct TotalEntries {
public:
    template <class BlockSizes>
    explicit TotalEntries(const BlockSizes& bsz)
        : max_(0)
        , n_  (0)
        , n2_ (0)
    {
        for (const auto& sz : bsz) {
            const auto n = static_cast<mwSignedIndex>(sz);

            max_  = std::max(max_, n);
            n_   += n;
            n2_  += n * n;
        }
    }

    mwSignedIndex max_;
    mwSignedIndex n_;
    mwSignedIndex n2_;
};


namespace {
    //  d        = multiSymmEig(A, sz)
    // [d, v]    = multiSymmEig(A, sz)
    // [d, v, w] = multiSymmEig(A, sz)
    bool args_ok(const int nlhs, const int nrhs, const mxArray* prhs[])
    {
        auto ok = nrhs == 2;

        ok = ok && ((nlhs == 1) || (nlhs == 2));
        ok = ok && (!mxIsEmpty(prhs[0]) && mxIsDouble(prhs[0]));
        ok = ok && (!mxIsEmpty(prhs[1]) &&
                    (mxIsDouble(prhs[1]) || mxIsInt32(prhs[1])));

        return ok;
    }

    std::vector<mwSignedIndex> blockSizes(const mxArray* SZ)
    {
        const auto n = mxGetNumberOfElements(SZ);

        std::vector<mwSignedIndex> sz;
        sz.reserve(n);

        if (mxIsInt32(SZ)) {
            const auto* bsz = static_cast<const int*>(mxGetData(SZ));

            std::transform(bsz, bsz + n, std::back_inserter(sz),
                           [](const int m)
                           {
                               return static_cast<mwSignedIndex>(m);
                           });
        }
        else {
            mxAssert (mxIsDouble(SZ), "Internal Error.");

            const auto* bsz = mxGetPr(SZ);

            std::transform(bsz, bsz + n, std::back_inserter(sz),
                           [](const double m)
                           {
                               return static_cast<mwSignedIndex>(m);
                           });
        }

        return sz;
    }

    template <class BlockSizes>
    bool sizesValid(const BlockSizes& bsz)
    {
        using BlockSize = typename BlockSizes::value_type;

        // All block sizes must be strictly positive.
        return std::accumulate(bsz.begin(), bsz.end(), true,
                               [](const bool b, const BlockSize& sz)
                               {
                                   return b && (sz > 0);
                               });
    }
} // Anonymous


//  d     = multiSymmEig(A, sz)
// [d, v] = multiSymmEig(A, sz)
void mexFunction(int nlhs, mxArray*       plhs[],
                 int nrhs, const mxArray* prhs[])
{
    if (args_ok(nlhs, nrhs, prhs)) {
        const auto bsz = blockSizes(prhs[1]);

        if (! sizesValid(bsz)) {
            mexErrMsgTxt("All block sizes must be strictly positive");
        }

        const auto calcEVect = nlhs > 1;
        const auto N         = TotalEntries(bsz);

        // Note: mxREAL because we assume symmetric matrices
        plhs[0] = mxCreateDoubleMatrix(N.n_, 1, mxREAL);
        auto* pr = mxGetPr(plhs[0]);

        double* v = nullptr;
        if (calcEVect) {
            plhs[1] = mxCreateDoubleMatrix(N.n2_, 1, mxREAL);
            v       = mxGetPr(plhs[1]);
        }

        auto        ev = EVProblem(N.max_);
        const auto* A  = mxGetPr(prhs[0]);

        for (const auto& sz : bsz) {
            if (ev.solve(sz, A)) {
                ev.eigenValues(pr);

                if (calcEVect) { ev.eigenVectors(v); }
            }

            const auto sz2 = sz * sz;

            A  += sz2;
            pr += sz;

            if (calcEVect) { v += sz2; }
        }
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
