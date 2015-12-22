#include <mex.h>
#include <lapack.h>

#include <algorithm>
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

    void eigenvalues(double* pr, double* pi) const
    {
        const auto& R = E_.RealPart_.Data();
        std::copy(R.begin(), R.end(), pr);

        const auto& I = E_.ImagPart_.Data();
        std::copy(I.begin(), I.end(), pi);
    }

    void rightEigenVectors(double* rev) const
    {
        const auto& V = E_.RightVect_.Data();
        std::copy(V.begin(), V.end(), rev);
    }

    void leftEigenVectors(double* lev) const
    {
        const auto& W = E_.LeftVect_.Data();
        std::copy(W.begin(), W.end(), lev);
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
        Array<> RealPart_;
        Array<> ImagPart_;

        Array<> RightVect_;
        Array<> LeftVect_;

        void resize(const mwSignedIndex n)
        {
            RealPart_.resize(n);
            ImagPart_.resize(n);

            RightVect_.resize(n, n);
            LeftVect_ .resize(n, n);
        }
    };

    struct ProblemCharacteristics {
        const char balance_{ 'B' };
        const char job_vr_ { 'V' };
        const char job_vl_ { 'V' };
        const char sense_  { 'B' };

        // Results of scaling the matrix.
        Array<> Scale_;

        // Reciprocal condition estimates
        Array<> RCondE_;        // Eigenvalues
        Array<> RCondV_;        // Right eigenvectors

        // Work arrays
        Array<>              Work_;
        Array<mwSignedIndex> IWork_;

        void resize(const mwSignedIndex n,
                    const mwSignedIndex lwork)
        {
            Scale_.resize(n);

            RCondE_.resize(n);
            RCondV_.resize(n);

            const auto One = static_cast<mwSignedIndex>(1);

            Work_ .resize(lwork);
            IWork_.resize(2 * (std::max(n, One) - 1));
        }
    };

    // System matrix.
    Array<> A_;

    // Numerical output
    EigenResults E_;

    //
    ProblemCharacteristics C_;

    void resize(const mwSignedIndex n)
    {
        A_.resize(n, n);
        E_.resize(n);
        C_.resize(n, optimalWorkSize(n));
    }

    mwSignedIndex optimalWorkSize(const mwSignedIndex n) const
    {
        mwSignedIndex N = n;

        Array<> work;  work.resize(1);

        double* A  = nullptr;  mwSignedIndex lda = n;
        double* wr = nullptr;
        double* wi = nullptr;
        double* vl = nullptr;  mwSignedIndex ldvl = n;
        double* vr = nullptr;  mwSignedIndex ldvr = n;

        double* scale  = nullptr;
        double* rconde = nullptr;
        double* rcondv = nullptr;

        mwSignedIndex  lwork = -1; // lwork query.
        mwSignedIndex* iwork = nullptr;

        double        abnrm;
        mwSignedIndex ilo, ihi, info = 0;

        auto Balance = C_.balance_;
        auto JobVL   = C_.job_vl_;
        auto JobVR   = C_.job_vr_;
        auto Sense   = C_.sense_;

        dgeevx(&Balance, &JobVL, &JobVR, &Sense,
               &N, A, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
               &ilo, &ihi, scale, &abnrm, rconde, rcondv, work.Data(),
               &lwork, iwork, &info);

        return static_cast<mwSignedIndex>(work.Data()[0]);
    }

    bool solve()
    {
        auto N    = A_.LD();

        auto* A   = A_.Data();
        auto  lda = A_.LD();

        auto* wr = E_.RealPart_.Data();
        auto* wi = E_.ImagPart_.Data();

        auto* vl   = E_.LeftVect_.Data();
        auto  ldvl = E_.LeftVect_.LD();

        auto* vr   = E_.RightVect_.Data();
        auto  ldvr = E_.RightVect_.LD();

        auto* scale  = C_.Scale_ .Data();
        auto* rconde = C_.RCondE_.Data();
        auto* rcondv = C_.RCondV_.Data();

        auto* work  = C_.Work_ .Data();
        auto  lwork = C_.Work_ .LD();
        auto* iwork = C_.IWork_.Data();

        double        abnrm;
        mwSignedIndex ilo, ihi, info = 0;

        auto Balance = C_.balance_;
        auto JobVL   = C_.job_vl_;
        auto JobVR   = C_.job_vr_;
        auto Sense   = C_.sense_;

        dgeevx(&Balance, &JobVL, &JobVR, &Sense,
               &N, A, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
               &ilo, &ihi, scale, &abnrm, rconde, rcondv, work,
               &lwork, iwork, &info);

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
    //  d        = multiEig(A, sz)
    // [d, v]    = multiEig(A, sz)
    // [d, v, w] = multiEig(A, sz)
    bool args_ok(const int nlhs, const int nrhs, const mxArray* prhs[])
    {
        auto ok = nrhs == 2;

        ok = ok && ((nlhs >= 1) && (nlhs <= 3));
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


//  d        = multiEig(A, sz)
// [d, v]    = multiEig(A, sz)
// [d, v, w] = multiEig(A, sz)
void mexFunction(int nlhs, mxArray*       plhs[],
                 int nrhs, const mxArray* prhs[])
{
    if (args_ok(nlhs, nrhs, prhs)) {
        const auto bsz = blockSizes(prhs[1]);

        if (! sizesValid(bsz)) {
            mexErrMsgTxt("All block sizes must be strictly positive");
        }

        const auto rightEV = nlhs > 1;
        const auto leftEV  = nlhs > 2;

        const auto N = TotalEntries(bsz);

        plhs[0] = mxCreateDoubleMatrix(N.n_, 1, mxCOMPLEX);
        auto* pr = mxGetPr(plhs[0]);
        auto* pi = mxGetPi(plhs[0]);

        double* v = nullptr;
        if (rightEV) {
            plhs[1] = mxCreateDoubleMatrix(N.n2_, 1, mxREAL);
            v       = mxGetPr(plhs[1]);
        }

        double* w = nullptr;
        if (leftEV) {
            plhs[2] = mxCreateDoubleMatrix(N.n2_, 1, mxREAL);
            w       = mxGetPr(plhs[2]);
        }

        auto        ev = EVProblem(N.max_);
        const auto* A  = mxGetPr(prhs[0]);

        for (const auto& sz : bsz) {
            if (ev.solve(sz, A)) {
                ev.eigenvalues(pr, pi);

                if (rightEV) { ev.rightEigenVectors(v); }
                if (leftEV)  { ev.leftEigenVectors (w); }
            }

            const auto sz2 = sz * sz;

            A  += sz2;
            pr += sz;
            pi += sz;

            if (rightEV) { v += sz2; }
            if (leftEV)  { w += sz2; }
        }
    }
    else {
        const std::string fn = mexFunctionName();
        std::ostringstream msg;

        msg << "Syntax:\n\t"
            << " d        = " << fn << "(A, sz) % or\n\t"
            << "[d, v]    = " << fn << "(A, sz) % or\n\t"
            << "[d, v, w] = " << fn << "(A, sz)";

        mexErrMsgTxt(msg.str().c_str());
    }
}
