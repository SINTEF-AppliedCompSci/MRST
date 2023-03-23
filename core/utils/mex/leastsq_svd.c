#include <stddef.h>
#include <string.h>

#include <mex.h>
#include <lapack.h>

#undef MAX
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

#if !defined(_WIN32)

#if !defined(dgelss)
#define dgelss dgelss_
#endif  /* !defined(dgelss) */

#endif  /* !defined(_WIN32) */

/* ======================================================================
 * Types
 * ====================================================================== */

/* Collection of arguments for DGELSS() */
struct LLSProblem {
    mwSignedIndex  M;
    mwSignedIndex  N;
    mwSignedIndex  NRHS;
    mwSignedIndex  LDA;
    mwSignedIndex  LDB;
    mwSignedIndex  LWORK;
    mwSignedIndex  RANK;
    mwSignedIndex  INFO;

    mwSignedIndex  work_cpty;

    double         RCOND;

    double        *A;
    double        *B;
    double        *S;
    double        *WORK;
};


/* ---------------------------------------------------------------------- */
static void
destroy_lls_structure(struct LLSProblem *lls)
/* ---------------------------------------------------------------------- */
{
    if (lls != NULL) {
        if (lls->A != NULL) { mxFree(lls->A); }

        mxFree(lls);
    }
}


/* ---------------------------------------------------------------------- */
static struct LLSProblem *
create_lls_structure(void)
/* ---------------------------------------------------------------------- */
{
    struct LLSProblem *lls, lls0 = { 0 };

    lls = mxMalloc(1 * sizeof *lls);
    if (lls != NULL) {
        *lls = lls0;

        /* Note: At least one WORK element needed to query DGELLS() for its
         * optimal LWORK value. */
        lls->A = mxMalloc(1 * sizeof *lls->A);

        if (lls->A == NULL) {
            destroy_lls_structure(lls);

            lls = NULL;
        }
        else {
            lls->work_cpty = 1;

            lls->B    = lls->A;
            lls->S    = lls->A;
            lls->WORK = lls->A;
        }
    }

    return lls;
}


/* ---------------------------------------------------------------------- */
static int
prepare_lls_problem(const mwSignedIndex  M,
                    const mwSignedIndex  N,
                    struct LLSProblem   *lls)
/* ---------------------------------------------------------------------- */
{
    mwSignedIndex opt_lwork;
    mwSignedIndex ABS_size;

    int     ret;
    void   *p;
    double *A;

    mxAssert (lls != NULL, "Internal Error");

    lls->M     = M;
    lls->N     = N;
    lls->NRHS  = 1;
    lls->LDA   = M;
    lls->LDB   = MAX(M, N);
    lls->LWORK = -1;            /* Query optimal LWORK */
    lls->RCOND = 1.0e-12;

    lls->INFO = 0;

    dgelss(&lls->M, &lls->N, &lls->NRHS,
           lls->A, &lls->LDA,
           lls->B, &lls->LDB,
           lls->S, &lls->RCOND, &lls->RANK,
           lls->WORK, &lls->LWORK, &lls->INFO);

    /* Total size requirement for ->A, ->B, and ->S collectively. */
    ABS_size = M*N + M + N;

    ret = 1;

    if ((opt_lwork = lls->WORK[0]) > MAX(1, lls->work_cpty - ABS_size))
    {
        p = mxRealloc(lls->A, (ABS_size + opt_lwork) * sizeof *lls->A);

        if (p != NULL) {
            lls->A         = p;
            lls->work_cpty = ABS_size + opt_lwork;
        }
        else {
            ret = 0;
        }
    }

    if (ret != 0) {
        lls->LWORK = lls->work_cpty - ABS_size;

        lls->B    = lls->A + (M * N);
        lls->S    = lls->B + lls->LDB;
        lls->WORK = lls->B + (M + N);
    }

    return ret;
}


/* ---------------------------------------------------------------------- */
static mwSignedIndex
solve_lls_problem(const double      *A,
                  const double      *b,
                  double            *x,
                  struct LLSProblem *lls)
/* ---------------------------------------------------------------------- */
{
    memcpy(lls->A, A, lls->M * lls->N * sizeof *lls->A);
    memcpy(lls->B, b, lls->M          * sizeof *lls->B);

    lls->INFO = 0;

    dgelss(&lls->M, &lls->N, &lls->NRHS,
           lls->A, &lls->LDA,
           lls->B, &lls->LDB,
           lls->S, &lls->RCOND, &lls->RANK,
           lls->WORK, &lls->LWORK, &lls->INFO);

    if (lls->INFO == 0) {
        memcpy(x, lls->B, lls->N * sizeof *x);
    }

    return lls->INFO;
}


/*
 * Syntax:
 *   x = leastsq_svd(A, b, rsz, csz)
 */
/* ---------------------------------------------------------------------- */
static int
args_ok(int nlhs, int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok;

    ok = (nrhs == 4) && (nlhs == 1);

    ok = ok &&   mxIsDouble(prhs[0]) ;
    ok = ok &&   mxIsDouble(prhs[1]) ;
    ok = ok && ( mxIsDouble(prhs[2]) || mxIsInt32(prhs[2]) );
    ok = ok && ( mxIsDouble(prhs[3]) || mxIsInt32(prhs[3]) );
    ok = ok && (mxGetNumberOfElements(prhs[2]) ==
                mxGetNumberOfElements(prhs[3]));

    return ok;
}


/* ---------------------------------------------------------------------- */
static mwSize *
extract_block_sizes(const mxArray *SZ)
/* ---------------------------------------------------------------------- */
{
    size_t  i, n;
    int    *pi;
    double *pd;

    mwSize *sz;

    n = mxGetNumberOfElements(SZ);

    sz = mxMalloc(n * sizeof *sz);

    if (sz != NULL) {
        if (mxIsDouble(SZ)) {
            pd = mxGetPr(SZ);

            for (i = 0; i < n; i++) { sz[ i ] = pd[ i ]; }
        }
        else {
            mxAssert (mxIsInt32(SZ), "Internal error.");

            pi = mxGetData(SZ);

            for (i = 0; i < n; i++) { sz[ i ] = pi[ i ]; }
        }
    }

    return sz;
}


/* ---------------------------------------------------------------------- */
static size_t
accumulate_size(const size_t nblk, const mwSize *sz)
/* ---------------------------------------------------------------------- */
{
    size_t blk, ret;

    ret = 0;
    for (blk = 0; blk < nblk; blk++) {
        ret += sz[blk];
    }

    return ret;
}


/* ---------------------------------------------------------------------- */
static mxArray *
solve_all_lls_problems(const double *A,
                       const double *b,
                       const size_t  n,
                       const mwSize *rsz,
                       const mwSize *csz)
/* ---------------------------------------------------------------------- */
{
    size_t             i;
    mwSize             M, N;
    struct LLSProblem *lls;

    double  *x;
    mxArray *X;

    lls = create_lls_structure();
    X   = mxCreateDoubleMatrix(accumulate_size(n, csz), 1, mxREAL);

    if (lls != NULL) {
        x = mxGetPr(X);

        for (i = 0; i < n; i++) {
            M = rsz[i];
            N = csz[i];

            if (prepare_lls_problem(M, N, lls) != 0) {
                /* Setup complete.  Affect solution. */
                solve_lls_problem(A, b, x, lls);
            }

            A += M * N;
            b += M    ;
            x +=     N;
        }
    }

    destroy_lls_structure(lls);

    return X;
}


/*
 * Solve a sequence of many, typically small, linear least squares (LLS)
 * problems using the LAPACK routine DGELSS (based on SVD factorisation of
 * the coefficient matrices, supports rank-deficient problems).  The matrix
 * entries, right-hand side and dimensions of each matrix is input, and the
 * entries of the system solutions are output.
 *
 * Syntax:
 *    x = leastsq_svd(A, b, rsz, csz)
 */
/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    char    errmsg[1023 + 1];
    size_t  nblocks;
    mwSize *rsz, *csz;

    enum { A = 0, b = 1, RSZ = 2, CSZ = 3 };

    if (args_ok(nlhs, nrhs, prhs)) {
        rsz = extract_block_sizes(prhs[ RSZ ]);
        csz = extract_block_sizes(prhs[ CSZ ]);

        if ((rsz != NULL) && (csz != NULL)) {
            /* Invert all matrices. */
            nblocks = mxGetNumberOfElements(prhs[ RSZ ]);

            plhs[0] = solve_all_lls_problems(mxGetPr(prhs[ A ]),
                                             mxGetPr(prhs[ b ]),
                                             nblocks, rsz, csz);

            mxFree(csz);  mxFree(rsz);
        }

        else {
            if (csz != NULL) { mxFree(csz); }
            if (rsz != NULL) { mxFree(rsz); }

            sprintf(errmsg,
                    "%s(): Failed to allocate sufficient memory to hold "
                    "matrix block sizes.", mexFunctionName());
            mexErrMsgTxt(errmsg);
        }
    }

    else {
        sprintf(errmsg,
                "Syntax is\n\t"
                "x = %s(A, b, rsz, csz);", mexFunctionName());
        mexErrMsgTxt(errmsg);
    }
}
