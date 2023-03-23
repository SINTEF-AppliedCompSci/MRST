#include <stddef.h>
#include <stdio.h>

#include <mex.h>

#include <lapack.h>

static int
args_ok(int nlhs, int nrhs, const mxArray *prhs[]);

static void
invert_dense_matrices(size_t nblocks, const mwSize *sz, mxArray *A);

static mwSize *
extract_block_sizes(const mxArray *sz);

/*
 * Invert a sequence of many, typically small, square matrices using
 * LAPACK routines DGETRF and DGETRI.  The matrix entries and
 * dimensions of each matrix is input, and the matrix entries of the
 * inverses are output.
 *
 * Syntax:
 *   B = invv(A, sz)
 */
/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    char    errmsg[1023 + 1];
    size_t  nblocks;
    mwSize *sz;

    if (args_ok(nlhs, nrhs, prhs)) {
        sz = extract_block_sizes(prhs[1]);

        if (sz != NULL) {
            /* Create result */
            plhs[0] = mxDuplicateArray(prhs[0]);

            /* Invert all matrices. */
            nblocks = mxGetNumberOfElements(prhs[1]);
            invert_dense_matrices(nblocks, sz, plhs[0]);

            mxFree(sz);
        }

        else {
            sprintf(errmsg,
                    "%s(): Failed to allocate sufficient memory to hold "
                    "matrix block sizes.", mexFunctionName());
            mexErrMsgTxt(errmsg);
        }
    }

    else {
        sprintf(errmsg,
                "Syntax is\n\t"
                "B = %s(A, sz);", mexFunctionName());
        mexErrMsgTxt(errmsg);
    }
}


/* ---------------------------------------------------------------------- */
static void
invert_dense_matrices(size_t nblocks, const mwSize *sz, mxArray *A)
/* ---------------------------------------------------------------------- */
{
    int    issue_warning;
    size_t block, j, offset;
    mwSize maxsize;

    mwSize M, N, LD, INFO, LWORK;

    mwSize *IPIV;
    double *WORK, *a;

    a = mxGetPr(A);

    /* Find maximal system size */
    maxsize = 0;
    for (block = 0; block < nblocks; block++) {
        if (sz[ block ] > maxsize) {
            maxsize = sz[ block ];
        }
    }

    LWORK = maxsize;
    IPIV  = mxMalloc(LWORK * sizeof *IPIV);
    WORK  = mxMalloc(LWORK * sizeof *WORK);

    if ((IPIV != NULL) && (WORK != NULL)) {

        issue_warning = 0;
        offset        = 0;

        for (block = 0; block < nblocks; block++) {
            M = N = LD = sz[ block ];

            /*  Factor matrix */
            dgetrf(&M, &N, a + offset, &LD, IPIV, &INFO);

            if (INFO != 0) {
                size_t N2 = N * N;

                for (j = 0; j < N2; j++) {
                    a[offset + j] = mxGetInf();
                }

                issue_warning = 1;
            }

            else {
                /* Invert matrix */
                dgetri(&N, a + offset, &LD, IPIV, WORK, &LWORK, &INFO);

                if (INFO != 0) {
                    size_t N2 = N * N;

                    for (j = 0; j < N2; j) {
                        a[offset + j] = mxGetInf();
                    }

                    issue_warning = 1;
                }
            }

            offset += N * N;
        }

        if (issue_warning) {
            mexWarnMsgIdAndTxt("InvertSparseBlockMatrix:Singular",
                               "Matrix block is singular "
                               "to working precision.");
        }

    }

    if (WORK != NULL) { mxFree(WORK); }
    if (IPIV != NULL) { mxFree(IPIV); }
}


/* ---------------------------------------------------------------------- */
static int
args_ok(int nlhs, int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok;

    ok = (nrhs == 2) && (nlhs <= 1);

    ok = ok &&   mxIsDouble (prhs[0]) ;
    ok = ok && (!mxIsComplex(prhs[1]));
    ok = ok && ( mxIsDouble (prhs[1]) || mxIsInt32(prhs[1]) );

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
