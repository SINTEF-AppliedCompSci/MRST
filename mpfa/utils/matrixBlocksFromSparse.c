#include <stddef.h>

#include <mex.h>

static int
args_ok(int nlhs, int nrhs, const mxArray *prhs[]);

static size_t *
extract_block_sizes(const mxArray *SZ);

static size_t
count_dense_block_entries(const size_t  nblk,
                          const size_t *rsz ,
                          const size_t *csz );

static mxArray *
get_blocks_from_sparse(const size_t   nblk,
                       const mwSize  *rsz ,
                       const mwSize  *csz ,
                       const mxArray *A   );

static mxArray *
compute_block_pointers(const size_t  nblk,
                       const size_t *rsz ,
                       const size_t *csz );

/*
 * Extract blocks, possibly non-square, on the diagonal of a sparse matrix.
 * This is a workaround since indexing of the form A(k) does not work in
 * MATLAB in the case of very large matrices.
 *
 * Syntax:
 *    s     = matrixBlocksFromSparse(S, rsz)
 *    s     = matrixBlocksFromSparse(S, rsz, csz)
 *   [s, p] = matrixBlocksFromSparse(...)
 */
/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    char errmsg[1023 + 1];

    size_t  nblk;
    size_t *rsz, *csz;

    if (args_ok(nlhs, nrhs, prhs)) {
        csz = rsz = extract_block_sizes(prhs[1]);

        if (nrhs == 3) {
            /* Rectangular matrix blocks:
             *
             *   s = matrixBlocksFromSparse(S, rsz, csz)
             */
            csz = extract_block_sizes(prhs[2]);
        }

        if ((rsz != NULL) && (csz != NULL)) {
            nblk    = mxGetNumberOfElements(prhs[1]);
            plhs[0] = get_blocks_from_sparse(nblk, rsz, csz, prhs[0]);

            if (nlhs == 2) {
                plhs[1] = compute_block_pointers(nblk, rsz, csz);
            }

            if (nrhs == 3) {
                mxFree(csz);
            }

            mxFree(rsz);
        }
        else {
            if ((nrhs == 3) && (csz != NULL)) { mxFree(csz); }
            if (rsz != NULL)                  { mxFree(rsz); }

            sprintf(errmsg,
                    "%s(): Failed to allocate sufficient memory "
                    "to hold copy of block size vector.",
                    mexFunctionName());

            mexErrMsgTxt(errmsg);
        }
    }
    else {
        sprintf(errmsg,
                "Calling sequence:\n\t"
                " Sbd       = %s(S, rsz)      %% or\n\t"
                " Sbd       = %s(S, rsz, csz) %% or\n\t"
                "[Sbd, pos] = %s(...)",
                mexFunctionName(),
                mexFunctionName(),
                mexFunctionName());

        mexErrMsgTxt(errmsg);
    }
}


/*
 *   s     = matrixBlocksFromSparse(S, rsz)
 *   s     = matrixBlocksFromSparse(S, rsz, csz)
 *  [s, p] = matrixBlocksFromSparse(...)
 */
/* ---------------------------------------------------------------------- */
static int
args_ok(int nlhs, int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok;

    ok = ((nlhs == 1) || (nlhs == 2)) &&
         ((nrhs == 2) || (nrhs == 3));

    ok = ok && !mxIsEmpty(prhs[0])
            && mxIsSparse(prhs[0]) && !mxIsComplex(prhs[0]);
    ok = ok && !mxIsEmpty(prhs[1])
            && (mxIsInt32(prhs[1]) || mxIsDouble(prhs[1]))
            && !mxIsComplex(prhs[1]);

    if (ok && (nrhs == 3)) {
        ok = ok && (mxIsInt32(prhs[2]) || mxIsDouble(prhs[2]))
                && (mxGetNumberOfElements(prhs[2]) ==
                    mxGetNumberOfElements(prhs[1]))
                && !mxIsComplex(prhs[2]);
    }

    return ok;
}


/* ---------------------------------------------------------------------- */
static size_t *
extract_block_sizes(const mxArray *SZ)
/* ---------------------------------------------------------------------- */
{
    size_t  nblk, b;

    int    *pi;
    double *pd;
    size_t *sz;

    nblk = mxGetNumberOfElements(SZ);
    sz   = mxMalloc(nblk * sizeof *sz);

    if (sz != NULL) {

        if (mxIsInt32(SZ)) {
            pi = mxGetData(SZ);

            for (b = 0; b < nblk; b++) {
                sz[ b ] = pi[ b ];
            }
        }
        else {
            mxAssert (mxIsDouble(SZ), "Internal error.");

            pd = mxGetPr(SZ);

            for (b = 0; b < nblk; b++) {
                sz[ b ] = pd[ b ];
            }
        }
    }

    return sz;
}


/* ---------------------------------------------------------------------- */
static size_t
count_dense_block_entries(const size_t  nblk,
                          const size_t *rsz ,
                          const size_t *csz )
/* ---------------------------------------------------------------------- */
{
    size_t b, sum_mxn;

    sum_mxn = 0;

    for (b = 0; b < nblk; b++) {
        sum_mxn += rsz[b] * csz[b];
    }

    return sum_mxn;
}


/* ---------------------------------------------------------------------- */
static mxArray *
get_blocks_from_sparse(const size_t   nblk,
                       const size_t  *rsz ,
                       const size_t  *csz ,
                       const mxArray *A   )
/* ---------------------------------------------------------------------- */
{
    size_t   b, m, n, I, J, sum_mxn;

    mwIndex  i, j, k, e;
    mwIndex *ir , *jc;
    double  *blk, *sa;

    mxArray *ret;

    sum_mxn = count_dense_block_entries(nblk, rsz, csz);

    /* ret = zeros([sum_mxn, 1]) */
    ret = mxCreateDoubleMatrix(sum_mxn, 1, mxREAL);

    ir  = mxGetIr( A );
    jc  = mxGetJc( A );
    sa  = mxGetPr( A );
    blk = mxGetPr(ret);

    I = J = k = 0;
    for (b = 0; b < nblk; b++) {
        m = rsz[ b ];
        n = csz[ b ];

        for (j = 0; j < n; j++) {
            for (k = jc[J + j], e = jc[J + j + 1]; k < e; k++) {
                i = ir[k];

                if ((i < I) || (i >= I + m)) {
                    /* Entry outside block.  Ignore. */
                    continue;
                }

                i -= I;

                blk[j*m + i] = sa[ k ];
            }
        }

        blk += m * n;
        I   += m    ;
        J   +=     n;
    }

    return ret;
}


/* ---------------------------------------------------------------------- */
static mxArray *
compute_block_pointers(const size_t  nblk,
                       const size_t *rsz ,
                       const size_t *csz )
/* ---------------------------------------------------------------------- */
{
    size_t   b;
    double  *p;
    mxArray *pos;

    /* pos = zeros([nblk + 1, 1]) */
    pos = mxCreateDoubleMatrix(nblk + 1, 1, mxREAL);
    p   = mxGetPr(pos);

    p[0] = 1; /* One-based indexing in MATLAB */

    for (b = 0; b < nblk; b++) {
        p[b + 1] = p[b] + (rsz[b] * csz[b]);
    }

    return pos;
}
