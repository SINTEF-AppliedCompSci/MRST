#include <stddef.h>

#include <mex.h>

static int
args_ok(int nlhs, int nrhs, const mxArray *prhs[]);

static void
extract_block_sizes(const mxArray *SZ     ,
                    size_t        *sum_sz2,
                    size_t        **sz    );
static void
get_blocks_from_sparse(size_t         nblk,
                       const mwSize  *sz  ,
                       const mxArray *A   ,
                       mxArray       *a   );

/*
 * Extract blocks on the diagonal of sparse matrix.  This is a workaround
 * since indexing of the form A(k) does not work in Matlab for very large
 * matrices.
 *
 * Syntax:
 *   s = matrixBlocksFromSparse(S, sz)
 */
/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    char errmsg[1023 + 1];

    size_t  nblk, sum_sz2;
    size_t *sz;

    if (args_ok(nlhs, nrhs, prhs)) {
        extract_block_sizes(prhs[1], &sum_sz2, &sz);

        if (sz != NULL) {
            plhs[0] = mxCreateDoubleMatrix(sum_sz2, 1, mxREAL);
            nblk    = mxGetNumberOfElements(prhs[1]);

            get_blocks_from_sparse(nblk, sz, prhs[0], plhs[0]);

            mxFree(sz);
        }
        else {
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
                "s = %s(S, sz)", mexFunctionName());
        mexErrMsgTxt(errmsg);
    }
}


/* ---------------------------------------------------------------------- */
static int
args_ok(int nlhs, int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok;

    ok = (nlhs == 1) && (nrhs == 2);

    ok = ok && !mxIsEmpty(prhs[0])
            && mxIsSparse(prhs[0]) && !mxIsComplex(prhs[0]);
    ok = ok && !mxIsEmpty(prhs[1])
            && (mxIsInt32(prhs[1]) || mxIsDouble(prhs[1]))
            && !mxIsComplex(prhs[1]);

    ok = ok && (mxGetM(prhs[0]) == mxGetN(prhs[0]));

    return ok;
}


/* ---------------------------------------------------------------------- */
static void
extract_block_sizes(const mxArray *SZ     ,
                    size_t        *sum_sz2,
                    size_t        **sz    )
/* ---------------------------------------------------------------------- */
{
    size_t  nblk, b;

    int    *pi;
    double *pd;
    size_t *sz0;

    nblk = mxGetNumberOfElements(SZ);
    sz0  = mxMalloc(nblk * sizeof *sz0);

    if (sz0 != NULL) {

        if (mxIsInt32(SZ)) {
            pi = mxGetData(SZ);

            *sum_sz2 = 0;
            for (b = 0; b < nblk; b++) {
                sz0[ b ]  = pi[ b ];
                *sum_sz2 += sz0[b] * sz0[b];
            }
        }
        else {
            mxAssert (mxIsDouble(SZ), "Internal error.");

            pd = mxGetPr(SZ);

            *sum_sz2 = 0;
            for (b = 0; b < nblk; b++) {
                sz0[ b ]  = pd[ b ];
                *sum_sz2 += sz0[b] * sz0[b];
            }
        }

        *sz = sz0;
    }
    else {
        *sz = NULL;
    }
}


/* ---------------------------------------------------------------------- */
static void
get_blocks_from_sparse(size_t         nblk,
                       const mwSize  *sz  ,
                       const mxArray *A   ,
                       mxArray       *a   )
/* ---------------------------------------------------------------------- */
{
    size_t   b, bsz, start;

    mwIndex  i, j, k;
    mwIndex *ir, *jc;
    double  *blk, *sa;

    ir  = mxGetIr(A);
    jc  = mxGetJc(A);
    sa  = mxGetPr(A);
    blk = mxGetPr(a);

    start = k = 0;
    for (b = 0; b < nblk; b++)
    {
        bsz = sz[ b ];

        for (j = 0; j < bsz; j++) {
            for (k = jc[start + j]; k < jc[start + j + 1]; k++) {

                i = ir[ k ] - start;

                blk[j*bsz + i] = sa[ k ];
            }
        }

        blk   += bsz * bsz;
        start += bsz;
    }
}
