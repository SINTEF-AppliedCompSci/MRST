#include <stddef.h>
#include <string.h>

#include <mex.h>

#include "partition.h"


/* ---------------------------------------------------------------------- */
static int
args_ok(int nlhs, int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int i, ok;

    ok = nlhs == 1;

    for (i = 0; ok && (i < nrhs); i++) {
        ok = mxIsDouble(prhs[i]) || mxIsInt32(prhs[i]);
    }

    return ok && (nrhs == 3);
}


/* ---------------------------------------------------------------------- */
static void
release_rhs(int *idx, int *fine_d, int *coarse_d)
/* ---------------------------------------------------------------------- */
{
    if (coarse_d != NULL) { mxFree(coarse_d); }
    if (fine_d   != NULL) { mxFree(fine_d);   }
    if (idx      != NULL) { mxFree(idx);      }
}


/* ---------------------------------------------------------------------- */
static void
copy_intvec(const mxArray *M_i, int *i)
/* ---------------------------------------------------------------------- */
{
    size_t e, ne;

    int    *pi;
    double *pd;

    ne = mxGetNumberOfElements(M_i);

    if (mxIsDouble(M_i)) {
        pd = mxGetPr(M_i);

        for (e = 0; e < ne; e++) { i[e] = pd[e]; }
    } else {
        pi = mxGetData(M_i);
        memcpy(i, pi, ne * sizeof *i);
    }
}


/* ---------------------------------------------------------------------- */
static void
extract_rhs(const mxArray *prhs[],
            int *nidx, int *ndims,
            int **idx, int **fine_d, int **coarse_d)
/* ---------------------------------------------------------------------- */
{
    int ni, nd, p, *i, *f, *c;

    ni = mxGetNumberOfElements(prhs[0]);
    nd = mxGetNumberOfElements(prhs[1]);

    if (mxGetNumberOfElements(prhs[2]) != (size_t) nd) {
        mexErrMsgTxt("'coarseDim' must have same number of "
                     "elements as 'fineDim'.");
    }

    i = mxMalloc(ni * sizeof *i);
    f = mxMalloc(nd * sizeof *f);
    c = mxMalloc(nd * sizeof *c);

    if ((i == NULL) || (f == NULL) || (c == NULL)) {
        release_rhs(i, f, c);

        *nidx = -1;    *ndims  = -1;
        *idx  = NULL;  *fine_d = NULL;  *coarse_d = NULL;
    } else {                      /* Adjust for 1-based indexing */
        copy_intvec(prhs[0], i);  for (p = 0; p < ni; p++) { i[p]--; }
        copy_intvec(prhs[1], f);
        copy_intvec(prhs[2], c);

        *nidx = ni;  *ndims  = nd;
        *idx  = i;   *fine_d = f;  *coarse_d = c;
    }
}


/* ---------------------------------------------------------------------- */
static void
assign_lhs(const int *p, mxArray *lhs)
/* ---------------------------------------------------------------------- */
{
    size_t e, ne;

    int    *pi;
    double *pd;

    ne = mxGetNumberOfElements(lhs);

    if (mxIsDouble(lhs)) {
        pd = mxGetPr(lhs);

        for (e = 0; e < ne; e++) { pd[e] = p[e] + 1; }
    } else {
        pi = mxGetData(lhs);

        for (e = 0; e < ne; e++) { pi[e] = p[e] + 1; }
    }
}


/*
 * p = mex_partition_ui(i, fineDim, coarseDim)
 */

/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int nidx, ndims, *idx, *p, *fine_d, *coarse_d;
    char errmsg[1023 + 1];

    if (args_ok(nlhs, nrhs, prhs)) {
        extract_rhs(prhs, &nidx, &ndims, &idx, &fine_d, &coarse_d);

        if (nidx > 0) {
            plhs[0] = mxDuplicateArray(prhs[0]);
            p       = mxMalloc(nidx * sizeof *p);

            partition_unif_idx(ndims, nidx, fine_d, coarse_d, idx, p);

            assign_lhs(p, plhs[0]);

            mxFree(p);    release_rhs(idx, fine_d, coarse_d);
        }
    } else {
        sprintf(errmsg,
                "Calling sequence is\n"
                "\tp = %s(ix, fineDim, coarseDim)\n",
                mexFunctionName());

        mexErrMsgTxt(errmsg);
    }
}
