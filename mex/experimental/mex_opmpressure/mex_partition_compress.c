#include <stddef.h>

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

    return ok && (nrhs == 1);
}


/* ---------------------------------------------------------------------- */
static void
free_aux_array(int *p)
/* ---------------------------------------------------------------------- */
{
    if (p != NULL) { mxFree(p); }
}


/* ---------------------------------------------------------------------- */
static int
allocate_aux_array(int n, int **p)
/* ---------------------------------------------------------------------- */
{
    int ret, *t;

    t = mxMalloc(n * sizeof *t);

    if (t == NULL) {
        ret = -1;
    } else {
        ret = n;
    }

    *p = t;

    return ret;
}


/* ---------------------------------------------------------------------- */
static void
extract_rhs(const mxArray *M_p, int *p)
/* ---------------------------------------------------------------------- */
{
    size_t e, ne;

    int    *pi;
    double *pd;

    ne = mxGetNumberOfElements(M_p);

    if (mxIsDouble(M_p)) {
        pd = mxGetPr(M_p);

        for (e = 0; e < ne; e++) { p[e] = pd[e] - 1; }
    } else {
        pi = mxGetData(M_p);

        for (e = 0; e < ne; e++) { p[e] = pi[e] - 1; }
    }
}


/* ---------------------------------------------------------------------- */
static void
assign_lhs(const int *p, mxArray *M_p)
/* ---------------------------------------------------------------------- */
{
    size_t e, ne;

    int    *pi;
    double *pd;

    ne = mxGetNumberOfElements(M_p);

    if (mxIsDouble(M_p)) {
        pd = mxGetPr(M_p);

        for (e = 0; e < ne; e++) { pd[e] = p[e] + 1; }
    } else {
        pi = mxGetData(M_p);

        for (e = 0; e < ne; e++) { pi[e] = p[e] + 1; }
    }
}


/*
 * p = mex_compress_partition(p)
 */

/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int  n, *p;
    char errmsg[1023 + 1];

    if (args_ok(nlhs, nrhs, prhs)) {
        plhs[0] = mxDuplicateArray(prhs[0]);
        n = mxGetNumberOfElements(prhs[0]);

        if (allocate_aux_array(n, &p) == n) {
            extract_rhs(prhs[0], p);

            partition_compress(n, p);

            assign_lhs(p, plhs[0]);

            free_aux_array(p);
        }
    } else {
        sprintf(errmsg,
                "Calling sequence is\n"
                "\tp2 = %s(p1)\n",
                mexFunctionName());

        mexErrMsgTxt(errmsg);
    }
}
