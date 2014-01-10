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

    return ok && (nrhs == 2);
}


/* ---------------------------------------------------------------------- */
static void
deallocate_aux_arrays(int *p, int *neigh)
/* ---------------------------------------------------------------------- */
{
    if (p     != NULL) { mxFree(p)    ; }
    if (neigh != NULL) { mxFree(neigh); }
}


/* ---------------------------------------------------------------------- */
static int
allocate_aux_arrays(int nc, int nneigh, int **p, int **neigh)
/* ---------------------------------------------------------------------- */
{
    int ret, *t, *n;

    t = mxMalloc(nc         * sizeof *t);
    n = mxMalloc(2 * nneigh * sizeof *n);

    if ((t == NULL) || (n == NULL)) {
        deallocate_aux_arrays(t, n);

        *p     = NULL;
        *neigh = NULL;

        ret = -1;
    } else {
        *p     = t;
        *neigh = n;

        ret = nc;
    }

    return ret;
}


/* ---------------------------------------------------------------------- */
static void
extract_neighbour(const mxArray *M_neigh, int *neigh)
/* ---------------------------------------------------------------------- */
{
    size_t e, ne;

    int    *pi;
    double *pd;

    ne = mxGetM(M_neigh);

    if (mxIsDouble(M_neigh)) {
        pd = mxGetPr(M_neigh);

        for (e = 0; e < ne; e++) {
            neigh[2*e + 0] = pd[e + 0*ne] - 1;
            neigh[2*e + 1] = pd[e + 1*ne] - 1;
        }
    } else {
        pi = mxGetData(M_neigh);

        for (e = 0; e < ne; e++) {
            neigh[2*e + 0] = pi[e + 0*ne] - 1;
            neigh[2*e + 1] = pi[e + 1*ne] - 1;
        }
    }
}


/* ---------------------------------------------------------------------- */
static void
extract_int_vec(const mxArray *M_v, int *v)
/* ---------------------------------------------------------------------- */
{
    size_t e, ne;

    int    *pi;
    double *pd;

    ne = mxGetNumberOfElements(M_v);

    if (mxIsDouble(M_v)) {
        pd = mxGetPr(M_v);

        for (e = 0; e < ne; e++) { v[e] = pd[e] - 1; }
    } else {
        pi = mxGetData(M_v);

        for (e = 0; e < ne; e++) { v[e] = pi[e] - 1; }
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
 * p = mex_compress_partition(p, N)
 */

/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int  nc, nneigh, nextra, *p, *neigh;
    char errmsg[1023 + 1];

    if (args_ok(nlhs, nrhs, prhs)) {
        plhs[0] = mxDuplicateArray(prhs[0]);
        nc      = mxGetNumberOfElements(prhs[0]);
        nneigh  = mxGetM(prhs[1]);

        mxAssert (mxGetN(prhs[1]) == 2,
                  "Neighbourship must be M-by-2");

        if (allocate_aux_arrays(nc, nneigh, &p, &neigh) == nc) {
            extract_int_vec(prhs[0], p);
            extract_neighbour(prhs[1], neigh);

            nextra = partition_split_disconnected(nc, nneigh, neigh, p);

            if (nextra == 0) {
                mexPrintf("Partition processing complete. No extra blocks.\n");

                assign_lhs(p, plhs[0]);
            } else if (nextra > 0) {
                mexPrintf("Partition processing complete, %d extra block%s.\n",
                          nextra, (nextra > 1) ? "s" : "");

                assign_lhs(p, plhs[0]);
            } else {
                mexPrintf("Partition processing failed. Insufficient memory?\n");
            }

            deallocate_aux_arrays(p, neigh);
        }
    } else {
        sprintf(errmsg,
                "Calling sequence is\n"
                "\tp2 = %s(p1, Neighbourship)\n",
                mexFunctionName());

        mexErrMsgTxt(errmsg);
    }
}
