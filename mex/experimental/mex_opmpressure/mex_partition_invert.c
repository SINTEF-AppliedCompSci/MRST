#include <stddef.h>
#include <string.h>

#include <mex.h>

#include "partition.h"


#define MAX(a,b) (((a) > (b)) ? (a) : (b))


/* ---------------------------------------------------------------------- */
static int
args_ok(int nlhs, int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok;

    ok =       (nlhs == 2) || (nlhs == 3);
    ok = ok && (nrhs == 1);
    ok = ok && (mxIsDouble(prhs[0]) || mxIsInt32(prhs[0]));

    return ok;
}


/* ---------------------------------------------------------------------- */
static void
extract_rhs(const mxArray *M_p, int *p)
/* ---------------------------------------------------------------------- */
{
    size_t  e, ne;

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
assign_int_vec(const int *v, mxArray *M_v)
/* ---------------------------------------------------------------------- */
{
    size_t e, ne;

    int    *pi;
    double *pd;

    ne = mxGetNumberOfElements(M_v);

    if (mxIsDouble(M_v)) {
        pd = mxGetPr(M_v);

        for (e = 0; e < ne; e++) { pd[e] = v[e] + 1; }
    } else {
        pi = mxGetData(M_v);

        for (e = 0; e < ne; e++) { pi[e] = v[e] + 1; }
    }
}


/* ---------------------------------------------------------------------- */
static int
max_block(int nc, const int *p)
/* ---------------------------------------------------------------------- */
{
    int c, m;

    m = -1;
    for (c = 0; c < nc; c++)
    {
        m = MAX(m, p[c]);
    }

    return m;
}


/* ---------------------------------------------------------------------- */
static void
adjust_one_based_idx(size_t n, int *v)
/* ---------------------------------------------------------------------- */
{
    size_t i;

    for (i = 0; i < n; i++) { v[i] += 1; }
}


/*
 * [pb2c, b2c]      = mex_partition_invert(p)
 * [pb2c, b2c, loc] = mex_partition_invert(p)
 */

/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int  nc, max_blk, *p, *pb2c, *b2c, *loc;
    char errmsg[1023 + 1];

    if (args_ok(nlhs, nrhs, prhs)) {
        nc = mxGetNumberOfElements(prhs[0]);
        p  = mxMalloc(nc * sizeof *p);

        extract_rhs(prhs[0], p);
        max_blk = max_block(nc, p);

        if (partition_allocate_inverse(nc, max_blk, &pb2c, &b2c)) {
            plhs[0] = mxCreateNumericMatrix(max_blk + 1 + 1, 1, mxINT32_CLASS, mxREAL);
            plhs[1] = mxCreateNumericMatrix(nc,              1, mxINT32_CLASS, mxREAL);
            if (nlhs == 3) {
                plhs[2] = mxCreateNumericMatrix(nc,          1, mxINT32_CLASS, mxREAL);
            }

            partition_invert(nc, p, pb2c, b2c);

            assign_int_vec(pb2c, plhs[0]);
            assign_int_vec(b2c , plhs[1]);

            if (nlhs == 3) {
                partition_localidx(max_blk + 1, pb2c, b2c,
                                   mxGetData(plhs[2]));
                adjust_one_based_idx(nc, mxGetData(plhs[2]));
            }
        } else {
            plhs[0] = mxCreateDoubleScalar(mxGetNaN());
            plhs[1] = mxCreateDoubleScalar(mxGetNaN());
            if (nlhs == 3) {
                plhs[2] = mxCreateDoubleScalar(mxGetNaN());
            }
        }

        partition_deallocate_inverse(pb2c, b2c);
        mxFree(p);              /* p != NULL guaranteed here */
    } else {
        sprintf(errmsg,
                "Calling sequence is\n"
                "\t[pb2c, b2c]      = %s(p) %% or\n"
                "\t[pb2c, b2c, loc] = %s(p)",
                mexFunctionName(), mexFunctionName());

        mexErrMsgTxt(errmsg);
    }
}
