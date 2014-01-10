#include <stddef.h>
#include <string.h>

#include <mex.h>

#include "coarse_sys.h"

#include "mrst_msmfe_support.h"

/* ---------------------------------------------------------------------- */
static int
verify_totmob(const mxArray *cs, const mxArray *totmob)
/* ---------------------------------------------------------------------- */
{
    int    ok;
    size_t nc;

    nc = mxGetNumberOfElements(mxGetField(cs, 0, "b2c"));

    ok = mxIsDouble(totmob) && (mxGetNumberOfElements(totmob) == nc);

    return ok;
}


/* ---------------------------------------------------------------------- */
static int
args_ok(int nlhs, int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok;

    ok = (nlhs == 1) && (nrhs == 2);

    if (ok) {
        ok =       verify_mex_cs(prhs[0]);
        ok = ok && verify_totmob(prhs[0], prhs[1]);
    }

    return ok;
}


/*
 * BI = mex_compute_coarse_BI(mex_cs, totmob)
 */

/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    double                 *totmob, *work;
    struct coarse_sys_meta  csys_meta;
    struct coarse_sys       csys;
    char                    errmsg[1023 + 1];

    if (args_ok(nlhs, nrhs, prhs)) {
        set_coarse_sys_meta(prhs[0], &csys_meta);
        set_coarse_sys     (prhs[0], &csys     );

        plhs[0] = mxCreateDoubleMatrix(csys_meta.sum_ndof2, 1, mxREAL);
        work    = mxMalloc(csys_meta.work_size * sizeof *work);

        if ((plhs[0] != NULL) && (work != NULL)) {
            csys.Binv = mxGetPr(plhs[0]);
            totmob    = mxGetPr(prhs[1]);

            coarse_sys_compute_Binv(csys_meta.nblocks   ,
                                    csys_meta.max_bcells,
                                    totmob              ,
                                    csys_meta.b2c_pos   ,
                                    csys_meta.b2c, &csys, work);

            mxFree(work);
        }
    } else {
        sprintf(errmsg,
                "Calling sequence is:\n\t"
                "BI = %s(mex_cs, totmob)",
                mexFunctionName());
        mexErrMsgTxt(errmsg);
    }
}
