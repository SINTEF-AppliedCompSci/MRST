#include <assert.h>
#include <string.h>
#include <mex.h>

#include "grid.h"
#include "mrst_api.h"
#include "mrst_objects.h"
#include "trans_tpfa.h"


/* ------------------------------------------------------------------ */
static int
args_ok(int nlhs, int nrhs, const mxArray *prhs[])
/* ------------------------------------------------------------------ */
{
    int grid_ok = 0, rock_ok = 0;

    if (nrhs > 0) {
        grid_ok = verify_mrst_grid(prhs[0]);
    }

    if (nrhs > 1) {
        /* rock must contain a field 'perm' */
        rock_ok = mxIsStruct(prhs[1]) &&
                  (mxGetFieldNumber(prhs[1], "perm") >= 0);
    }

    return (nlhs == 1) && grid_ok && rock_ok;
}


/*
 * T = mex_compute_trans(G, rock)
 */

/* ------------------------------------------------------------------ */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ------------------------------------------------------------------ */
{
    size_t  totconn;
    grid_t *g;
    double *perm;
    char    errmsg[1023 + 1];

    double *htrans;

    if (args_ok(nlhs, nrhs, prhs)) {
        g       = mrst_grid(prhs[0]);
        perm    = NULL;
	totconn = 0;

        if (g != NULL) {
            perm    = mrst_perm(g->dimensions, prhs[1]);
            totconn = g->cell_facepos[ g->number_of_cells ];
        }        

        /* Create return values */
        plhs[0] = mxCreateDoubleMatrix(totconn, 1, mxREAL);

        /* Compute one-sided transmissibilities for all cells */
        htrans  = mxGetPr(plhs[0]); /* mxCreateDoubleMatrix == ZEROS */

        if (perm != NULL) {
            tpfa_htrans_compute(g, perm, htrans);
        }

        mxFree(perm);
        free_mrst_grid(g);
    } else {
        sprintf(errmsg,
                "Calling sequence is\n"
                "\tT = %s(G, rock)\n",
                mexFunctionName());

        mexErrMsgTxt(errmsg);
    }
}
