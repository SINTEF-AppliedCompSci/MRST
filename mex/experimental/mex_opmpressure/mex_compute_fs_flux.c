#include <stddef.h>
#include <string.h>

#include <mex.h>

#include "grid.h"
#include "coarse_conn.h"
#include "coarse_sys.h"

#include "mrst_api.h"
#include "mrst_msmfe_support.h"


/* ---------------------------------------------------------------------- */
static int
args_ok(int nlhs, int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok;

    ok = (nlhs == 1) && (nrhs == 3);

    if (ok) {
        ok =       verify_mrst_grid(prhs[0]);
        ok = ok && verify_mex_cs(prhs[1]);
        ok = ok && !mxIsEmpty(prhs[2]) && mxIsDouble(prhs[2]);
    }

    return ok;
}


/*
 * flux = mex_compute_fs_flux(G, mex_cs, block_hflux)
 */
/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    size_t                 ngconn;
    grid_t                 *g;
    double                 *v_c;
    double                 *work;

    struct coarse_topology ct;
    struct coarse_sys_meta csys_meta;
    struct coarse_sys      csys;

    char errmsg[1023 + 1];

    if (args_ok(nlhs, nrhs, prhs)) {
        set_coarse_sys_meta(prhs[1], &csys_meta);
        set_coarse_sys     (prhs[1], &csys     );

        g   = mrst_grid(prhs[0]);
        v_c = mxGetPr  (prhs[2]);

        ngconn = g->cell_facepos[ g->number_of_cells ];
        work = mxMalloc(ngconn * sizeof *work);

        if ((g != NULL) && (work != NULL)) {
            plhs[0] = mxCreateDoubleMatrix(g->number_of_faces, 1, mxREAL);

            /* This is, generally, very dangerous.  However, the
             * coarse_sys_compute_fs_flux() function only accesses the
             * 'nblocks' field and it is better to have a reproducible
             * crash if *that* ever changes... */
            ct.nblocks    = csys_meta.nblocks;
            ct.nfaces     = 0;
            ct.neighbours = NULL;
            ct.blkfacepos = NULL;   ct.blkfaces = NULL;
            ct.subfacepos = NULL;   ct.subfaces = NULL;

            coarse_sys_compute_fs_flux(g, &ct, &csys,
                                       csys_meta.b2c_pos, csys_meta.b2c,
                                       v_c, mxGetPr(plhs[0]), work);

            mxFree(work);
            free_mrst_grid(g);
        }
    } else {
        sprintf(errmsg,
                "Calling sequence is:\n\t"
                "flux = %s(G, mex_cs, block_hflux)",
                mexFunctionName());
        mexErrMsgTxt(errmsg);
    }
}
