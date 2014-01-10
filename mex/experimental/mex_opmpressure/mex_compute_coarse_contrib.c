#include <stddef.h>
#include <string.h>

#include <mex.h>

#include "partition.h"
#include "coarse_sys.h"


/* ---------------------------------------------------------------------- */
static int
args_ok(int nlhs, int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int i, ok;

    ok = (nlhs == 2) && (nrhs == 6);
    ok = ok && mxIsDouble(prhs[0]);
    ok = ok && mxIsDouble(prhs[1]);

    for (i = 2; ok && (i < nrhs); i++) {
        ok = mxIsDouble(prhs[i]) || mxIsInt32(prhs[i]);
    }

    return ok;
}


/* ---------------------------------------------------------------------- */
static void
extract_int_vec(const mxArray *M_v, int *v)
/* ---------------------------------------------------------------------- */
{
    size_t i, n;

    int    *pi;
    double *pd;

    n = mxGetNumberOfElements(M_v);

    if (mxIsDouble(M_v)) {
        pd = mxGetPr(M_v);

        for (i = 0; i < n; i++) { v[i] = pd[i]; }
    } else {
        pi = mxGetData(M_v);

        for (i = 0; i < n; i++) { v[i] = pi[i]; }
    }
}


/* ---------------------------------------------------------------------- */
static void
add_constant(size_t n, int a, int *v)
/* ---------------------------------------------------------------------- */
{
    size_t i;

    for (i = 0; i < n; i++) { v[i] += a; }
}


/* ---------------------------------------------------------------------- */
static int
max_element(size_t n, int *v)
/* ---------------------------------------------------------------------- */
{
    int    ret;
    size_t i;

    ret = v[0];
    for (i = 1; i < n; i++) {
        ret = (v[i] > ret) ? v[i] : ret;
    }

    return ret;
}


/* ---------------------------------------------------------------------- */
static int
max_diff(size_t n, int *p)
/* ---------------------------------------------------------------------- */
{
    int    d, ret;
    size_t i;

    ret = -1;
    for (i = 0; i < n; i++) {
        d   = p[i + 1] - p[i];
        ret = (d > ret) ? d : ret;
    }

    return ret;
}


/* ---------------------------------------------------------------------- */
static int
extract_partition(const mxArray *M_p, int *nc, int *nb, int **p)
/* ---------------------------------------------------------------------- */
{
    int ret;

    *nc = mxGetNumberOfElements(M_p);
    *p  = mxMalloc(*nc * sizeof **p);

    if (*p != NULL) {
        extract_int_vec(M_p, *p);
        add_constant(*nc, -1, *p);

        *nb = max_element(*nc, *p) + 1;

        ret = 1;
    } else {
        ret = 0;
    }

    return ret;
}


/*
 * [cell_ip, Binv] = mex_compute_coarse_contrib(BIf, Psi, p, pconn, ...
 *                                              dof_pos, blk_ncf)
 */
/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok, nc, nb, b, i, tmp, max_bcells, max_nconn, max_nbpair;
    int *p, *pconn, *blk_ncf, *b2c_pos, *b2c;

    size_t totblkdof;

    char errmsg[1023 + 1];

    double *BI, *work, *totmob;
    
    struct coarse_sys sys;

    if (args_ok(nlhs, nrhs, prhs)) {
        BI        = mxGetPr(prhs[0]);
        sys.basis = mxGetPr(prhs[1]);
        ok = extract_partition(prhs[2], &nc, &nb, &p);

        if (ok) {
            ok          = partition_allocate_inverse(nc, nb - 1,
                                                     &b2c_pos, &b2c);

            pconn       = mxMalloc((nc + 1) * sizeof *pconn      );
            blk_ncf     = mxMalloc( nb      * sizeof *blk_ncf    );
            sys.dof_pos = mxMalloc((nb + 1) * sizeof *sys.dof_pos);
            sys.bf_pos  = mxCalloc( nb + 1  , sizeof *sys.bf_pos );
            sys.ip_pos  = mxCalloc( nb + 1  , sizeof *sys.ip_pos );

            if (ok && (pconn != NULL) && (blk_ncf != NULL) &&
                (sys.dof_pos != NULL) && (sys.bf_pos != NULL) &&
                (sys.ip_pos != NULL)) {
                partition_invert(nc, p, b2c_pos, b2c);

                extract_int_vec(prhs[3], pconn      );
                extract_int_vec(prhs[4], sys.dof_pos);
                extract_int_vec(prhs[5], blk_ncf    );

                add_constant(nc + 1, -1, pconn);
                add_constant(nb + 1, -1, sys.dof_pos);

                max_nconn  = max_diff(nc, pconn);
                max_bcells = max_diff(nb, b2c_pos);
                max_nbpair = max_diff(nb, sys.dof_pos);
                max_nbpair = max_nbpair * (max_nbpair + 1) / 2;

                totblkdof = 0;
                for (b = 1; b <= nb; b++) {
                    tmp = sys.dof_pos[b] - sys.dof_pos[b - 1];

                    totblkdof += tmp * tmp;
                    
                    sys.bf_pos[b] = sys.bf_pos[b - 1] + tmp*blk_ncf[b - 1];

                    tmp = tmp * (tmp + 1) / 2;
                    sys.ip_pos[b]  = sys.ip_pos[b - 1] + tmp*(b2c_pos[b] - b2c_pos[b - 1]);
                }

                plhs[0] = mxCreateDoubleMatrix(sys.ip_pos[nb], 1, mxREAL);
                plhs[1] = mxCreateDoubleMatrix(totblkdof     , 1, mxREAL);
                work    = mxMalloc((max_bcells + max_nbpair*(max_nbpair + 1)/2) *
                                   sizeof *work);
                totmob  = mxMalloc(nc * sizeof *totmob);

                if ((plhs[0] != NULL) && (plhs[1] != NULL) &&
                    (work != NULL) && (totmob != NULL)) {
                    sys.cell_ip = mxGetPr(plhs[0]);
                    sys.Binv    = mxGetPr(plhs[1]);

                    for (i = 0; i < nc; i++) { totmob[i] = 1.0; }

                    coarse_sys_compute_cell_ip(nc, max_nconn, nb, pconn, BI,
                                               b2c_pos, b2c, &sys);
                    coarse_sys_compute_Binv(nb, max_bcells, totmob,
                                            b2c_pos, b2c, &sys, work);
                }

                if (totmob != NULL) { mxFree(totmob); }
                if (work   != NULL) { mxFree(work);   }
            }

            if (sys.ip_pos  != NULL) { mxFree(sys.ip_pos);  }
            if (sys.bf_pos  != NULL) { mxFree(sys.bf_pos);  }
            if (sys.dof_pos != NULL) { mxFree(sys.dof_pos); }
            if (blk_ncf     != NULL) { mxFree(blk_ncf);     }
            if (pconn       != NULL) { mxFree(pconn);       }

            partition_deallocate_inverse(b2c_pos, b2c);
        }

        if (p != NULL) { mxFree(p); }
    } else {
        sprintf(errmsg,
                "Calling sequence is\n\t"
                "[cell_ip, Binv] = %s(BIf, Psi, p, pconn, dof_pos, blk_ncf)",
                mexFunctionName());

        mexErrMsgTxt(errmsg);
    }
}
