#include <string.h>

#include <mex.h>
#define MAT_SIZE_T mwSignedIndex

#include "hybsys.h"

#define MAX(a,b) (((a) > (b)) ? (a) : (b))

/* ---------------------------------------------------------------------- */
static int
verify_args(int nlhs, int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok;

    ok =       (nlhs == 5) && (nrhs == 3);
    ok = ok && mxIsDouble(prhs[0]);
    ok = ok && (mxIsDouble(prhs[1]) || mxIsInt32(prhs[1]));
    ok = ok && (mxGetNumberOfElements(prhs[0]) >
                mxGetNumberOfElements(prhs[1]));

    return ok;
}


/* ---------------------------------------------------------------------- */
static int
count_cellconn(int nc, const int *pconn)
/* ---------------------------------------------------------------------- */
{
    int c, nconn, max_nconn;

    max_nconn = 0;

    for (c = 0; c < nc; c++) {
        nconn = pconn[c + 1] - pconn[c];

        max_nconn = MAX(max_nconn, nconn);
    }

    return max_nconn;
}


/* ---------------------------------------------------------------------- */
static void
deallocate_aux_arrays(int *pconn, double *src, double *gflux,
                      double *BIV, double *P)
/* ---------------------------------------------------------------------- */
{
    /* Apparently mxFree() makes no guarantee regarding NULL arguments.
     * Institute a belt-and-suspenders approach to releasing resources.
     */
    if (P     != NULL) { mxFree(P);     }
    if (BIV   != NULL) { mxFree(BIV);   }
    if (gflux != NULL) { mxFree(gflux); }
    if (src   != NULL) { mxFree(src);   }
    if (pconn != NULL) { mxFree(pconn); }
}


/* ---------------------------------------------------------------------- */
static int
allocate_aux_arrays(int nc, int nconn_tot,
                    double **src, double **gflux,
                    double **BIV, double **P)
/* ---------------------------------------------------------------------- */
{
    int    ret;
    double *s, *g, *B, *p;

    s = mxMalloc(nc        * sizeof *s);
    g = mxMalloc(nconn_tot * sizeof *g);
    B = mxMalloc(nconn_tot * sizeof *B);
    p = mxMalloc(nc        * sizeof *p);

    if ((s == NULL) || (g == NULL) || (B == NULL) || (p == NULL)) {
        deallocate_aux_arrays(NULL, s, g, B, p);

        *src   = NULL;
        *gflux = NULL;
        *BIV   = NULL;
        *P     = NULL;

        ret = 0;
    } else {
        *src   = s;
        *gflux = g;
        *BIV   = B;
        *P     = p;

        ret = 1;
    }

    return ret;
}


/* ---------------------------------------------------------------------- */
static int
get_pconn(const mxArray *M_pconn, int **pconn)
/* ---------------------------------------------------------------------- */
{
    int ret;

    size_t e, ne;

    int    *pi;
    double *pd;

    ne = mxGetNumberOfElements(M_pconn);

    *pconn = mxMalloc(ne * sizeof **pconn);

    if (*pconn != NULL) {
        if (mxIsDouble(M_pconn)) {
            pd = mxGetPr(M_pconn);

            for (e = 0; e < ne; e++) { (*pconn)[e] = pd[e] - 1; }
        } else {
            pi = mxGetData(M_pconn);

            for (e = 0; e < ne; e++) { (*pconn)[e] = pi[e] - 1; };
        }

        ret = ne - 1;
    } else {
        ret = -1;
    }

    return ret;
}


/*
 * [S, r, F1, F2, L] = mex_schur_comp(BI, connPos, conns)
 */

/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok, nc, nconn_tot, max_nconn, sum_nconn2;
    int p2, c, i, nconn, *pconn;
    double *Binv, *ptr1, *ptr2, *src, *gpress, *BIV, *P;
    struct hybsys *sys;

    ok = verify_args(nlhs, nrhs, prhs);

    if (ok) {
        nc = get_pconn(prhs[1], &pconn);

        nconn_tot = pconn[nc];
        max_nconn = count_cellconn(nc, pconn);

        allocate_aux_arrays(nc, nconn_tot, &src, &gpress, &BIV, &P);

        sum_nconn2 = mxGetNumberOfElements(prhs[0]);
        plhs[0]    = mxCreateDoubleMatrix(sum_nconn2, 1, mxREAL);
        plhs[1]    = mxCreateDoubleMatrix(nconn_tot,  1, mxREAL);
        plhs[2]    = mxCreateDoubleMatrix(nconn_tot,  1, mxREAL);
        plhs[3]    = mxCreateDoubleMatrix(nconn_tot,  1, mxREAL);
        plhs[4]    = mxCreateDoubleMatrix(nc,         1, mxREAL);

        sys = hybsys_allocate_unsymm(max_nconn, nc, nconn_tot);
        hybsys_init(max_nconn, sys);

        for (i = 0; i < nc; i++)        { src[i]    = 0.0;   /* No sources */
                                          P[i]      = 0.0; } /* No accumulation */
        for (i = 0; i < nconn_tot; i++) { gpress[i] = 0.0;   /* No gravity */
                                          BIV[i]    = 0.0; } /* No V^2 */

        Binv = mxGetPr(prhs[0]);

        hybsys_schur_comp_unsymm(nc, pconn, Binv, BIV, P, sys);

        ptr1 = mxGetPr(plhs[0]);
        ptr2 = mxGetPr(plhs[1]);
        p2 = 0;
        for (c = 0; c < nc; c++) {
            nconn = pconn[c + 1] - pconn[c];

            hybsys_cellcontrib_unsymm(c, nconn, pconn[c], p2,
                                      gpress, src, Binv, sys);

            memcpy(ptr1 + p2      , sys->S, nconn * nconn * sizeof *ptr1);
            memcpy(ptr2 + pconn[c], sys->r, nconn         * sizeof *ptr2);

            p2 += nconn * nconn;
        }

        ptr1 = mxGetPr(plhs[2]);
        memcpy(ptr1, sys->F1, nconn_tot * sizeof *ptr1);

        ptr1 = mxGetPr(plhs[3]);
        memcpy(ptr1, sys->F2, nconn_tot * sizeof *ptr1);

        ptr1 = mxGetPr(plhs[4]);
        memcpy(ptr1, sys->L , nc        * sizeof *ptr1);

        hybsys_free(sys);
        deallocate_aux_arrays(pconn, src, gpress, BIV, P);
    }
}
