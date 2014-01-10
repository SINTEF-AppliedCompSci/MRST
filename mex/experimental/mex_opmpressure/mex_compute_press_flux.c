#include <string.h>

#include "mex.h"

#include "mrst_objects.h"

#include "hybsys.h"


#define MAX(a,b) (((a) > (b)) ? (a) : (b))

/* ---------------------------------------------------------------------- */
static int
verify_src(const mxArray *src)
/* ---------------------------------------------------------------------- */
{
    int ok;

    ok = mxIsEmpty(src);        /* Support empty src */

    if (!ok) {
        ok =        mxIsStruct(src);
        ok = ok && (mxGetFieldNumber(src, "cell") >= 0);
        ok = ok && (mxGetFieldNumber(src, "rate") >= 0) &&
            mxIsDouble(mxGetField(src, 0, "rate"));
    }

    return ok;
}


/* [v, p] = mex_compute_press_flux(BI, pi, connPos, conns, F, L, src) */
/* ---------------------------------------------------------------------- */
static int
verify_args(int nlhs, int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok;

    ok =       (nlhs == 2) && (nrhs == 7);
    ok = ok && mxIsDouble(prhs[0]);                          /* BI */
    ok = ok && mxIsDouble(prhs[1]);                          /* lam */
    ok = ok && (mxIsDouble(prhs[2]) || mxIsInt32(prhs[2]));  /* connPos */
    ok = ok && (mxIsDouble(prhs[3]) || mxIsInt32(prhs[3]));  /* conns */
    ok = ok && mxIsDouble(prhs[4]);                          /* F */
    ok = ok && mxIsDouble(prhs[5]);                          /* L */
    ok = ok && (verify_src(prhs[6]) || mxIsDouble(prhs[6])); /* src */

    return ok;
}


/* ---------------------------------------------------------------------- */
static void
count_conns(const mxArray *M_pconn, int *nc, int *max_nconn, int *nconn_tot)
/* ---------------------------------------------------------------------- */
{
    int    c, nconn, *pi;
    double           *pd;

    *nc = mxGetNumberOfElements(M_pconn) - 1;

    *max_nconn = *nconn_tot = 0;

    if (mxIsDouble(M_pconn)) {
        pd = mxGetPr(M_pconn);

        for (c = 0; c < *nc; c++) {
            nconn = pd[c + 1] - pd[c];
            *max_nconn  = MAX(*max_nconn, nconn);
            *nconn_tot += nconn;
        }
    } else {
        pi = mxGetData(M_pconn);

        for (c = 0; c < *nc; c++) {
            nconn = pi[c + 1] - pi[c];
            *max_nconn  = MAX(*max_nconn, nconn);
            *nconn_tot += nconn;
        }
    }
}


/* ---------------------------------------------------------------------- */
static void
deallocate_aux_arrays(int *pconn, int *conn,
                      double *gflux, double *work)
/* ---------------------------------------------------------------------- */
{
    /* Apparently mxFree() makes no guarantee regarding NULL arguments.
     * Institute a belt-and-suspenders approach to releasing resources.
     */
    if (work  != NULL) { mxFree(work);  }
    if (gflux != NULL) { mxFree(gflux); }
    if (conn  != NULL) { mxFree(conn);  }
    if (pconn != NULL) { mxFree(pconn); }
}


/* ---------------------------------------------------------------------- */
static int
allocate_aux_arrays(int max_nconn, int nc, int nconn_tot,
                    int **pconn, int **conn,
                    double **gflux, double **work)
/* ---------------------------------------------------------------------- */
{
    int    ret, *p, *c;
    double *g, *w;

    p = mxMalloc((nc + 1)  * sizeof *p);
    c = mxMalloc(nconn_tot * sizeof *c);
    g = mxMalloc(nconn_tot * sizeof *g);
    w = mxMalloc(max_nconn * sizeof *w);

    if ((p == NULL) || (c == NULL) || (g == NULL) || (w == NULL)) {
        deallocate_aux_arrays(p, c, g, w);

        *pconn = NULL;
        *conn  = NULL;
        *gflux = NULL;
        *work  = NULL;

        ret = 0;
    } else {
        *pconn = p;
        *conn  = c;
        *gflux = g;
        *work  = w;

        ret = 1;
    }

    return ret;
}


/* ---------------------------------------------------------------------- */
static void
copy_M_int_vector(const mxArray *M_a, int *a)
/* ---------------------------------------------------------------------- */
{
    size_t nel, i;

    int    *pi;
    double *pd;

    nel = mxGetNumberOfElements(M_a);

    if (mxIsDouble(M_a)) {
        pd = mxGetPr(M_a);

        for (i = 0; i < nel; i++) { a[i] = pd[i] - 1; }
    } else {
        pi = mxGetData(M_a);

        for (i = 0; i < nel; i++) { a[i] = pi[i] - 1; }
    }
}


/* ---------------------------------------------------------------------- */
static void
get_pconn(const mxArray *M_pconn, int *pconn)
/* ---------------------------------------------------------------------- */
{
    copy_M_int_vector(M_pconn, pconn);
}


/* ---------------------------------------------------------------------- */
static void
get_conn(const mxArray *M_conn, int *conn)
/* ---------------------------------------------------------------------- */
{
    copy_M_int_vector(M_conn, conn);
}


/* ---------------------------------------------------------------------- */
static void
set_source_term(size_t nc, const mxArray *Q, double *q)
/* ---------------------------------------------------------------------- */
{
    size_t c;
    double *mq;

    if (mxIsEmpty(Q)) {
        for (c = 0; c < nc; c++) {
            q[c] = 0.0;
        }
    } else if (mxIsStruct(Q) && verify_src(Q)) {
        mq = mrst_src(nc, Q);

        if (mq != NULL) {
            memcpy(q, mq, nc * sizeof *mq);
        } else {
            for (c = 0; c < nc; c++) {
                q[c] = 0.0;
            }
        }

        mrst_src_deallocate(mq);
    } else {
        if (mxGetNumberOfElements(Q) != nc) {
            mexErrMsgTxt("Source terms must be specified in each "
                         "cell if provided as a single vector of "
                         "type DOUBLE.");
        }

        mq = mxGetPr(Q);
        memcpy(q, mq, nc * sizeof *mq);
    }
}
        
/*
 * [v, p] = mex_compute_press_flux(BI, pi, connPos, conns, F, L, src)
 */

/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok, i, nc, max_nconn, nconn_tot, *pconn, *conn;
    double *gflux, *work, *Binv, *pi, *ptr;

    struct hybsys *sys;

    char errmsg[1023 + 1];

    ok = verify_args(nlhs, nrhs, prhs);

    if (ok) {
        count_conns(prhs[2], &nc, &max_nconn, &nconn_tot);

        allocate_aux_arrays(max_nconn, nc, nconn_tot,
                            &pconn, &conn, &gflux, &work);

        plhs[0] = mxCreateDoubleMatrix(nconn_tot, 1, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(nc,        1, mxREAL);

        sys = hybsys_allocate_symm(max_nconn, nc, nconn_tot);
        hybsys_init(max_nconn, sys);

        set_source_term(nc, prhs[6], sys->q);
        
        ptr = mxGetPr(prhs[4]);
        memcpy(sys->F1, ptr, nconn_tot * sizeof *sys->F1);

        ptr = mxGetPr(prhs[5]);
        memcpy(sys->L , ptr, nc        * sizeof *sys->L);

        get_pconn(prhs[2], pconn);
        get_conn (prhs[3], conn);

        Binv = mxGetPr(prhs[0]);
        pi   = mxGetPr(prhs[1]);

        for (i = 0; i < nconn_tot; i++) { gflux[i] = 0.0; } /* No gravity */

        hybsys_compute_press_flux(nc, pconn, conn, gflux, Binv, sys,
                                  pi, mxGetPr(plhs[1]), mxGetPr(plhs[0]),
                                  work);

        hybsys_free(sys);
        deallocate_aux_arrays(pconn, conn, gflux, work);
    } else {
        sprintf(errmsg,
                "Calling sequence is:\n\t"
                "[v, p] = %s(BI, pi, connPos, conns, F, L, src)",
                mexFunctionName());
        mexErrMsgTxt(errmsg);
    }
}
