#include <string.h>

#include <mex.h>
#define MAT_SIZE_T mwSignedIndex

#include "call_umfpack.h"
#include "hybsys.h"
#include "mrst_objects.h"

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


/*
 * [S, r, F, L, q] = mex_schur_comp_symm(BI, connPos, conns, src)
 */
/* ---------------------------------------------------------------------- */
static int
verify_args(int nlhs, int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok;

    ok =       (nlhs == 5) && (nrhs == 4);
    ok = ok && mxIsDouble(prhs[0]);
    ok = ok && (mxIsDouble(prhs[1]) || mxIsInt32(prhs[1]));
    ok = ok && (mxGetNumberOfElements(prhs[0]) >
                mxGetNumberOfElements(prhs[1]));
    ok = ok && verify_src (prhs[3]);                          /* src */    

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
deallocate_aux_arrays(int *pconn, double *gflux)
/* ---------------------------------------------------------------------- */
{
    /* Apparently mxFree() makes no guarantee regarding NULL arguments.
     * Institute a belt-and-suspenders approach to releasing resources.
     */
    if (gflux != NULL) { mxFree(gflux); }
    if (pconn != NULL) { mxFree(pconn); }
}


/* ---------------------------------------------------------------------- */
static int
allocate_aux_arrays(int nconn_tot, double **gflux)
/* ---------------------------------------------------------------------- */
{
    int    ret;
    double *g;
   
    g = mxMalloc(nconn_tot * sizeof *g);

    if ((g == NULL)) {
        deallocate_aux_arrays(NULL,  g);
     
        *gflux = NULL;

        ret = 0;
    } else {
        *gflux = g;

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


#if defined(ASSEMBLE_AND_SOLVE_UMFPACK) && ASSEMBLE_AND_SOLVE_UMFPACK
/* ---------------------------------------------------------------------- */
static void
get_conn(const mxArray *M_conn, int *conn)
/* ---------------------------------------------------------------------- */
{
    size_t nel, i;

    int    *pi;
    double *pd;

    nel = mxGetNumberOfElements(M_conn);

    if (mxIsDouble(M_conn)) {
        pd = mxGetPr(M_conn);

        for (i = 0; i < nel; i++) { conn[i] = pd[i] - 1; }
    } else {
        pi = mxGetData(M_conn);

        for (i = 0; i < nel; i++) { conn[i] = pi[i] - 1; }
    }
}

/* ---------------------------------------------------------------------- */
static int
get_number_of_faces(int nc, int *pconn, int *conn)
/* ---------------------------------------------------------------------- */
{
   int N, i, nf;

   N = pconn[nc];
   nf = 0;
   for(i=0; i<N; ++i)
   {
      nf = MAX(nf, conn[i]+1);
   }

   return nf;
}
#endif


/*
 * [S, r, F, L, q] = mex_schur_comp_symm(BI, connPos, conns, src)
 */

/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok, nc, nconn_tot, max_nconn, sum_nconn2;
    int p2, c, i, nconn, *pconn;
    double *Binv, *ptr1, *ptr2, *src, *gpress;
    struct hybsys *sys;

#if defined(ASSEMBLE_AND_SOLVE_UMFPACK) && ASSEMBLE_AND_SOLVE_UMFPACK
    int *conn, nf;
    double *b, *x;
    struct Sparse A;
#endif

    char errmsg[1023 + 1];

    ok = verify_args(nlhs, nrhs, prhs);

    if (ok) {
        nc = get_pconn(prhs[1], &pconn);

        nconn_tot = pconn[nc];
        max_nconn = count_cellconn(nc, pconn);

        /* allocate_aux_arrays(nc, nconn_tot, &src, &gpress); */
        allocate_aux_arrays(nconn_tot, &gpress);
        
        sum_nconn2 = mxGetNumberOfElements(prhs[0]);
        plhs[0]    = mxCreateDoubleMatrix(sum_nconn2, 1, mxREAL); /* S */
        plhs[1]    = mxCreateDoubleMatrix(nconn_tot,  1, mxREAL); /* r */
        plhs[2]    = mxCreateDoubleMatrix(nconn_tot,  1, mxREAL); /* F */
        plhs[3]    = mxCreateDoubleMatrix(nc,         1, mxREAL); /* L */
        plhs[4]    = mxCreateDoubleMatrix(nc,         1, mxREAL); /* q */

        sys = hybsys_allocate_symm(max_nconn, nc, nconn_tot);
        hybsys_init(max_nconn, sys);
 
        src    = mrst_src(nc, prhs[3]); 
        for (i = 0; i < nconn_tot; i++) { gpress[i] = 0.0; } /* No gravity */

        Binv = mxGetPr(prhs[0]);

        hybsys_schur_comp_symm(nc, pconn, Binv, sys);

        ptr1 = mxGetPr(plhs[0]); /* S */
        ptr2 = mxGetPr(plhs[1]); /* r */
        p2 = 0;
        for (c = 0; c < nc; c++) {
            nconn = pconn[c + 1] - pconn[c];

            hybsys_cellcontrib_symm(c, nconn, pconn[c], p2,
                                    gpress, src, Binv, sys);

            memcpy(ptr1 + p2      , sys->S, nconn * nconn * sizeof *ptr1);
            memcpy(ptr2 + pconn[c], sys->r, nconn         * sizeof *ptr2);

            p2 += nconn * nconn;
        }

#if defined(ASSEMBLE_AND_SOLVE_UMFPACK) && ASSEMBLE_AND_SOLVE_UMFPACK
        conn = mxMalloc(mxGetNumberOfElements(prhs[2]) * sizeof *conn);
        get_conn(prhs[2], conn);

        nf          = get_number_of_faces(nc, pconn, conn);
        hybsys_assemble(nc, nf, pconn, conn, ptr, sys->r, &A, &b);
        x = mxMalloc(A.n * sizeof *x);

        callMWUMFPACK(A.n, A.ia, A.ja, A.sa, b, x);

        for (i = 0; i < nf; ++i) mexPrintf("x[%d] = %f\n", i, x[i]);
        mxFree(x);  free(b);
        free(A.ia); free(A.ja); free(A.sa);

        mxFree(conn);
#endif

        ptr1 = mxGetPr(plhs[2]); /* F */
        memcpy(ptr1, sys->F1, nconn_tot * sizeof *ptr1);

        ptr1 = mxGetPr(plhs[3]); /* L */
        memcpy(ptr1, sys->L , nc        * sizeof *ptr1);

        ptr1 = mxGetPr(plhs[4]); /* q */
        memcpy(ptr1, sys->q , nc        * sizeof *ptr1);

        hybsys_free(sys);
        deallocate_aux_arrays(pconn, gpress);
        mrst_src_deallocate (src);  
    } else {
        sprintf(errmsg,
                "Calling sequence is:\n\t"
                "[S, r, F, L, q] = %s(BI, connPos, conns, src)",
                mexFunctionName());
        mexErrMsgTxt(errmsg);
    }
}
