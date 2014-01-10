#include <stddef.h>
#include <string.h>

#include <mex.h>

#include "grid.h"

#include "mrst_api.h"
#include "mrst_objects.h"
#include "call_umfpack.h"

#include "flow_bc.h"
#include "well.h"

#include "sparse_sys.h"

#include "compr_quant.h"
#include "trans_tpfa.h"
#include "cfs_tpfa.h"


/* ---------------------------------------------------------------------- */
static int
verify_state(const mxArray *state)
/* ---------------------------------------------------------------------- */
{
    return mxIsStruct(state) &&
        (mxGetFieldNumber(state, "pressure") >= 0);
}


/* ------------------------------------------------------------------ */
static int
verify_faces_structure(mxArray *faces)
/* ------------------------------------------------------------------ */
{
    /* Shallow structural inspection only.  Assume valid fields... */
    int ok;

    ok =       (mxGetFieldNumber(faces, "neighbors") >= 0);
    ok = ok && (mxGetFieldNumber(faces, "areas"    ) >= 0);
    ok = ok && (mxGetFieldNumber(faces, "normals"  ) >= 0);
    ok = ok && (mxGetFieldNumber(faces, "centroids") >= 0);

    return ok;
}


/* ------------------------------------------------------------------ */
static int
verify_cells_structure(mxArray *cells)
/* ------------------------------------------------------------------ */
{
    /* Shallow structural inspection only.  Assume valid fields... */
    int ok;

    ok =       (mxGetFieldNumber(cells, "facePos"  ) >= 0);
    ok = ok && (mxGetFieldNumber(cells, "faces"    ) >= 0);
    ok = ok && (mxGetFieldNumber(cells, "volumes"  ) >= 0);
    ok = ok && (mxGetFieldNumber(cells, "centroids") >= 0);

    return ok;
}


/* ---------------------------------------------------------------------- */
static int
verify_grid(const mxArray *G)
/* ---------------------------------------------------------------------- */
{
    int nodes_ok = 0, faces_ok = 0, cells_ok = 0, field_no;

    mxArray *pm;

    if (mxIsStruct(G)) {
        nodes_ok = mxGetFieldNumber(G, "nodes") >= 0;

        field_no = mxGetFieldNumber(G, "faces");
        faces_ok = field_no >= 0;
        if (faces_ok) {
            pm = mxGetFieldByNumber(G, 0, field_no);
            faces_ok = verify_faces_structure(pm);
        }

        field_no = mxGetFieldNumber(G, "cells");
        cells_ok = field_no >= 0;
        if (cells_ok) {
            pm = mxGetFieldByNumber(G, 0, field_no);
            cells_ok = verify_cells_structure(pm);
        }
    }

    return nodes_ok && faces_ok && cells_ok;
}


/* ---------------------------------------------------------------------- */
static int
verify_htrans(const mxArray *G, const mxArray *T)
/* ---------------------------------------------------------------------- */
{
    int    ok;
    size_t ncf;

    ncf = getNumberOfCellFaces(G);

    ok  =       !mxIsEmpty(T) && mxIsDouble(T);
    ok  = ok && (mxGetNumberOfElements(T) == ncf);

    return ok;
}


/* ---------------------------------------------------------------------- */
static int
verify_trans(const mxArray *G, const mxArray *ft)
/* ---------------------------------------------------------------------- */
{
    int    ok;
    size_t nf;

    nf = getNumberOfFaces(G);

    ok =       !mxIsEmpty(ft) && mxIsDouble(ft);
    ok = ok && (mxGetNumberOfElements(ft) == nf);

    return ok;
}


/*
 * fq.np      = number of phases
 * fq.A.c     = (R / B)_c
 * fq.A.f     = (R / B)_f
 * fq.ct      = total compressibility (per cell)
 * fq.perfmob = phase mobility per well perforation
 * fq.vd      = volume discrepancy term (per cell)
 * fq.pmobf   = phase mobility per face
 */
/* ---------------------------------------------------------------------- */
static int
verify_fluid_quant(const mxArray *G, const mxArray *fq)
/* ---------------------------------------------------------------------- */
{
    int    ok;
    size_t nc, nf, np;

    const mxArray *A, *Ac, *Af;

    ok =        !mxIsEmpty(fq) && mxIsStruct(fq);
    ok = ok && (mxGetFieldNumber(fq, "np"      ) >= 0);
    ok = ok && (mxGetFieldNumber(fq, "A"       ) >= 0);
    ok = ok && (mxGetFieldNumber(fq, "ct"      ) >= 0);
    ok = ok && (mxGetFieldNumber(fq, "perfmob" ) >= 0);
    ok = ok && (mxGetFieldNumber(fq, "vd"      ) >= 0);
    ok = ok && (mxGetFieldNumber(fq, "pmobf"   ) >= 0);

    if (ok) {
        nc = getNumberOfCells(G);
        nf = getNumberOfFaces(G);

        A  = mxGetField(fq, 0, "A");
        ok = ok && !mxIsEmpty(A) && mxIsStruct(A);

        Ac = mxGetField(A, 0, "c");
        Af = mxGetField(A, 0, "f");

        ok = ok && (Ac != NULL) && (Af != NULL);

        if (ok) {
            mxAssert (mxGetNumberOfElements(mxGetField(fq, 0, "np")) == 1,
                      "Number of phases, fq.np, must be scalar.");

            np = (size_t) mxGetScalar(mxGetField(fq, 0, "np"));

            ok = ok && mxIsDouble(Ac) && mxIsDouble(Af);

            ok = ok && (mxGetNumberOfElements(Ac) == np * np * nc);
            ok = ok && (mxGetNumberOfElements(Af) == np * np * nf);
        }
    }

    return ok;
}


/* ---------------------------------------------------------------------- */
static int
verify_grav_terms(const mxArray *G, const mxArray *fq, const mxArray *Zf)
/* ---------------------------------------------------------------------- */
{
    int    ok;
    size_t nf, np;

    ok = !mxIsEmpty(Zf) && mxIsDouble(Zf);

    if (ok) {
        nf = getNumberOfFaces(G);
        np = (size_t) mxGetScalar(mxGetField(fq, 0, "np"));

        ok = mxGetNumberOfElements(Zf) == nf * np;
    }

    return ok;
}


/* ---------------------------------------------------------------------- */
static int
verify_timestep(const mxArray *dt)
/* ---------------------------------------------------------------------- */
{
    int ok;

    ok =       !mxIsEmpty(dt) && mxIsDouble(dt);
    ok = ok && (mxGetNumberOfElements(dt) == 1);
    ok = ok && (mxGetScalar(dt) > 0);

    return ok;
}


/* ---------------------------------------------------------------------- */
static int
verify_well(const mxArray *W)
/* ---------------------------------------------------------------------- */
{
    int ok;
    ok = mxIsEmpty(W);          /* Support empty well */

    if (!ok) {
        ok =        mxIsStruct(W);
        ok = ok && (mxGetFieldNumber(W, "cells") >= 0);
        ok = ok && (mxGetFieldNumber(W, "type" ) >= 0);
        ok = ok && (mxGetFieldNumber(W, "val"  ) >= 0);
        ok = ok && (mxGetFieldNumber(W, "WI"   ) >= 0);
        ok = ok && (mxGetFieldNumber(W, "dZ"   ) >= 0);
    }

    return ok;
}


/* ---------------------------------------------------------------------- */
static int
verify_bc(const mxArray *bc)
/* ---------------------------------------------------------------------- */
{
    int ok;

    ok = mxIsEmpty(bc);         /* Support empty bc */

    if (!ok) {
        ok =        mxIsStruct(bc);
        ok = ok && (mxGetFieldNumber(bc, "face" ) >= 0);
        ok = ok && (mxGetFieldNumber(bc, "type" ) >= 0) &&
            mxIsCell(mxGetField(bc, 0, "type"));
        ok = ok && (mxGetFieldNumber(bc, "value") >= 0) &&
            mxIsDouble(mxGetField(bc, 0, "value"));
    }

    return ok;
}


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
 * [x, Tf, Gf, wbhp, wflux] = mex_cfs_tpfa(x, G, T, ft, fq,
 *                                         Zf, dt, p0, pvol,
 *                                         W, bc, src)
 *
 * fq.np    = number of phases
 * fq.A.c   = (R / B)_c
 * fq.A.f   = (R / B)_f
 * fq.ct    = total compressibility (per cell)
 * fq.Lt    = total mobility (per cell)
 * fq.vd    = volume discrepancy (per cell)
 * fq.pmobf = phase mobility per face */
/* ---------------------------------------------------------------------- */
static int
args_ok(int nlhs, int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok;

    ok = (nlhs == 5) && (nrhs == 12);

    ok = ok && verify_state      (prhs[0]);                   /* x */
    ok = ok && verify_grid       (prhs[1]);                   /* G */
    ok = ok && verify_htrans     (prhs[1], prhs[2]);          /* T */
    ok = ok && verify_trans      (prhs[1], prhs[3]);          /* ft */
    ok = ok && verify_fluid_quant(prhs[1], prhs[4]);          /* fq */
    ok = ok && verify_grav_terms (prhs[1], prhs[4], prhs[5]); /* Zf */
    ok = ok && verify_timestep   (prhs[6]);                   /* dt */
    ok = ok && mxIsDouble        (prhs[7]);                   /* p0 */
    ok = ok && mxIsDouble        (prhs[8]);                   /* pvol */
    ok = ok && verify_well       (prhs[9]);                   /* W */
    ok = ok && verify_bc         (prhs[10]);                  /* bc */
    ok = ok && verify_src        (prhs[11]);                  /* src */

    return ok;
}


/* ---------------------------------------------------------------------- */
static void
set_mrst_cq(const mxArray *M_fq, struct compr_quantities *cq)
/* ---------------------------------------------------------------------- */
{
    const mxArray *A;

    cq->nphases   = mxGetScalar(mxGetField(M_fq, 0, "np"   ));
    A             =             mxGetField(M_fq, 0, "A"    ) ;
    cq->Ac        = mxGetPr    (mxGetField(A   , 0, "c"    ));
    cq->Af        = mxGetPr    (mxGetField(A   , 0, "f"    ));
    cq->totcompr  = mxGetPr    (mxGetField(M_fq, 0, "ct"   ));
    cq->voldiscr  = mxGetPr    (mxGetField(M_fq, 0, "vd"   ));
    cq->phasemobf = mxGetPr    (mxGetField(M_fq, 0, "pmobf"));
}


/* ---------------------------------------------------------------------- */
static void
set_effective_well_params(well_t *W, const mxArray *M_fq,
                          struct completion_data *wconn)
/* ---------------------------------------------------------------------- */
{
    const mxArray *A;

    A               =         mxGetField(M_fq, 0, "A"      ) ;
    wconn->A        = mxGetPr(mxGetField(A   , 0, "p"      ));
    wconn->gpot     = mxGetPr(mxGetField(M_fq, 0, "wdp"    ));
    wconn->phasemob = mxGetPr(mxGetField(M_fq, 0, "perfmob"));
}



/*
 * [x, Tf, Gf, wbhp, wflux] = mex_cfs_tpfa(x, G, T, ft, fq,
 *                                         Zf, dt, p0, pvol,
 *                                         W, bc, src)
 *
 * fq.A.c   = (R/B)_c
 * fq.A.f   = (R/B)_f
 * fq.ct    = total compressibility (per cell)
 * fq.pmobf = phase mobility per face
 */
/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    size_t     nc, nconn_tot;
    grid_t    *g;
    double    *src, *cpress, *fpress, *fflux;
    double     dt, *trans, *gravcap_f, *cpress0, *porevol;

    int        fld_no;

    char errmsg[1023 + 1];

    flowbc_t                *bc;
    struct mrst_well        *W;
    well_t                  *wdesc;
    well_control_t          *wctrl;

    double *WI, *wdp, *wbhp, *wflux;

    struct compr_quantities  cq;
    struct completion_data   wconn;
    struct cfs_tpfa_data    *h;
    mxArray                 *x_fpress;

    if (args_ok(nlhs, nrhs, prhs)) {
        set_mrst_cq(prhs[4], &cq);

        g   = mrst_grid(prhs[1]);
        W   = NULL;
        bc  = NULL;
        src = NULL;

        if (g != NULL) {
            W   = mrst_well  (                    prhs[9]);
            bc  = mrst_flowbc(g->number_of_faces, prhs[10]);
            src = mrst_src   (g->number_of_cells, prhs[11]);
        }

        if (W == NULL) {
            wdesc = NULL;  wctrl = NULL;  WI = NULL;  wdp = NULL;

            plhs[3] = mxCreateDoubleMatrix(0, 1, mxREAL);
            plhs[4] = mxCreateDoubleMatrix(0, 1, mxREAL);
        } else {
            wdesc    = W->wdesc;
            wctrl    = W->wctrl;
            wconn.WI = W->wdata->WI;

            set_effective_well_params(wdesc, prhs[4], &wconn);

            plhs[3] = mxCreateDoubleMatrix(wdesc->number_of_wells, 1, mxREAL);
            plhs[4] = mxCreateDoubleMatrix(mrst_well_count_totperf(prhs[9]), 1, mxREAL);
        }

        if ((g != NULL) && (bc != NULL) && (src != NULL)) {
            nc        = g->number_of_cells;
            nconn_tot = g->cell_facepos[nc];

            plhs[0]   = mxDuplicateArray    (prhs[0]);
            plhs[1]   = mxCreateDoubleMatrix(cq.nphases, g->number_of_faces, mxREAL);
            plhs[2]   = mxDuplicateArray    (plhs[1]);

            if ((fld_no = mxGetFieldNumber(plhs[0], "facePressure")) < 0) {
                fld_no   = mxAddField(plhs[0], "facePressure");
                x_fpress = mxCreateDoubleMatrix(g->number_of_faces, 1, mxREAL);
                mxSetFieldByNumber(plhs[0], 0, fld_no, x_fpress);
            }

            h = cfs_tpfa_construct(g, wdesc, cq.nphases);

            if (h != NULL) {
                dt        = mxGetScalar(prhs[6]);
                trans     = mxGetPr    (prhs[3]);
                gravcap_f = mxGetPr    (prhs[5]);
                cpress0   = mxGetPr    (prhs[7]);
                porevol   = mxGetPr    (prhs[8]);

                cfs_tpfa_assemble(g, dt, wdesc, bc, src, &cq, trans,
                                  gravcap_f, wctrl, &wconn,
                                  cpress0, porevol, h);

                callMWUMFPACK(h->A->m, h->A->ia, h->A->ja,
                              h->A->sa, h->b, h->x);

                cpress = mxGetPr(mxGetField(plhs[0], 0, "pressure"));
                fpress = mxGetPr(mxGetFieldByNumber(plhs[0], 0, fld_no));
                fflux  = mxGetPr(mxGetField(plhs[0], 0, "flux"));
                wbhp   = mxGetPr(plhs[3]);
                wflux  = mxGetPr(plhs[4]);

                cfs_tpfa_press_flux(g, bc, wdesc, cq.nphases, trans,
                                    cq.phasemobf, gravcap_f, &wconn, h,
                                    cpress, fflux, wbhp, wflux);

                cfs_tpfa_fpress(g, bc, cq.nphases, mxGetPr(prhs[2]) /* htrans */,
                                cq.phasemobf, gravcap_f, cpress, fflux, fpress);

                cfs_tpfa_retrieve_masstrans(g, cq.nphases, h, mxGetPr(plhs[1]));
                cfs_tpfa_retrieve_gravtrans(g, cq.nphases, h, mxGetPr(plhs[2]));
            }

            cfs_tpfa_destroy(h);
            mrst_src_deallocate(src);
            mrst_flowbc_deallocate(bc);
            mrst_well_deallocate(W);
        }

        free_mrst_grid(g);
    } else {
        sprintf(errmsg,
                "Calling sequence is\n"
                "\t[x, Tf, Gf, wbhp, wflux] = %s(x, G, T, ft, fq, Zf,"
                " dt, p0, pvol, W, bc, src)",
                mexFunctionName());
        mexErrMsgTxt(errmsg);
    }
}
