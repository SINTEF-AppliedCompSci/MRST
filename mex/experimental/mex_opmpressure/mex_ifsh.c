#include <stddef.h>
#include <string.h>

#include <mex.h>

#include "mrst_objects.h"

#include "grid.h"
#include "mrst_api.h"
#include "call_umfpack.h"

#include "flow_bc.h"
#include "well.h"

#include "mimetic.h"

#include "sparse_sys.h"

#include "fsh.h"

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
verify_rock(const mxArray *rock)
/* ---------------------------------------------------------------------- */
{
    return mxIsStruct(rock) &&
        (mxGetFieldNumber(rock, "perm") >= 0);
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


/* [x, wbhp, wflux] = mex_ifsh(x, G, rock, W, bc, src) */
/* ---------------------------------------------------------------------- */
static int
args_ok(int nlhs, int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok;

    ok = (nlhs == 3) && (nrhs == 6);

    ok = ok && verify_state(prhs[0]);
    ok = ok && verify_grid (prhs[1]);
    ok = ok && verify_rock (prhs[2]);
    ok = ok && verify_well (prhs[3]);
    ok = ok && verify_bc   (prhs[4]);
    ok = ok && verify_src  (prhs[5]);

    return ok;
}


struct disc_data {
    int    *ncf;
    double *Binv, *gpress, *totmob, *omega;
    double *ddata;
};


/* ---------------------------------------------------------------------- */
static void
deallocate_disc_data(struct disc_data *data)
/* ---------------------------------------------------------------------- */
{
    if (data != NULL) {
        if (data->ddata != NULL) { mxFree(data->ddata); }
        if (data->ncf   != NULL) { mxFree(data->ncf);   }

        mxFree(data);
    }
}


/* ---------------------------------------------------------------------- */
static struct disc_data *
allocate_disc_data(grid_t *g, struct fsh_data *h)
/* ---------------------------------------------------------------------- */
{
    size_t nc, ngconn, ngconn2, ddata_sz;
    struct disc_data *new;

    new = mxMalloc(1 * sizeof *new);

    if (new != NULL) {
        nc      = g->number_of_cells;
        ngconn  = g->cell_facepos[nc];
        ngconn2 = h->sum_ngconn2;

        ddata_sz  = ngconn2;    /* Binv */
        ddata_sz += ngconn;     /* gpress */
        ddata_sz += 2 * nc;     /* totmob + omega */

        new->ncf   = mxMalloc(ngconn   * sizeof *new->ncf);
        new->ddata = mxMalloc(ddata_sz * sizeof *new->ddata);

        if ((new->ncf == NULL) || (new->ddata == NULL)) {
            deallocate_disc_data(new);
            new = NULL;
        } else {
            new->Binv   = new->ddata;
            new->gpress = new->Binv   + ngconn2;
            new->totmob = new->gpress + ngconn;
            new->omega  = new->totmob + nc;
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
static void
count_gconn(grid_t *g, int *ncf)
/* ---------------------------------------------------------------------- */
{
    size_t i, nc;

    nc = g->number_of_cells;

    for (i = 0; i < nc; i++) {
        ncf[i] = g->cell_facepos[i + 1] - g->cell_facepos[i];
    }
}


/* ---------------------------------------------------------------------- */
static void
compute_gpress(grid_t *g, double *grav, double *gpress)
/* ---------------------------------------------------------------------- */
{
    size_t d, i, j, c, nc;

    double *cc, *fc;

    nc = g->number_of_cells;
    d  = g->dimensions;
    i  = 0;

    for (c = 0; c < nc; c++) {
        cc = g->cell_centroids + (c * d);

        for (; i < g->cell_facepos[c + 1]; i++) {
            fc = g->face_centroids + (g->cell_faces[i] * d);

            gpress[i] = 0.0;
            for (j = 0; j < d; j++) {
                gpress[i] += grav[j] * (fc[j] - cc[j]);
            }
        }
    }
}


/* ---------------------------------------------------------------------- */
static void
set_mobility(grid_t *g, double *totmob)
/* ---------------------------------------------------------------------- */
{
    size_t c, nc;

    nc = g->number_of_cells;

    for (c = 0; c < nc; c++) { totmob[c] = 1.0; }
}


/* ---------------------------------------------------------------------- */
static void
set_omega(grid_t *g, double *omega)
/* ---------------------------------------------------------------------- */
{
    size_t c, nc;

    nc = g->number_of_cells;

    for (c = 0; c < nc; c++) { omega[c] = 0.0; }
}


/*
 * [x, wbhp, wflux] = mex_ifsh(x, G, rock, W, bc, src)
 */
/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    double    grav[3] = { 0.0 };

    grid_t    *g;
    flowbc_t  *bc;
    double    *src, *perm;

    double *cpress, *fflux;

    struct mrst_well *W;
    well_t           *wdesc;
    well_control_t   *wctrl;

    double *WI, *wdp, *wbhp, *wflux;

    char errmsg[1023 + 1];

    struct fsh_data  *h;
    struct disc_data *disc_data;

    if (args_ok(nlhs, nrhs, prhs)) {
        plhs[0] = mxDuplicateArray(prhs[0]);

        g  = mrst_grid(prhs[1]);
        bc = NULL;  W = NULL;  src = NULL;  perm = NULL;

        if (g != NULL) {
            W    = mrst_well  (                    prhs[3]);
            bc   = mrst_flowbc(g->number_of_faces, prhs[4]);
            src  = mrst_src   (g->number_of_cells, prhs[5]);
            perm = mrst_perm  (g->dimensions,      prhs[2]);
        }

        if (W == NULL) {
            wdesc = NULL;  wctrl = NULL;  WI = NULL;  wdp = NULL;

            plhs[1] = mxCreateDoubleMatrix(0, 1, mxREAL);
            plhs[2] = mxCreateDoubleMatrix(0, 1, mxREAL);
        } else {
            wdesc = W->wdesc;
            wctrl = W->wctrl;
            WI    = W->wdata->WI;
            wdp   = W->wdata->wdp;

            plhs[1] = mxCreateDoubleMatrix(W->wdesc->number_of_wells, 1, mxREAL);
            plhs[2] = mxCreateDoubleMatrix(mrst_well_count_totperf(prhs[3]), 1, mxREAL);
        }

        if ((g != NULL) && (bc != NULL) &&
            (src != NULL) && (perm != NULL)) {
            h = ifsh_construct(g, wdesc);
            disc_data = NULL;

            if (h != NULL) {
                disc_data = allocate_disc_data(g, h);
            }

            if (disc_data != NULL) {
                count_gconn   (g, disc_data->ncf);
                compute_gpress(g, grav, disc_data->gpress);
                set_mobility  (g, disc_data->totmob);
                set_omega     (g, disc_data->omega);

                mim_ip_simple_all(g->number_of_cells, g->dimensions,
                                  h->max_ngconn, g->cell_facepos, g->cell_faces,
                                  g->face_cells, g->face_centroids,
                                  g->face_normals, g->face_areas,
                                  g->cell_centroids, g->cell_volumes, perm,
                                  disc_data->Binv);

                csrmatrix_zero(h->A);
                vector_zero(h->A->m, h->b);

                ifsh_assemble(bc, src,
                              disc_data->Binv, disc_data->gpress,
                              wctrl, WI, wdp,
                              h);

                callMWUMFPACK(h->A->m, h->A->ia, h->A->ja, h->A->sa, h->b, h->x);

                cpress = mxGetPr(mxGetField(plhs[0], 0, "pressure"));
                fflux  = mxGetPr(mxGetField(plhs[0], 0, "flux"));

                wbhp   = mxGetPr(plhs[1]);
                wflux  = mxGetPr(plhs[2]);

                fsh_press_flux(g, disc_data->Binv, disc_data->gpress,
                               h, cpress, fflux, wbhp, wflux);
            }

            deallocate_disc_data(disc_data);
            fsh_destroy(h);

            mrst_flowbc_deallocate(bc);
            mrst_perm_deallocate  (perm);
            mrst_src_deallocate   (src);
        }

        mrst_well_deallocate(W);
        free_mrst_grid(g);
    } else {
        sprintf(errmsg,
                "Calling sequence is\n"
                "\t[x, wbhp, wflux] = %s(x, G, rock, W, bc, src)",
                mexFunctionName());
        mexErrMsgTxt(errmsg);
    }
}
