#include <stddef.h>
#include <string.h>

#include <mex.h>

#include "mrst_api.h"
#include "call_umfpack.h"

#include "grid.h"
#include "flow_bc.h"
#include "well.h"

#include "compr_quant.h"
#include "mimetic.h"

#include "sparse_sys.h"

#include "fsh.h"

/* ---------------------------------------------------------------------- */
static int
verify_state(const mxArray *state)
/* ---------------------------------------------------------------------- */
{
    return !mxIsEmpty(state) && mxIsStruct(state) &&
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
verify_fluid_quant(const mxArray *fq)
/* ---------------------------------------------------------------------- */
{
    int ok;

    ok =        !mxIsEmpty(fq) && mxIsStruct(fq);
    ok = ok && (mxGetFieldNumber(fq, "totmob")   >= 0);
    ok = ok && (mxGetFieldNumber(fq, "totcompr") >= 0);
    ok = ok && (mxGetFieldNumber(fq, "zeta")     >= 0);
    ok = ok && (mxGetFieldNumber(fq, "porevol")  >= 0); /* quirk */

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

#if 0
    /* Discount well for the time being */
    if (!ok) {
        ok =        mxIsStruct(W);
        ok = ok && (mxGetFieldNumber(W, "cells") >= 0);
        ok = ok && (mxGetFieldNumber(W, "type" ) >= 0);
        ok = ok && (mxGetFieldNumber(W, "val"  ) >= 0);
        ok = ok && (mxGetFieldNumber(W, "WI"   ) >= 0);
        ok = ok && (mxGetFieldNumber(W, "dZ"   ) >= 0);
    }
#endif

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
 * [x, wbhp, wflux] = mex_fsh(x, G, rock, ...
 *                            fq, p0, dt, ...
 *                            W, bc, src)
 */
/* ---------------------------------------------------------------------- */
static int
args_ok(int nlhs, int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok;

    ok = (nlhs == 3) && (nrhs == 9);

    ok = ok && verify_state      (prhs[0]); /* x */
    ok = ok && verify_grid       (prhs[1]); /* G */
    ok = ok && verify_rock       (prhs[2]); /* rock */
    ok = ok && verify_fluid_quant(prhs[3]); /* fq */
    ok = ok && mxIsDouble        (prhs[4]); /* p0 */
    ok = ok && verify_timestep   (prhs[5]); /* dt */
    ok = ok && verify_well       (prhs[6]); /* W */
    ok = ok && verify_bc         (prhs[7]); /* bc */
    ok = ok && verify_src        (prhs[8]); /* src */

    return ok;
}


/* ---------------------------------------------------------------------- */
static int *
get_int_vector(const mxArray *v)
/* ---------------------------------------------------------------------- */
{
    size_t i, n;

    int *ret, *pi;
    double *pd;

    n = mxGetNumberOfElements(v);

    ret = mxMalloc(n * sizeof *ret);

    if (ret != NULL) {
        if (mxIsDouble(v)) {
            pd = mxGetPr(v);

            for (i = 0; i < n; i++) { ret[i] = pd[i]; }
        } else {
            mxAssert (mxIsInt32(v),
                      "Only DOUBLE (flint) and INT32 types "
                      "supported for indices");

            pi = mxGetData(v);

            memcpy(ret, pi, n * sizeof *ret);
        }
    }

    return ret;
}


/* ---------------------------------------------------------------------- */
static void
free_mrst_grid(grid_t *g)
/* ---------------------------------------------------------------------- */
{
    if (g != NULL) {
        if (g->cell_centroids != NULL) { mxFree(g->cell_centroids); }
        if (g->face_normals   != NULL) { mxFree(g->face_normals);   }
        if (g->face_centroids != NULL) { mxFree(g->face_centroids); }

        if (g->cell_facepos   != NULL) { mxFree(g->cell_facepos);   }
        if (g->cell_faces     != NULL) { mxFree(g->cell_faces);     }
        if (g->face_cells     != NULL) { mxFree(g->face_cells);     }

        mxFree(g);
    }
}


/* ---------------------------------------------------------------------- */
static grid_t *
mrst_grid(const mxArray *G)
/* ---------------------------------------------------------------------- */
{
    int     copy_ok;
    grid_t *g;

    g = mxMalloc(1 * sizeof *g);

    if (g != NULL) {
        copy_ok  = (g->face_cells     = getFaceCellNeighbors(G)) != NULL;
        copy_ok += (g->cell_faces     = getCellFaces        (G)) != NULL;
        copy_ok += (g->cell_facepos   = getCellFacePos      (G)) != NULL;

        copy_ok += (g->face_centroids = getFaceCentroids    (G)) != NULL;
        copy_ok += (g->face_normals   = getFaceNormals      (G)) != NULL;
        copy_ok += (g->cell_centroids = getCellCentroids    (G)) != NULL;

        if (copy_ok != 6) {
            free_mrst_grid(g);
            g = NULL;
        } else {
            g->dimensions      = getNumberOfDimensions(G);
            g->number_of_cells = getNumberOfCells     (G);
            g->number_of_faces = getNumberOfFaces     (G);

            g->cell_volumes    = getCellVolumes(G);
            g->face_areas      = getFaceAreas  (G);
        }
    }

    return g;
}


/* ---------------------------------------------------------------------- */
static flowbc_t *
mrst_flowbc(size_t nf, const mxArray *BC)
/* ---------------------------------------------------------------------- */
{
    char type_str[] = "pressure";
    flowbc_t *bc;

    int    *bcf;
    double *bcv;

    mxArray *bc_type, *bct;
    mwSize i, n;

    bc = allocate_flowbc(nf);

    if (bc != NULL) {
        if (!mxIsEmpty(BC)) {
            n = mxGetNumberOfElements(mxGetField(BC, 0, "face"));

            bcf = get_int_vector(mxGetField(BC, 0, "face"));
            bct = mxGetField(BC, 0, "type");
            bcv = mxGetPr(mxGetField(BC, 0, "value"));

            if (bcf != NULL) {
                for (i = 0; i < n; i++) {
                    bc_type = mxGetCell(bct, i);

                    mxAssert (mxIsChar(bc_type),
                              "bc.type must be cell array of strings.");

                    mxGetString(bc_type, type_str, sizeof type_str);

                    if      (strcmp(type_str, "pressure") == 0)
                    {
                        bc->type [bcf[i] - 1] = PRESSURE;
                        bc->bcval[bcf[i] - 1] = bcv[i];
                    }
                    else if (strcmp(type_str, "flux"    ) == 0)
                    {
                        bc->type [bcf[i] - 1] = FLUX;
                        bc->bcval[bcf[i] - 1] = bcv[i];
                    }
                }

                mxFree(bcf);
            }
        }
    }

    return bc;
}


/* ---------------------------------------------------------------------- */
static double *
mrst_src(size_t nc, const mxArray *SRC)
/* ---------------------------------------------------------------------- */
{
    size_t i, n;
    int    *cix;
    double *crat, *src;

    src = mxMalloc(nc * sizeof *src);

    if (src != NULL) {
        for (i = 0; i < nc; i++) { src[i] = 0.0; }

        if (!mxIsEmpty(SRC)) {
            cix  = get_int_vector(mxGetField(SRC, 0, "cell"));
            crat = mxGetPr(mxGetField(SRC, 0, "rate"));

            if (cix != NULL) {
                n = mxGetNumberOfElements(mxGetField(SRC, 0, "cell"));

                for (i = 0; i < n; i++) {
                    src[cix[i] - 1] = crat[i];
                }

                mxFree(cix);
            }
        }
    }

    return src;
}


struct well_data {
    double *WI, *wdp;
};

struct mrst_well {
    well_t           *wdesc;
    well_control_t   *wctrl;
    struct well_data *wdata;
};


/* ---------------------------------------------------------------------- */
static void
well_descriptor_deallocate(well_t *wdesc)
/* ---------------------------------------------------------------------- */
{
    if (wdesc != NULL) {
        if (wdesc->well_cells   != NULL) { mxFree(wdesc->well_cells)  ; }
        if (wdesc->well_connpos != NULL) { mxFree(wdesc->well_connpos); }

        mxFree(wdesc);
    }
}


/* ---------------------------------------------------------------------- */
static well_t *
well_descriptor_allocate(size_t nw, size_t nperf)
/* ---------------------------------------------------------------------- */
{
    well_t *new;

    new = mxMalloc(1 * sizeof *new);

    if (new != NULL) {
        new->well_connpos = mxCalloc(nw + 1, sizeof *new->well_connpos);
        new->well_cells   = mxMalloc(nperf * sizeof *new->well_cells);

        if ((new->well_connpos == NULL) || (new->well_cells == NULL)) {
            well_descriptor_deallocate(new);
            new = NULL;
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
static void
well_control_deallocate(well_control_t *wctrl)
/* ---------------------------------------------------------------------- */
{
    if (wctrl != NULL) {
        if (wctrl->target != NULL) { mxFree(wctrl->target); }
        if (wctrl->ctrl   != NULL) { mxFree(wctrl->ctrl)  ; }
        if (wctrl->type   != NULL) { mxFree(wctrl->type)  ; }

        mxFree(wctrl);
    }
}


/* ---------------------------------------------------------------------- */
static well_control_t *
well_control_allocate(size_t nw)
/* ---------------------------------------------------------------------- */
{
    well_control_t *new;

    new = mxMalloc(1 * sizeof *new);

    if (new != NULL) {
        new->type   = mxMalloc(nw * sizeof *new->type);
        new->ctrl   = mxMalloc(nw * sizeof *new->ctrl);
        new->target = mxMalloc(nw * sizeof *new->target);

        if ((new->type == NULL) || (new->ctrl == NULL) ||
            (new->target == NULL)) {
            well_control_deallocate(new);
            new = NULL;
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
static void
well_data_deallocate(struct well_data *wdata)
/* ---------------------------------------------------------------------- */
{
    if (wdata != NULL) {
        if (wdata->wdp != NULL) { mxFree(wdata->wdp); }
        if (wdata->WI  != NULL) { mxFree(wdata->WI) ; }

        mxFree(wdata);
    }
}


/* ---------------------------------------------------------------------- */
static struct well_data *
well_data_allocate(size_t nperf)
/* ---------------------------------------------------------------------- */
{
    struct well_data *new;

    new = mxMalloc(1 * sizeof *new);

    if (new != NULL) {
        new->WI  = mxMalloc(nperf * sizeof *new->WI);
        new->wdp = mxMalloc(nperf * sizeof *new->wdp);

        if ((new->WI == NULL) || (new->wdp == NULL)) {
            well_data_deallocate(new);
            new = NULL;
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
static void
mrst_well_deallocate(struct mrst_well *W)
/* ---------------------------------------------------------------------- */
{
    if (W != NULL) {
        well_data_deallocate      (W->wdata);
        well_control_deallocate   (W->wctrl);
        well_descriptor_deallocate(W->wdesc);

        mxFree(W);
    }
}


/* ---------------------------------------------------------------------- */
static struct mrst_well *
mrst_well_allocate(size_t nw, size_t nperf)
/* ---------------------------------------------------------------------- */
{
    struct mrst_well *new;

    new = mxMalloc(1 * sizeof *new);

    if (new != NULL) {
        new->wdesc = well_descriptor_allocate(nw, nperf);
        new->wctrl = well_control_allocate   (nw);
        new->wdata = well_data_allocate      (nperf);

        if ((new->wdesc == NULL) || (new->wctrl == NULL) ||
            (new->wdata == NULL)) {
            mrst_well_deallocate(new);
            new = NULL;
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
static size_t
mrst_well_count_totperf(const mxArray *W)
/* ---------------------------------------------------------------------- */
{
    mwIndex fld_no;
    size_t nw, w, totperf;

    nw = mxGetNumberOfElements(W);

    fld_no = mxGetFieldNumber(W, "cells");

    totperf = 0;
    for (w = 0; w < nw; w++) {
        totperf += mxGetNumberOfElements(mxGetFieldByNumber(W, w, fld_no));
    }

    return totperf;
}


/* ---------------------------------------------------------------------- */
static void
assign_int_vector(const mxArray *M_v, int *v)
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

        memcpy(v, pi, n * sizeof *v);
    }
}


/* ---------------------------------------------------------------------- */
static void
assign_double_vector(const mxArray *M_v, double *v)
/* ---------------------------------------------------------------------- */
{
    size_t i, n;

    int    *pi;
    double *pd;

    n = mxGetNumberOfElements(M_v);

    if (mxIsDouble(M_v)) {
        pd = mxGetPr(M_v);

        memcpy(v, pd, n * sizeof *v);
    } else {
        pi = mxGetData(M_v);

        for (i = 0; i < n; i++) { v[i] = pi[i]; }
    }
}


/* ---------------------------------------------------------------------- */
static void
mrst_well_set_descriptor(const mxArray *W, well_t *wdesc)
/* ---------------------------------------------------------------------- */
{
    int      i;
    mwIndex  fld;
    size_t   w, nw;
    mxArray *cells;

    nw = mxGetNumberOfElements(W);
    fld = mxGetFieldNumber(W, "cells");

    for (w = 0; w < nw; w++) {
        cells = mxGetFieldByNumber(W, w, fld);

        wdesc->well_connpos[w + 1] = wdesc->well_connpos[w] +
            (int) mxGetNumberOfElements(cells);

        assign_int_vector(cells,
                          wdesc->well_cells +
                          wdesc->well_connpos[w]);
    }

    for (i = 0; i < wdesc->well_connpos[nw]; i++) {
        wdesc->well_cells[i] -= 1; /* 1-based indexing in M */
    }

    wdesc->number_of_wells = (int) nw;
}


/* ---------------------------------------------------------------------- */
static void
mrst_well_set_control(const mxArray *W, well_control_t *wctrl)
/* ---------------------------------------------------------------------- */
{
    size_t w, nw;
    char   ctrl_type[] = "rate";

    mxArray *target,  *type;
    mwIndex  trgt_fld, typ_fld;

    nw = mxGetNumberOfElements(W);

    typ_fld  = mxGetFieldNumber(W, "type");
    trgt_fld = mxGetFieldNumber(W, "val" );

    for (w = 0; w < nw; w++) {
        type   = mxGetFieldByNumber(W, w, typ_fld);
        target = mxGetFieldByNumber(W, w, trgt_fld);

        mxAssert (mxIsChar(type), "'W.type' field must be CHAR string.");
        mxGetString(type, ctrl_type, sizeof ctrl_type);

        if (strcmp(ctrl_type, "bhp") == 0) {
            wctrl->ctrl[w] = BHP;
        } else {
            mxAssert (strcmp(ctrl_type, "rate") == 0,
                      "'W.type' must be 'bhp' or 'rate'.");
            wctrl->ctrl[w] = RATE;
        }

        mxAssert (mxGetNumberOfElements(target),
                  "'W.val' must be a scalar value.");
        wctrl->target[w] = mxGetPr(target)[0];
    }
}


/* ---------------------------------------------------------------------- */
static void
mrst_well_set_data(const mxArray *W, struct well_data *wdata)
/* ---------------------------------------------------------------------- */
{
    mwIndex WI_fld, wdp_fld;
    size_t  w, nw, p;

    mxArray *WI, *wdp;

    nw = mxGetNumberOfElements(W);

    WI_fld  = mxGetFieldNumber(W, "WI");
    wdp_fld = mxGetFieldNumber(W, "dZ");

    p = 0;
    for (w = 0; w < nw; w++) {
        WI  = mxGetFieldByNumber(W, w, WI_fld);
        wdp = mxGetFieldByNumber(W, w, wdp_fld);

        assign_double_vector(WI , wdata->WI  + p);
        assign_double_vector(wdp, wdata->wdp + p);

        p += mxGetNumberOfElements(WI);
    }
}


/* ---------------------------------------------------------------------- */
static struct mrst_well *
mrst_well(const mxArray *W)
/* ---------------------------------------------------------------------- */
{
    size_t nw, nperf;
    struct mrst_well *new;

    if (!mxIsEmpty(W)) {
        nw = mxGetNumberOfElements(W);
        nperf = mrst_well_count_totperf(W);

        new = mrst_well_allocate(nw, nperf);

        if (new != NULL) {
            mrst_well_set_descriptor(W, new->wdesc);
            mrst_well_set_control   (W, new->wctrl);
            mrst_well_set_data      (W, new->wdata);
        }
    } else {
        new = NULL;
    }

    return new;
}


/* ---------------------------------------------------------------------- */
static double *
mrst_perm(int d, const mxArray *rock)
/* ---------------------------------------------------------------------- */
{
    return getPermeability(mxGetField(rock, 0, "perm"), d);
}


struct disc_data {
    int    *ncf;
    double *Binv, *Biv, *gpress, *P;
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
        ddata_sz += ngconn;     /* Biv */
        ddata_sz += ngconn;     /* gpress */
        ddata_sz += nc;         /* P */

        new->ncf   = mxMalloc(ngconn   * sizeof *new->ncf);
        new->ddata = mxMalloc(ddata_sz * sizeof *new->ddata);

        if ((new->ncf == NULL) || (new->ddata == NULL)) {
            deallocate_disc_data(new);
            new = NULL;
        } else {
            new->Binv   = new->ddata;
            new->Biv    = new->Binv   + ngconn2;
            new->gpress = new->Biv    + ngconn;
            new->P      = new->gpress + ngconn;
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
compute_compr_terms(grid_t           *g,
                    double            dt,
                    const double     *fflux,
                    const double     *p0,
                    const mxArray    *fq,
                    double           *src,
                    struct disc_data *d)
/* ---------------------------------------------------------------------- */
{
    const mxArray *zeta;
    const mxArray *ct;
    const mxArray *pv;

    zeta = mxGetField(fq, 0, "zeta");
    ct   = mxGetField(fq, 0, "totcompr");
    pv   = mxGetField(fq, 0, "porevol");

    compr_flux_term(g, fflux, mxGetPr(zeta), d->Biv);

    compr_accum_term(g->number_of_cells, dt,
                     mxGetPr(pv), mxGetPr(ct), d->P);

    compr_src_add_press_accum(g->number_of_cells, p0, d->P, src);
}


/* ---------------------------------------------------------------------- */
static void
set_mobility_effects(grid_t *g, const mxArray *fq, double *Binv)
/* ---------------------------------------------------------------------- */
{
    const mxArray *totmob;

    totmob = mxGetField(fq, 0, "totmob");

    mim_ip_mobility_update(g->number_of_cells, g->cell_facepos,
                           mxGetPr(totmob), Binv, Binv);
}


/*
 * [x, wbhp, wflux] = mex_ifsh(x, G, rock, fq, p0, dt, W, bc, src)
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
    double     dt;
    double    *src, *perm, *p0;

    double *cpress, *fflux;

    struct mrst_well *W;
    well_t           *wdesc;
    well_control_t   *wctrl;

    double *WI, *wdp, *wbhp, *wflux;

    struct fsh_data  *h;
    struct disc_data *disc_data;

    if (args_ok(nlhs, nrhs, prhs)) {
        plhs[0] = mxDuplicateArray(prhs[0]);

        g  = mrst_grid(prhs[1]);
        bc = NULL;  W = NULL;  src = NULL;  perm = NULL;
        p0 = mxGetPr(prhs[4]);
        dt = mxGetScalar(prhs[5]);

        if (g != NULL) {
            W    = mrst_well  (                    prhs[6]);
            bc   = mrst_flowbc(g->number_of_faces, prhs[7]);
            src  = mrst_src   (g->number_of_cells, prhs[8]);
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
            plhs[2] = mxCreateDoubleMatrix(mrst_well_count_totperf(prhs[6]), 1, mxREAL);
        }

        if ((g != NULL) && (bc != NULL) &&
            (src != NULL) && (perm != NULL)) {
            h = cfsh_construct(g, wdesc);
            disc_data = NULL;

            if (h != NULL) {
                disc_data = allocate_disc_data(g, h);
            }

            if (disc_data != NULL) {
                cpress = mxGetPr(mxGetField(plhs[0], 0, "pressure"));
                fflux  = mxGetPr(mxGetField(plhs[0], 0, "flux"));

                count_gconn   (g,       disc_data->ncf);
                compute_gpress(g, grav, disc_data->gpress);

                compute_compr_terms(g, dt, fflux, p0, prhs[3], src, disc_data);

                mim_ip_simple_all(g->number_of_cells, g->dimensions,
                                  h->max_ngconn,
                                  g->cell_facepos, g->cell_faces,
                                  g->face_cells, g->face_centroids,
                                  g->face_normals, g->face_areas,
                                  g->cell_centroids, g->cell_volumes, perm,
                                  disc_data->Binv);

                set_mobility_effects(g, prhs[3], disc_data->Binv);

                csrmatrix_zero(h->A);
                vector_zero(h->A->m, h->b);

                cfsh_assemble(bc, src, disc_data->Binv, disc_data->Biv,
                              disc_data->P, disc_data->gpress,
                              wctrl, WI, NULL, wdp,
                              h);

                callMWUMFPACK(h->A->m, h->A->ia, h->A->ja, h->A->sa, h->b, h->x);

                wbhp   = mxGetPr(plhs[1]);
                wflux  = mxGetPr(plhs[2]);

                fsh_press_flux(g, disc_data->Binv, disc_data->gpress,
                               h, cpress, fflux, wbhp, wflux);
            }

            deallocate_disc_data(disc_data);
            fsh_destroy(h);
            deallocate_flowbc(bc);    mxFree(perm);    mxFree(src);
        }

        mrst_well_deallocate(W);
        free_mrst_grid(g);
    } else {
    }
}
