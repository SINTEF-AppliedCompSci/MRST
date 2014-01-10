#include <string.h>

#include <mex.h>

#include "grid.h"
#include "mrst_api.h"

#include "mrst_objects.h"


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

        mxAssert (mxGetNumberOfElements(target) == 1,
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


/* ======================================================================
 * Public interface below separator
 * ====================================================================== */

/* ---------------------------------------------------------------------- */
/* Create a flow boundary condition structure to support a model with
 * 'nf' faces (i.e., G.faces.num) from an MRST 'bc' structure (addBC). */
/* ---------------------------------------------------------------------- */
flowbc_t *
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
/* Release flowbc_t memory resources allocated in mrst_flowbc() */
/* ---------------------------------------------------------------------- */
void
mrst_flowbc_deallocate(flowbc_t *bc)
/* ---------------------------------------------------------------------- */
{
    deallocate_flowbc(bc);
}


/* ---------------------------------------------------------------------- */
/* Create a source vector from an MRST 'src' object (addSource) */
/* ---------------------------------------------------------------------- */
double *
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


/* ---------------------------------------------------------------------- */
/* Release memory resources allocated in mrst_src() */
/* ---------------------------------------------------------------------- */
void
mrst_src_deallocate(double *src)
/* ---------------------------------------------------------------------- */
{
    if (src != NULL) { mxFree(src); }
}


/* ---------------------------------------------------------------------- */
/* Extract full permeability tensor from an MRST 'rock' object. */
/* ---------------------------------------------------------------------- */
double *
mrst_perm(int d, const mxArray *rock)
/* ---------------------------------------------------------------------- */
{
    return getPermeability(mxGetField(rock, 0, "perm"), d);
}


/* ---------------------------------------------------------------------- */
/* Release memory resources allocated in mrst_perm() */
/* ---------------------------------------------------------------------- */
void
mrst_perm_deallocate(double *perm)
/* ---------------------------------------------------------------------- */
{
    if (perm != NULL) { mxFree(perm); }
}


/* ---------------------------------------------------------------------- */
/* Count total number of perforations in all wells represented by MRST
 * well structure 'W'. */
/* ---------------------------------------------------------------------- */
size_t
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
/* Build (slight abstraction) of an OPM well management structure from
 * an MRST well structure, W, (addWell).  Return NULL in case of
 * allocation failure or if there are no wells in the model
 * (isempty(W)). */
/* ---------------------------------------------------------------------- */
struct mrst_well *
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
/* Release memory resource allocated in mrst_well() */
/* ---------------------------------------------------------------------- */
void
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
