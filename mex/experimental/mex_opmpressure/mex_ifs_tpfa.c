#include <stddef.h>
#include <string.h>

#include <mex.h>

#include "mrst_api.h"
#include "call_umfpack.h"

#include "grid.h"
#include "flow_bc.h"
#include "well.h"

#include "sparse_sys.h"

#include "trans_tpfa.h"
#include "ifs_tpfa.h"

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


/* x = mex_ifs_tpfa(x, G, rock, src) */
/* ---------------------------------------------------------------------- */
static int
args_ok(int nlhs, int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok;

    ok = (nlhs == 1) && (nrhs == 4);

    ok = ok && verify_state(prhs[0]);
    ok = ok && verify_grid (prhs[1]);
    ok = ok && verify_rock (prhs[2]);
    ok = ok && verify_src  (prhs[3]);

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


/* ---------------------------------------------------------------------- */
static double *
mrst_perm(int d, const mxArray *rock)
/* ---------------------------------------------------------------------- */
{
    return getPermeability(mxGetField(rock, 0, "perm"), d);
}


struct disc_data {
    double *htrans, *trans, *gpress, *totmob;
    double *ddata;
};


/* ---------------------------------------------------------------------- */
static void
deallocate_disc_data(struct disc_data *data)
/* ---------------------------------------------------------------------- */
{
    if (data != NULL) {
        if (data->ddata != NULL) { mxFree(data->ddata); }

        mxFree(data);
    }
}


/* ---------------------------------------------------------------------- */
static struct disc_data *
allocate_disc_data(grid_t *g)
/* ---------------------------------------------------------------------- */
{
    size_t nc, nf, ngconn, ddata_sz;
    struct disc_data *new;

    new = mxMalloc(1 * sizeof *new);

    if (new != NULL) {
        nc      = g->number_of_cells;
        nf      = g->number_of_faces;
        ngconn  = g->cell_facepos[nc];

        ddata_sz  = ngconn;     /* htrans */
        ddata_sz += nf;         /* trans */
        ddata_sz += ngconn;     /* gpress */
        ddata_sz += nc;         /* totmob */

        new->ddata = mxMalloc(ddata_sz * sizeof *new->ddata);

        if (new->ddata == NULL) {
            deallocate_disc_data(new);
            new = NULL;
        } else {
            new->htrans = new->ddata;
            new->trans  = new->htrans + ngconn;
            new->gpress = new->trans  + nf;
            new->totmob = new->gpress + ngconn;
        }
    }

    return new;
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


/*
 * x = mex_ifs_tpfa(x, G, rock, src)
 */
/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    grid_t    *g;
    double    *src, *perm;

    double *cpress, *fflux;

    double grav[] = { 0.0, 0.0, 10.0 };

    char errmsg[1023 + 1];

    struct ifs_tpfa_data *h;
    struct disc_data *disc_data;

    if (args_ok(nlhs, nrhs, prhs)) {
        plhs[0] = mxDuplicateArray(prhs[0]);

        g  = mrst_grid(prhs[1]);
        src = NULL;  perm = NULL;

        if (g != NULL) {
            src  = mrst_src   (g->number_of_cells, prhs[3]);
            perm = mrst_perm  (g->dimensions,      prhs[2]);
        }

        if ((g != NULL) && (src != NULL) && (perm != NULL)) {
            h = ifs_tpfa_construct(g);
            disc_data = allocate_disc_data(g);

            if (disc_data != NULL) {
                set_mobility  (g, disc_data->totmob);

                compute_gpress(g, grav, disc_data->gpress);

                tpfa_htrans_compute(g, perm, disc_data->htrans);
                tpfa_eff_trans_compute(g, disc_data->totmob,
                                       disc_data->htrans,
                                       disc_data->trans);

                ifs_tpfa_assemble(g, disc_data->trans, src, disc_data->gpress,
                                  h);

                callMWUMFPACK(h->A->m, h->A->ia, h->A->ja, h->A->sa, h->b, h->x);

                cpress = mxGetPr(mxGetField(plhs[0], 0, "pressure"));
                fflux  = mxGetPr(mxGetField(plhs[0], 0, "flux"));

                ifs_tpfa_press_flux(g, disc_data->trans, h, cpress, fflux);
            }

            deallocate_disc_data(disc_data);
            ifs_tpfa_destroy(h);
            mxFree(perm);    mxFree(src);
        }

        free_mrst_grid(g);
    } else {
        sprintf(errmsg,
                "Calling sequence is\n"
                "\tx = %s(x, G, rock, src)",
                mexFunctionName());
        mexErrMsgTxt(errmsg);
    }
}
