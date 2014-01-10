#include <stddef.h>
#include <string.h>

#include <mex.h>

#include "flow_bc.h"
#include "grid.h"
#include "mrst_api.h"
#include "ifsh_ms.h"
#include "sparse_sys.h"
#include "call_umfpack.h"


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


/* x = mex_ifsh(x, G, rock, p, src) */
/* ---------------------------------------------------------------------- */
static int
args_ok(int nlhs, int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok;

    ok = (nlhs == 1) && (nrhs == 5);

    ok = ok && verify_state(prhs[0]);
    ok = ok && verify_grid (prhs[1]);
    ok = ok && verify_rock (prhs[2]);
    ok = ok && (mxIsDouble (prhs[3]) || mxIsInt32(prhs[3]));
    ok = ok && verify_src  (prhs[4]);

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
static double *
mrst_perm(int d, const mxArray *rock)
/* ---------------------------------------------------------------------- */
{
    return getPermeability(mxGetField(rock, 0, "perm"), d);
}


/* ---------------------------------------------------------------------- */
static int *
mrst_partition(const mxArray *M_p)
/* ---------------------------------------------------------------------- */
{
    size_t i, n;
    int    *p;

    n = mxGetNumberOfElements(M_p);

    p = mxMalloc(n * sizeof *p);

    if (p != NULL) {
        assign_int_vector(M_p, p);

        for (i = 0; i < n; i++) { p[i] -= 1; }
    }

    return p;
}


struct disc_data {
    double *totmob;
};


/* ---------------------------------------------------------------------- */
static void
deallocate_disc_data(struct disc_data *data)
/* ---------------------------------------------------------------------- */
{
    if (data != NULL) {
        if (data->totmob != NULL) { mxFree(data->totmob); }

        mxFree(data);
    }
}


/* ---------------------------------------------------------------------- */
static struct disc_data *
allocate_disc_data(grid_t *g)
/* ---------------------------------------------------------------------- */
{
    size_t nc;
    struct disc_data *new;

    new = mxMalloc(1 * sizeof *new);

    if (new != NULL) {
        nc = g->number_of_cells;
        new->totmob = mxMalloc(nc * sizeof *new->totmob);

	if (new->totmob == NULL) {
            deallocate_disc_data(new);
            new = NULL;
	}
    }

    return new;
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
local_solver(struct CSRMatrix *A, double *b, double *x)
/* ---------------------------------------------------------------------- */
{
    callMWUMFPACK(A->m, A->ia, A->ja, A->sa, b, x);
}


/*
 * x = mex_ifsh_ms(x, G, rock, p, src)
 */
/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    grid_t    *g;
    int       *p;
    double    *src, *perm;

    double *cpress, *fflux;

    struct disc_data    *disc;
    struct ifsh_ms_data *h;

    if (args_ok(nlhs, nrhs, prhs)) {
        plhs[0] = mxDuplicateArray(prhs[0]);

        g    = mrst_grid(prhs[1]);
        perm = NULL;    p = NULL;    src = NULL;

        if (g != NULL) {
            perm = mrst_perm  (g->dimensions,      prhs[2]);
            p    = mrst_partition(                 prhs[3]);
            src  = mrst_src   (g->number_of_cells, prhs[4]);
            disc = allocate_disc_data(g);
        }

        if ((g != NULL) && (src != NULL) && (p != NULL) && (perm != NULL)) {
            set_mobility(g, disc->totmob);

            h = ifsh_ms_construct(g, p, perm, src,
                                  disc->totmob, local_solver);

            if (h != NULL) {
                ifsh_ms_assemble(src, disc->totmob, h);

                callMWUMFPACK(h->A->m, h->A->ia, h->A->ja, h->A->sa, h->b, h->x);

                cpress = mxGetPr(mxGetField(plhs[0], 0, "pressure"));
                fflux  = mxGetPr(mxGetField(plhs[0], 0, "flux"));

                ifsh_ms_press_flux(g, h, cpress, fflux);
            }

            ifsh_ms_destroy(h);
            mxFree(src);  mxFree(p);  mxFree(perm);
        }

        deallocate_disc_data(disc);
        free_mrst_grid(g);
    }
}
