#include <assert.h>
#include <string.h>
#include <mex.h>

#include "grid.h"
#include "mrst_api.h"
#include "mrst_objects.h"
#include "mimetic.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))


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


/* ------------------------------------------------------------------ */
static int
verify_grid_structure(const mxArray *G)
/* ------------------------------------------------------------------ */
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


/* ------------------------------------------------------------------ */
static int
args_ok(int nlhs, int nrhs, const mxArray *prhs[])
/* ------------------------------------------------------------------ */
{
    int grid_ok = 0, rock_ok = 0;

    if (nrhs > 0) {
        grid_ok = verify_grid_structure(prhs[0]);
    }

    if (nrhs > 1) {
        /* rock must contain a field 'perm' */
        rock_ok = mxIsStruct(prhs[1]) &&
                  (mxGetFieldNumber(prhs[1], "perm") >= 0);
    }

    return (nlhs == 1) && grid_ok && rock_ok;
}


/* ---------------------------------------------------------------------- */
static void
count_connections(grid_t *g, int *max_nconn, int *sum_nconn2)
/* ---------------------------------------------------------------------- */
{
    int c, nconn;

    *max_nconn = *sum_nconn2 = 0;

    for (c = 0; c < g->number_of_cells; c++) {
        nconn = g->cell_facepos[c + 1] - g->cell_facepos[c];

        *max_nconn   = MAX(*max_nconn, nconn);
        *sum_nconn2 += nconn * nconn;
    }
}


/*
 * BI = mex_ip_simple(G, rock)
 */

/* ------------------------------------------------------------------ */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ------------------------------------------------------------------ */
{
    int     max_nconn, sum_nconn2;
    grid_t *g;
    double *perm;
    char    errmsg[1023 + 1];

    double *Binv;

    if (args_ok(nlhs, nrhs, prhs)) {
        g          = mrst_grid(prhs[0]);
        perm       = NULL;
        sum_nconn2 = 0;

        if (g != NULL) {
            perm = mrst_perm (g->dimensions, prhs[1]);
            count_connections(g, &max_nconn, &sum_nconn2);
        }
        
        /* Create return values */
        plhs[0] = mxCreateDoubleMatrix(sum_nconn2, 1, mxREAL);

        /* Compute IP for all cells (reservoir) */
        Binv    = mxGetPr(plhs[0]); /* mxCreateDoubleMatrix == ZEROS */

        if (perm != NULL) {
            /* mim_ip_simple_all() requires zeroed target array... */
            mim_ip_simple_all(g->number_of_cells, g->dimensions,
                              max_nconn, g->cell_facepos, g->cell_faces,
                              g->face_cells, g->face_centroids,
                              g->face_normals, g->face_areas,
                              g->cell_centroids, g->cell_volumes,
                              perm, Binv);
        }

        mxFree(perm);
        free_mrst_grid(g);
    } else {
        sprintf(errmsg,
                "Calling sequence is\n"
                "\tBI = %s(G, rock)\n",
                mexFunctionName());

        mexErrMsgTxt(errmsg);
    }
}
