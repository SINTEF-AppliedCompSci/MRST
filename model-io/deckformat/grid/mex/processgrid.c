/*===========================================================================
//
// File: processgrid.c
//
// Created: Fri Jun 19 08:46:53 2009
//
// Author: Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
//
// $Date: 2012-10-03 10:08:03 +0200 (Wed, 03 Oct 2012) $
//
// $Revision: 1005 $
//
//===========================================================================*/

/*
  Copyright 2009, 2010, 2011, 2012 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010, 2011, 2012 Statoil ASA.

  This file is part of The Open Reservoir Simulator Project (OpenRS).

  OpenRS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OpenRS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/

/* Copyright 2009, 2010 SINTEF ICT, Applied Mathematics. */
/* Mex gateway by Jostein R. Natvig, SINTEF ICT.   */


#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <mex.h>

#include "preprocess.h"
#include "mxgrdecl.h"


/* ---------------------------------------------------------------------- */
static mxArray *
allocate_nodes(size_t nnodes)
/* ---------------------------------------------------------------------- */
{
    size_t      nflds;
    const char *fields[] = { "num", "coords" };

    mxArray *nodes, *num, *coords;

    nflds  = sizeof(fields) / sizeof(fields[0]);

    nodes  = mxCreateStructMatrix(1, 1, nflds, fields);

    num    = mxCreateDoubleScalar(nnodes);
    coords = mxCreateDoubleMatrix(nnodes, 3, mxREAL);

    if ((nodes != NULL) && (num != NULL) && (coords != NULL)) {
        mxSetField(nodes, 0, "num"   , num);
        mxSetField(nodes, 0, "coords", coords);
    } else {
        if (coords != NULL) { mxDestroyArray(coords); }
        if (num    != NULL) { mxDestroyArray(num);    }
        if (nodes  != NULL) { mxDestroyArray(nodes);  }

        nodes = NULL;
    }

    return nodes;
}


/* ---------------------------------------------------------------------- */
static void
fill_nodes(mxArray *nodes, struct processed_grid *grid)
/* ---------------------------------------------------------------------- */
{
    size_t i, j, nnodes;
    double *coords;

    nnodes = grid->number_of_nodes;
    coords = mxGetPr(mxGetField(nodes, 0, "coords"));

    for (j = 0; j < 3; ++j) {
        for (i = 0; i < nnodes; ++i) {
            *coords++ = grid->node_coordinates[3*i + j];
        }
    }
}


/* ---------------------------------------------------------------------- */
static mxArray *
allocate_faces(size_t nf, size_t nfacenodes)
/* ---------------------------------------------------------------------- */
{
    size_t      nflds;
    const char *fields[] = { "num", "neighbors", "nodePos", "nodes", "tag" };

    mxArray *faces, *num, *neighbors, *nodePos, *nodes, *tag;

    nflds = sizeof(fields) / sizeof(fields[0]);

    faces     = mxCreateStructMatrix(1, 1, nflds, fields);

    num       = mxCreateDoubleScalar (nf);
    neighbors = mxCreateNumericMatrix(nf        , 2, mxINT32_CLASS, mxREAL);
    nodePos   = mxCreateNumericMatrix(nf + 1    , 1, mxINT32_CLASS, mxREAL);
    nodes     = mxCreateNumericMatrix(nfacenodes, 1, mxINT32_CLASS, mxREAL);
    tag       = mxCreateNumericMatrix(nf        , 1, mxINT32_CLASS, mxREAL);

    if ((faces != NULL) && (num != NULL) && (neighbors != NULL) &&
        (nodePos != NULL) && (nodes != NULL) && (tag != NULL)) {
        mxSetField(faces, 0, "num"      , num);
        mxSetField(faces, 0, "neighbors", neighbors);
        mxSetField(faces, 0, "nodePos"  , nodePos);
        mxSetField(faces, 0, "nodes"    , nodes);
        mxSetField(faces, 0, "tag"      , tag);
    } else {
        if (tag       != NULL) { mxDestroyArray(tag);       }
        if (nodes     != NULL) { mxDestroyArray(nodes);     }
        if (nodePos   != NULL) { mxDestroyArray(nodePos);   }
        if (neighbors != NULL) { mxDestroyArray(neighbors); }
        if (num       != NULL) { mxDestroyArray(num);       }
        if (faces     != NULL) { mxDestroyArray(faces);     }

        faces = NULL;
    }

    return faces;
}


/* ---------------------------------------------------------------------- */
static void
fill_faces(mxArray *faces, struct processed_grid *grid)
/* ---------------------------------------------------------------------- */
{
    size_t i, f, nf, nfn;

    int *pi;

    nf = grid->number_of_faces;

    /* Fill faces.neighbors */
    pi = mxGetData(mxGetField(faces, 0, "neighbors"));
    for (i = 0; i < 2; i++) {
        for (f = 0; f < nf; f++) {
            /* Add one for one-based indexing in M */
            *pi++ = grid->face_neighbors[2*f + i] + 1;
        }
    }

    /* Fill faces.nodePos */
    pi = mxGetData(mxGetField(faces, 0, "nodePos"));
    for (i = 0; i <= nf; i++) { pi[i] = grid->face_ptr[i] + 1; }

    /* Fill faces.nodes */
    pi  = mxGetData(mxGetField(faces, 0, "nodes"));
    nfn = grid->face_ptr[nf];  /* Total number of face nodes */
    for (i = 0; i < nfn; i++) { pi[i] = grid->face_nodes[i] + 1; }

    /* Fill faces.tag */
    pi = mxGetData(mxGetField(faces, 0, "tag"));
    for (f = 0; f < nf; f++) { pi[f] = grid->face_tag[f] + 1; }
}


/* ---------------------------------------------------------------------- */
static size_t
count_halffaces(size_t nf, const int *neighbors)
/* ---------------------------------------------------------------------- */
{
    int    c1, c2;
    size_t nhf, f;

    for (f = nhf = 0; f < nf; f++) {
        c1 = neighbors[2*f + 0];
        c2 = neighbors[2*f + 1];

        nhf += c1 >= 0;
        nhf += c2 >= 0;
    }

    return nhf;
}


/* ---------------------------------------------------------------------- */
static mxArray *
allocate_cells(size_t nc, size_t ncf)
/* ---------------------------------------------------------------------- */
{
    size_t      nflds;
    const char *fields[] = { "num", "facePos", "faces", "indexMap" };

    mxArray *cells, *num, *facePos, *faces, *indexMap;

    nflds = sizeof(fields) / sizeof(fields[0]);

    cells    = mxCreateStructMatrix(1, 1, nflds, fields);

    num      = mxCreateDoubleScalar (nc);
    facePos  = mxCreateNumericMatrix(nc + 1, 1, mxINT32_CLASS, mxREAL);
    faces    = mxCreateNumericMatrix(ncf   , 2, mxINT32_CLASS, mxREAL);
    indexMap = mxCreateNumericMatrix(nc    , 1, mxINT32_CLASS, mxREAL);

    if ((cells != NULL) && (num != NULL) && (facePos != NULL) &&
        (faces != NULL) && (indexMap != NULL)) {
        mxSetField(cells, 0, "num"     , num     );
        mxSetField(cells, 0, "facePos" , facePos );
        mxSetField(cells, 0, "faces"   , faces   );
        mxSetField(cells, 0, "indexMap", indexMap);
    } else {
        if (indexMap != NULL) { mxDestroyArray(indexMap); }
        if (faces    != NULL) { mxDestroyArray(faces);    }
        if (facePos  != NULL) { mxDestroyArray(facePos);  }
        if (num      != NULL) { mxDestroyArray(num);      }
        if (cells    != NULL) { mxDestroyArray(cells);    }

        cells = NULL;
    }

    return cells;
}


/* ---------------------------------------------------------------------- */
static void
fill_cells(mxArray *cells, struct processed_grid *grid)
/* ---------------------------------------------------------------------- */
{
    size_t c, nc, f, nf, i;

    int c1, c2, cf_tag, nhf;
    int *pi1, *pi2;

    nc = grid->number_of_cells;
    nf = grid->number_of_faces;

    /* Simultaneously fill cells.facePos and cells.faces by transposing the
     * neighbours mapping. */
    pi1 = mxGetData(mxGetField(cells, 0, "facePos"));
    pi2 = mxGetData(mxGetField(cells, 0, "faces"  ));
    for (i = 0; i < nc + 1; i++) { pi1[i] = 0; }

    /* 1) Count connections (i.e., faces per cell). */
    for (f = 0; f < nf; f++) {
        c1 = grid->face_neighbors[2*f + 0];
        c2 = grid->face_neighbors[2*f + 1];

        if (c1 >= 0) { pi1[c1 + 1] += 1; }
        if (c2 >= 0) { pi1[c2 + 1] += 1; }
    }

    /* 2) Define start pointers (really, position *end* pointers at start). */
    for (c = 1; c <= nc; c++) {
        pi1[0] += pi1[c];
        pi1[c]  = pi1[0] - pi1[c];
    }

    /* 3) Fill connection structure whilst advancing end pointers. */
    nhf    = pi1[0];
    pi1[0] = 0;

    mxAssert (((size_t) nhf) == mxGetM(mxGetField(cells, 0, "faces")),
              "Number of half faces (SIZE(cells.faces,1)) incorrectly "
              "determined earlier.");

    for (f = 0; f < nf; f++) {
        cf_tag = 2*grid->face_tag[f] + 1;             /* [1, 3, 5] */
        c1     = grid->face_neighbors[2*f + 0];
        c2     = grid->face_neighbors[2*f + 1];

        if (c1 >= 0) {
            pi2[ pi1[ c1 + 1 ] + 0*nhf ] = f + 1;
            pi2[ pi1[ c1 + 1 ] + 1*nhf ] = cf_tag + 1;  /* out */

            pi1[ c1 + 1 ] += 1;
        }
        if (c2 >= 0) {
            pi2[ pi1[ c2 + 1 ] + 0*nhf ] = f + 1;
            pi2[ pi1[ c2 + 1 ] + 1*nhf ] = cf_tag + 0;  /* in */

            pi1[ c2 + 1 ] += 1;
        }
    }

    /* Finally, adjust pointer array for one-based indexing in M. */
    for (i = 0; i < nc + 1; i++) { pi1[i] += 1; }

    /* Fill cells.indexMap.  Note that despite the name, 'local_cell_index'
     * really *is* the (zero-based) indexMap of the 'processed_grid'. */
    pi1 = mxGetData(mxGetField(cells, 0, "indexMap"));
    for (c = 0; c < nc; c++) { pi1[c] = grid->local_cell_index[c] + 1; }
}


/* ---------------------------------------------------------------------- */
static mxArray *
allocate_grid(struct processed_grid *grid, const char *func)
/* ---------------------------------------------------------------------- */
{
    size_t nflds, nhf;
    const char *fields[] = { "nodes", "faces", "cells",
                             "type", "cartDims", "griddim" };

    mxArray *G, *nodes, *faces, *cells;
    mxArray *type, *typestr, *cartDims, *griddim;

    nflds    = sizeof(fields) / sizeof(fields[0]);
    nhf      = count_halffaces(grid->number_of_faces, grid->face_neighbors);

    G        = mxCreateStructMatrix(1, 1, nflds, fields);

    nodes    = allocate_nodes(grid->number_of_nodes);
    faces    = allocate_faces(grid->number_of_faces,
                              grid->face_ptr[ grid->number_of_faces ]);
    cells    = allocate_cells(grid->number_of_cells, nhf);
    type     = mxCreateCellMatrix(1, 1);
    typestr  = mxCreateString(func);
    cartDims = mxCreateDoubleMatrix(1, 3, mxREAL);
    griddim  = mxCreateDoubleScalar(3);

    if ((G        != NULL) && (nodes   != NULL) && (faces   != NULL) &&
        (cells    != NULL) && (type    != NULL) && (typestr != NULL) &&
        (cartDims != NULL) && (griddim != NULL)) {
        mxSetCell(type, 0, typestr);

        mxSetField(G, 0, "nodes"   , nodes   );
        mxSetField(G, 0, "faces"   , faces   );
        mxSetField(G, 0, "cells"   , cells   );
        mxSetField(G, 0, "type"    , type    );
        mxSetField(G, 0, "cartDims", cartDims);
        mxSetField(G, 0, "griddim" , griddim );
    } else {
        if (griddim  != NULL) { mxDestroyArray(griddim);  }
        if (cartDims != NULL) { mxDestroyArray(cartDims); }
        if (typestr  != NULL) { mxDestroyArray(typestr);  }
        if (type     != NULL) { mxDestroyArray(type);     }
        if (cells    != NULL) { mxDestroyArray(cells);    }
        if (faces    != NULL) { mxDestroyArray(faces);    }
        if (nodes    != NULL) { mxDestroyArray(nodes);    }
        if (G        != NULL) { mxDestroyArray(G);        }

        G = NULL;
    }

    return G;
}


/* ---------------------------------------------------------------------- */
static void
fill_grid(mxArray *G, struct processed_grid *grid)
/* ---------------------------------------------------------------------- */
{
    double *pr;

    pr    = mxGetPr(mxGetField(G, 0, "cartDims"));
    pr[0] = grid->dimensions[0];
    pr[1] = grid->dimensions[1];
    pr[2] = grid->dimensions[2];

    fill_nodes(mxGetField(G, 0, "nodes"), grid);
    fill_faces(mxGetField(G, 0, "faces"), grid);
    fill_cells(mxGetField(G, 0, "cells"), grid);
}


/* ---------------------------------------------------------------------- */
static int
args_ok(int nlhs, int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok;

    ok =       ((nlhs == 0) || (nlhs == 1));
    ok = ok && ((nrhs == 1) || (nrhs == 2));

    ok = ok && !mxIsEmpty(prhs[0]);
    ok = ok && mxIsStruct(prhs[0]);

    ok = ok && (mxGetFieldNumber(prhs[0], "cartDims") >= 0);
    ok = ok && (mxGetFieldNumber(prhs[0], "COORD"   ) >= 0);
    ok = ok && (mxGetFieldNumber(prhs[0], "ZCORN"   ) >= 0);

    if (ok && (nrhs == 2)) {
        ok = mxIsDouble(prhs[1]) && (mxGetNumberOfElements(prhs[1]) == 1);
    }

    return ok;
}


/* ---------------------------------------------------------------------- */
static double
define_tolerance(int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    double tol;

    tol = 0.0;

    if (nrhs == 2) {
        tol = mxGetScalar(prhs[1]);
    }

    return tol;
}


/* G = processgrid(grdecl)
   G = processgrid(grdecl, tolerance)
*/
/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    double                tolerance;
    char                  errmsg[1023 + 1];
    struct grdecl         grdecl;
    struct processed_grid g;

    if (args_ok(nlhs, nrhs, prhs)) {
        mx_init_grdecl(&grdecl, prhs[0]);
        tolerance = define_tolerance(nrhs, prhs);

        process_grdecl(&grdecl, tolerance, &g);

        plhs[0] = allocate_grid(&g, mexFunctionName());

        if (plhs[0] != NULL) {
            fill_grid(plhs[0], &g);
        } else {
            /* Failed to create grid structure.  Return empty. */
            plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
        }

        free_processed_grid(&g);
    } else {
        sprintf(errmsg,
                "Calling sequence is\n\t"
                "G = %s(grdecl)\t%%or\n\t"
                "G = %s(grdecl, tolerance)\n"
                "The 'grdecl' must be a valid structure with fields\n"
                "\t'cartDims', 'COORD', 'ZCORN'",
                mexFunctionName(), mexFunctionName());
        mexErrMsgTxt(errmsg);
    }
}

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
