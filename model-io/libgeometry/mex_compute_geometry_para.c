/*
 * Copyright 2010 (c) SINTEF ICT, Applied Mathematics.
 * Jostein R. Natvig <Jostein.R.Natvig at sintef.no>
 */
#include <stdlib.h>
#include <stdio.h>

#include <mex.h>

#include "grid.h"
#include "mrst_api.h"

#include "geometry.h"


/* ------------------------------------------------------------------ */
static int
verify_faces_structure(mxArray *faces)
/* ------------------------------------------------------------------ */
{
    /* Shallow structural inspection only.  Assume valid fields... */
    int ok;

    ok =       (mxGetFieldNumber(faces, "neighbors") >= 0);
    ok = ok && (mxGetFieldNumber(faces, "nodePos"  ) >= 0);
    ok = ok && (mxGetFieldNumber(faces, "nodes"    ) >= 0);

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
args_ok(int nlhs, int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    return (nlhs == 5) && (nrhs == 2) && verify_grid(prhs[0]);
}


/* ---------------------------------------------------------------------- */
static void
compute_geometry(const mxArray *G ,
                 const int      d ,
                 const int      nc,
                 const int      nf,
                 double        *a ,
                 double        *fc,
                 double        *fn,
                 double        *cc,
                 double        *cv)
/* ---------------------------------------------------------------------- */
{
    int    *nodepos, *facenodes, *cellfaces, *facepos, *neighbors;
    double *coords;

    mxAssert(d == 3, "Sorry, only support for 3D grids.");

    /* Grid topology: */
    nodepos   = getFaceNodePos(G);
    facenodes = getFaceNodes(G);
    cellfaces = getCellFaces(G);
    facepos   = getCellFacePos(G);
    neighbors = getFaceCellNeighbors(G);
    coords    = getNodeCoordinates(G);



    /* Compute face geometry */
    compute_face_geometry(d, coords, nf, nodepos, facenodes, fn, fc, a);


    /* Compute cell geometry */
    compute_cell_geometry(d, coords, nodepos, facenodes, neighbors,
                          fn, fc, nc, facepos, cellfaces, cc, cv);


    /* Clean up */
    mxFree(coords);   mxFree(facenodes);
    mxFree(nodepos);  mxFree(cellfaces);  mxFree(facepos);
}


/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int nc, nf, nd;
    int flag;
    double dtmp;
    const mxArray *G;
    mxArray       *fa, *fc, *fn, *cc, *cv;

    char errmsg[1023 + 1];

    if (args_ok(nlhs, nrhs, prhs)) {
        G  = prhs[0];
	flag = (int)  mxGetScalar(prhs[1]);
	/* mexPrintf("Using %i threads\n",flag); */
	omp_set_num_threads(flag);
        nc = getNumberOfCells(G);
        nf = getNumberOfFaces(G);
        nd = getNumberOfDimensions(G);

        fa = mxCreateDoubleMatrix(nf,  1, mxREAL); /* Face area */
        fc = mxCreateDoubleMatrix(nd, nf, mxREAL); /* Face centroid */
        fn = mxCreateDoubleMatrix(nd, nf, mxREAL); /* Face normal */
        cc = mxCreateDoubleMatrix(nd, nc, mxREAL); /* Cell centroid */
        cv = mxCreateDoubleMatrix(nc,  1, mxREAL); /* Cell volume */

        compute_geometry(G, nd, nc, nf,
                         mxGetPr(fa), mxGetPr(fc), mxGetPr(fn),
                         mxGetPr(cc), mxGetPr(cv));

        plhs[0] = fa;
        plhs[1] = fc;
        plhs[2] = fn;
        plhs[3] = cc;
        plhs[4] = cv;
    } else {
        sprintf(errmsg,
                "Calling sequence is\n\t"
                "[fa, fc, fn, cc, cv] = %s(G)",
                mexFunctionName());
        mexErrMsgTxt(errmsg);
    }
}
