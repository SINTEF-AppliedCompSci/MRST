#include <stddef.h>

#include <mex.h>

#include "grid.h"
#include "mrst_api.h"
#include "coarse_conn.h"


/* ---------------------------------------------------------------------- */
static int
args_ok(int nlhs, int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok;

    ok =       ((nlhs == 1) && ((nrhs == 2) || (nrhs == 3)));
    ok = ok && (mxIsStruct(prhs[0]));
    ok = ok && (mxIsDouble(prhs[1]) || mxIsInt32(prhs[1]));
    ok = ok && ((nrhs == 2) || (mxIsNumeric(prhs[2]) &&
                                (mxGetNumberOfElements(prhs[2]) == 1)));

    return ok;
}


/* ---------------------------------------------------------------------- */
static int *
extract_partition(size_t n, const mxArray *M_p)
/* ---------------------------------------------------------------------- */
{
    size_t i;
    int    *p;

    int    *pi;
    double *pd;

    n = mxGetNumberOfElements(M_p);
    p = mxMalloc(n * sizeof *p);

    if (p != NULL) {
        if (mxIsDouble(M_p)) {
            pd = mxGetPr(M_p);

            for (i = 0; i < n; i++) { p[i] = pd[i] - 1; }
        } else {
            pi = mxGetData(M_p);

            for (i = 0; i < n; i++) { p[i] = pi[i] - 1; }
        }
    }

    return p;
}


/* ---------------------------------------------------------------------- */
static void
assign_int_vec(int *v, mxArray *M_v)
/* ---------------------------------------------------------------------- */
{
    size_t i, n;

    int    *pi;
    double *pd;

    n = mxGetNumberOfElements(M_v);

    if (mxIsDouble(M_v)) {
        pd = mxGetPr(M_v);

        for (i = 0; i < n; i++) { pd[i] = v[i] + 1; }
    } else {
        pi = mxGetData(M_v);

        for (i = 0; i < n; i++) { pi[i] = v[i] + 1; }
    }
}


/* ---------------------------------------------------------------------- */
static mxArray *
generate_coarse_cells(struct coarse_topology *topo)
/* ---------------------------------------------------------------------- */
{
    mxArray *cells;
    mxArray *fld;

    const char *fields[] = {"num", "facePos", "faces"};
    int  nfields   = sizeof(fields) / sizeof(fields[0]);

    mwSize dims[] = { 1, 1 };
    mwSize ndims  = sizeof(dims) / sizeof(dims[0]);

    cells = mxCreateStructArray(ndims, dims, nfields, fields);

    if (cells != NULL) {
        fld = mxCreateDoubleScalar(topo->nblocks);
        if (fld != NULL) {
            mxSetField(cells, 0, "num", fld);
        }

        fld = mxCreateNumericMatrix(topo->nblocks + 1, 1,
                                    mxINT32_CLASS, mxREAL);

        if (fld != NULL) {
            assign_int_vec(topo->blkfacepos, fld);
            mxSetField(cells, 0, "facePos", fld);
        }

        fld = mxCreateNumericMatrix(topo->blkfacepos[topo->nblocks], 1,
                                    mxINT32_CLASS, mxREAL);
        if (fld != NULL) {
            assign_int_vec(topo->blkfaces, fld);
            mxSetField(cells, 0, "faces", fld);
        }
    }

    return cells;
}


/* ---------------------------------------------------------------------- */
static mxArray *
generate_coarse_faces(struct coarse_topology *topo)
/* ---------------------------------------------------------------------- */
{
    mxArray *faces;
    mxArray *fld;

    int    fld_no;
    size_t f;
    int   *pi;

    const char *fields[] = {"num", "neighbors", "tag"};
    int  nfields   = sizeof(fields) / sizeof(fields[0]);

    mwSize dims[] = { 1, 1 };
    mwSize ndims  = sizeof(dims) / sizeof(dims[0]);

    faces = mxCreateStructArray(ndims, dims, nfields, fields);

    if (faces != NULL) {
        fld = mxCreateDoubleScalar(topo->nfaces);
        if (fld != NULL) {
            mxSetField(faces, 0, "num", fld);
        }

        fld = mxCreateNumericMatrix(topo->nfaces, 2, mxINT32_CLASS, mxREAL);
        if (fld != NULL) {
            pi = mxGetData(fld);
            for (f = 0; f < (size_t)topo->nfaces; f++) {
                pi[f + 0*topo->nfaces] = topo->neighbours[2*f + 0] + 1;
                pi[f + 1*topo->nfaces] = topo->neighbours[2*f + 1] + 1;
            }

            mxSetField(faces, 0, "neighbors", fld);
        }

        fld = mxCreateNumericMatrix(topo->nfaces, 1, mxINT32_CLASS, mxREAL);
        if (fld != NULL) {
            pi = mxGetData(fld);
            for (f = 0; f < (size_t)topo->nfaces; f++) {
                pi[f] = 0;
            }

            mxSetField(faces, 0, "tag", fld);
        }

        if (topo->subfacepos != NULL) {
            fld_no = mxAddField(faces, "subfacePos");
            if (fld_no >= 0) {
                fld = mxCreateNumericMatrix(topo->nfaces + 1, 1,
                                            mxINT32_CLASS, mxREAL);
                if (fld != NULL) {
                    assign_int_vec(topo->subfacepos, fld);
                    mxSetFieldByNumber(faces, 0, fld_no, fld);
                }
            }
        }

        if (topo->subfaces != NULL) {
            fld_no = mxAddField(faces, "subfaces");
            if (fld_no >= 0) {
                fld = mxCreateNumericMatrix(topo->subfacepos[topo->nfaces], 1,
                                            mxINT32_CLASS, mxREAL);
                if (fld != NULL) {
                    assign_int_vec(topo->subfaces, fld);
                    mxSetFieldByNumber(faces, 0, fld_no, fld);
                }
            }
        }
    }

    return faces;
}


/* ---------------------------------------------------------------------- */
static mxArray *
generate_coarse_grid(struct coarse_topology *topo,
                     const  mxArray         *p)
/* ---------------------------------------------------------------------- */
{
    mxArray *CG;
    mxArray *fld;

    const char *fields[] = {"cells", "faces", "partition"};
    int  nfields   = sizeof(fields) / sizeof(fields[0]);

    mwSize dims[]  = { 1, 1 };
    mwSize ndims   = sizeof(dims) / sizeof(dims[0]);

    CG = mxCreateStructArray(ndims, dims, nfields, fields);

    if (CG != NULL) {
        fld = generate_coarse_cells(topo);
        if (fld != NULL) {
            mxSetField(CG, 0, "cells", fld);
        }

        fld = generate_coarse_faces(topo);
        if (fld != NULL) {
            mxSetField(CG, 0, "faces", fld);
        }

        fld = mxDuplicateArray(p);
        if (fld != NULL) {
            mxSetField(CG, 0, "partition", fld);
        }
    }

    return CG;
}


/*
 * CG = mex_generate_coarsegrid(G, p)
 * CG = mex_generate_coarsegrid(G, p, expct_nconn)
 */

/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int                     nc, nf, *p, *fneighbours;
    int                     expct_nconn;
    char                    errmsg[1023 + 1];
    struct coarse_topology *topo;

    if (args_ok(nlhs, nrhs, prhs)) {
        nc = getNumberOfCells(prhs[0]);
        nf = getNumberOfFaces(prhs[0]);
        p  = extract_partition(nc, prhs[1]);

        fneighbours = getFaceCellNeighbors(prhs[0]);

        if ((p != NULL) && (fneighbours != NULL)) {
            if (nrhs == 2) {
                expct_nconn = 0;
            } else {
                expct_nconn = mxGetScalar(prhs[2]);
            }

            if (expct_nconn < 0) {
                expct_nconn = 0;
            }

            topo = coarse_topology_create(nc, nf, expct_nconn,
                                          p, fneighbours);

            plhs[0] = generate_coarse_grid(topo, prhs[1]);

            coarse_topology_destroy(topo);

            mxFree(fneighbours);  mxFree(p);
        }
    } else {
        sprintf(errmsg,
                "Calling sequence is\n"
                "\tCG = %s(G, p) %% or\n"
                "\tCG = %s(G, p, expected_nconn)",
                mexFunctionName(), mexFunctionName());

        mexErrMsgTxt(errmsg);
    }
}
