/* "API" */

#include <stddef.h>

#include <mex.h>
#include "grid.h"
#include "mrst_api.h"

/* ------------------------------------------------------------------ */
static mxArray*
getField(const mxArray *a, const char *field, const char *subfield)
/* ------------------------------------------------------------------ */
{
    mxArray *b,*c=NULL;
    if (subfield)
    {
        b = mxGetField(a , 0, field);
        if (b!=NULL)
        {
            c = mxGetField(b, 0, subfield);
        }
    }
    else
    {
        c = mxGetField(a , 0, field);        
    }
    if (c==NULL)
    {
        mexErrMsgTxt("Internal error in mrst_api.c:getField:"
                     " attempting to access nonexistant field.");
    }

    return c;
}

/* ------------------------------------------------------------------ */
static int *
extractIntMatrix(const mxArray *a)
/* ------------------------------------------------------------------ */
{
    size_t  n = mxGetNumberOfElements(a);
    int    *q = mxMalloc(n * sizeof *q);
    if (q != NULL) {
        if (mxIsInt32(a)) {
            int *p = mxGetData(a);
            size_t i;
            for (i=0; i<n; ++i) {
                q[i] = p[i]-1;
            }
        } else if(mxIsDouble(a)) {
            double *p = mxGetPr(a);
            size_t i;
            for (i=0; i<n; ++i) {
                mxAssert ((1 <= p[i]) && (p[i] <= INT_MAX),
                          "Matrix entry exceeds INT_MAX");
                q[i] = p[i]-1;
            }
        }
    }

    return q;
}

/* ------------------------------------------------------------------ */
static int *
extractIntMatrixTranspose(const mxArray *a)
/* ------------------------------------------------------------------ */
{
    int *p = extractIntMatrix(a);
    size_t M = mxGetM(a);
    size_t N = mxGetN(a);

    int *q = mxMalloc(M * N * sizeof *q);
    if (q != NULL) {
        size_t i,j;
        for(i=0; i<M; ++i) {
            for(j=0; j<N; ++j) {
                q[i*N+j] = p[i+M*j];
            }
        }
    }
    mxFree(p);
    return q;
}

/* ------------------------------------------------------------------ */
static double *
extractDoubleMatrixTranspose(const mxArray *a)
/* ------------------------------------------------------------------ */
{
    size_t M  = mxGetM(a);
    size_t N  = mxGetN(a);
    double *q = mxMalloc(M * N * sizeof *q);
    if (q != NULL) {
        double *p = mxGetPr(a);
        size_t i,j;
        for(i=0; i<M; ++i) {
            for(j=0; j<N; ++j) {
                q[i*N+j] = p[i+M*j];
            }
        }
    }
    return q;
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


/* ------------------------------------------------------------------ */
int
getNumberOfDimensions(const mxArray *G)
/* ------------------------------------------------------------------ */
{
    return mxGetN(getField(G, "nodes", "coords"));
}

/* ------------------------------------------------------------------ */
void
getLocal2GlobalCellMap(const mxArray *G)
/* ------------------------------------------------------------------ */
{
    (void)G;
    mxAssert(0, "Not implemented!");
}

/* ------------------------------------------------------------------ */
int
getNumberOfNodes(const mxArray *G)
/* ------------------------------------------------------------------ */
{
    return mxGetM(getField(G, "nodes", "coords"));
}

/* ------------------------------------------------------------------ */
double *
getNodeCoordinates(const mxArray *G)
/* ------------------------------------------------------------------ */
{
    return extractDoubleMatrixTranspose(getField(G, "nodes", "coords"));
#if 0
    mxArray *p1, *p2;

    p1 = mxGetField(G , 0, "nodes" );
    p2 = mxGetField(p1, 0, "coords");

    const int n = getNumberOfNodes(G);
    const int d = getNumberOfDimensions(G);

    double *v = mxMalloc(n * d * sizeof *v);
    if (v != NULL) {

        double *tmp = mxGetPr(p2);
        int i,j;
        for (i=0; i<n; ++i) {
            for(j=0; j<d; ++j) {
                v[d*i+j] = tmp[i + n*j];
            }
        }
    }
    return v;
#endif
}

/* ------------------------------------------------------------------ */
int
getNumberOfFaces(const mxArray *G)
/* ------------------------------------------------------------------ */
{
    return mxGetNumberOfElements(getField(G, "faces", "nodePos"))-1;
}

/* ------------------------------------------------------------------ */
int *
getFaceNodePos(const mxArray *G)
/* ------------------------------------------------------------------ */
{
    return extractIntMatrix(getField(G, "faces", "nodePos"));
}

/* ------------------------------------------------------------------ */
int
getNumberOfFaceNodes(const mxArray *G)
/* ------------------------------------------------------------------ */
{
    return mxGetNumberOfElements(getField(G, "faces", "nodes"));
}

/* ------------------------------------------------------------------ */
int *
getFaceNodes(const mxArray *G)
/* ------------------------------------------------------------------ */
{
    return extractIntMatrix(getField(G, "faces", "nodes"));
}


/* ------------------------------------------------------------------ */
int *
getFaceCellNeighbors(const mxArray *G)
/* ------------------------------------------------------------------ */
{
    return extractIntMatrixTranspose(getField(G, "faces", "neighbors"));
}

/* ------------------------------------------------------------------ */
double *
getFaceAreas(const mxArray *G)
/* ------------------------------------------------------------------ */
{
    mxArray *arr = getField(G, "faces", "areas");
    double  *values = NULL;
    if (arr!=NULL)
    {
      values = mxGetPr(arr);
    }
    return values;
}

/* ------------------------------------------------------------------ */
double *
getFaceNormals(const mxArray *G)
/* ------------------------------------------------------------------ */
{
    mxArray *arr = getField(G, "faces", "normals");
    double  *values = NULL;
    if (arr!=NULL)
    {
      values = extractDoubleMatrixTranspose(arr);
    }
    return values;
}

/* ------------------------------------------------------------------ */
double *
getFaceCentroids(const mxArray *G)
/* ------------------------------------------------------------------ */
{
    mxArray *arr = getField(G, "faces", "centroids");
    double  *values = NULL;
    if (arr!=NULL)
    {
      values = extractDoubleMatrixTranspose(arr);
    }
    return values;
}

/* ------------------------------------------------------------------ */
int
getNumberOfCells(const mxArray *G)
/* ------------------------------------------------------------------ */
{
    return mxGetNumberOfElements(getField(G, "cells", "facePos"))-1;
}

/* ------------------------------------------------------------------ */
int *getCellFacePos(const mxArray *G)
/* ------------------------------------------------------------------ */
{
    return extractIntMatrix(getField(G, "cells", "facePos"));
}

/* ------------------------------------------------------------------ */
int getNumberOfCellFaces(const mxArray *G)
/* ------------------------------------------------------------------ */
{
    return mxGetM(getField(G, "cells", "faces"));
}

/* ------------------------------------------------------------------ */
int *getCellFaces(const mxArray *G)
/* ------------------------------------------------------------------ */
{
    return extractIntMatrix(getField(G, "cells", "faces"));
}

/* ------------------------------------------------------------------ */
double *
getCellVolumes(const mxArray *G)
/* ------------------------------------------------------------------ */
{
    mxArray *arr = getField(G, "cells", "volumes");
    double  *values = NULL;
    if (arr!=NULL)
    {
      values = mxGetPr(arr);
    }
    return values;
}



/* ------------------------------------------------------------------ */
double *
getCellCentroids(const mxArray *G)
/* ------------------------------------------------------------------ */
{
    mxArray *arr = getField(G, "cells", "centroids");
    double  *values = NULL;
    if (arr!=NULL)
    {
      values = extractDoubleMatrixTranspose(arr);
    }
    return values;
}



/*
 *
 *
 *
 */


/* ------------------------------------------------------------------ */
static int
allocate_perm_data(int ncells, int d, double **K)
/* ------------------------------------------------------------------ */
{
    int     ret;
    size_t  alloc_sz;
    double *perm;

    alloc_sz  = d * d;
    alloc_sz *= ncells;

    perm = mxMalloc(alloc_sz * sizeof *perm);

    if (perm != NULL) {
        *K = perm;
        ret = 1;
    } else {
        *K = NULL;
        ret = 0;
    }

    return ret;
}


/* ------------------------------------------------------------------ */
double *
getPermeability(const mxArray *perm, int d)
/* ------------------------------------------------------------------ */
{
    int    ncomp, alloc_ok;
    size_t ncells, c, i, off;

    double *k, *tensor;

    ncells = mxGetM(perm);
    ncomp  = mxGetN(perm);

    alloc_ok = allocate_perm_data(ncells, d, &k);

    if (alloc_ok) {
        for (i = 0; i < (size_t) ((ncells * d) * d); i++) {
            k[i] = 0.0;
        }

        tensor = mxGetPr(perm);

        if (ncomp == 1) {
            /* Isotropic (scalar) tensor */
            for (c = off = 0; c < ncells; c++, off += d * d) {
                for (i = 0; i < (size_t) d; i++) {
                    k[i*(d + 1) + off] = tensor[c];
                }
            }
        } else if (ncomp == d) {
            /* Diagonal tensor */
            for (c = off = 0; c < ncells; c++, off += d * d) {
                for (i = 0; i < (size_t) d; i++) {
                    k[i*(d + 1) + off] = tensor[c + i*ncells];
                }
            }
        } else if (d == 2) {
            /* Full 2D tensor */
            mxAssert (ncomp == 3, "");

            for (c = off = 0; c < ncells; c++, off += d * d) {
                k[0 + off] = tensor[c + 0*ncells];
                k[1 + off] = tensor[c + 1*ncells];
                k[2 + off] = tensor[c + 1*ncells];
                k[3 + off] = tensor[c + 2*ncells];
            }
        } else {
            /* Full 3D tensor */
            mxAssert ((d == 3) && (ncomp == 6), "");

            for (c = off = 0; c < ncells; c++, off += d * d) {
                k[0 + off] = tensor[c + 0*ncells];
                k[1 + off] = tensor[c + 1*ncells];
                k[2 + off] = tensor[c + 2*ncells];

                k[3 + off] = tensor[c + 1*ncells];
                k[4 + off] = tensor[c + 3*ncells];
                k[5 + off] = tensor[c + 4*ncells];

                k[6 + off] = tensor[c + 2*ncells];
                k[7 + off] = tensor[c + 4*ncells];
                k[8 + off] = tensor[c + 5*ncells];
            }
        }
    }
    return k;
}

/* ---------------------------------------------------------------------- */
int
verify_mrst_grid(const mxArray *G)
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
void
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
grid_t *
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
            g->dimensions       = getNumberOfDimensions(G);
            g->number_of_cells  = getNumberOfCells     (G);
            g->number_of_faces  = getNumberOfFaces     (G);

            g->cell_volumes     = getCellVolumes(G);
            g->face_areas       = getFaceAreas  (G);

            g->face_nodepos     = NULL; /* Not accessed */
            g->face_nodes       = NULL; /* Not accessed */
            g->node_coordinates = NULL; /* Not accessed */
        }
    }

    return g;
}

/* ---------------------------------------------------------------------- */
grid_t *
mrst_grid_topology(const mxArray *G)
/* ---------------------------------------------------------------------- */
{
    int     topo_ok;
    grid_t *g;

    g = mxMalloc(1 * sizeof *g);

    if (g != NULL) {
        /* Topology description should never be NULL */
        topo_ok  = (g->face_cells     = getFaceCellNeighbors(G)) != NULL;
        topo_ok += (g->cell_faces     = getCellFaces        (G)) != NULL;
        topo_ok += (g->cell_facepos   = getCellFacePos      (G)) != NULL;

        /* Put the (purportedly) allocated geometry fields in a
         * consistent initial state lest free_mrst_grid() fail
         * spectacularly (e.g., free()ing a random pointer) in case of
         * memory allocation failure above (i.e., if topo_ok != 3). */

        g->face_centroids = NULL; /* Not accessed */
        g->face_normals   = NULL; /* Not accessed */
        g->cell_centroids = NULL; /* Not accessed */

        if (topo_ok != 3) {
            free_mrst_grid(g);
            g = NULL;
        } else {
            g->dimensions       = getNumberOfDimensions(G);
            g->number_of_cells  = getNumberOfCells     (G);
            g->number_of_faces  = getNumberOfFaces     (G);

            g->cell_volumes     = NULL; /* Not accessed */
            g->face_areas       = NULL; /* Not accessed */

            g->face_nodepos     = NULL; /* Not accessed */
            g->face_nodes       = NULL; /* Not accessed */
            g->node_coordinates = NULL; /* Not accessed */
        }
    }

    return g;
}

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
