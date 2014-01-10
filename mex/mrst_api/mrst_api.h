#ifndef MRST_API_H_INCLUDED
#define MRST_API_H_INCLUDED

/* We preserve the old typedef here, since it is gone from grid.h.
 * This header is not standalone, grid.h (and more) must be included
 * before mrst_api.h
 */
typedef struct UnstructuredGrid grid_t;

#ifdef __cplusplus
extern "C" {
#endif

/* 
 *  "API" to MRST grid : implements access to raw C vectors.
 */

int  getNumberOfDimensions (const mxArray *G);
void getLocal2GlobalCellMap(const mxArray *G);

/* Node coordinates */
int     getNumberOfNodes      (const mxArray *G);
double *getNodeCoordinates(const mxArray *G); /* copy */

/* Face topology */
int     getNumberOfFaces      (const mxArray *G);
int     getNumberOfFaceNodes  (const mxArray *G);
int    *getFaceNodePos        (const mxArray *G);    /* copy */
int    *getFaceNodes          (const mxArray *G);    /* copy */
int    *getFaceCellNeighbors  (const mxArray *G);    /* copy */

/* Face geometry */
double *getFaceAreas          (const mxArray *G);
double *getFaceNormals        (const mxArray *G);    /* copy */
double *getFaceCentroids      (const mxArray *G);    /* copy */

/* Cell topology */
int     getNumberOfCells      (const mxArray *G);
int     getNumberOfCellFaces  (const mxArray *G);
int    *getCellFacePos        (const mxArray *G);    /* copy */
int    *getCellFaces          (const mxArray *G);    /* copy */

/* Cell geometry */
double *getCellVolumes        (const mxArray *G);
double *getCellCentroids      (const mxArray *G);    /* copy */


/* Rock properties */
double *
getPermeability(const mxArray *perm, int d);    /* copy */


/* Grid stuff */
grid_t *mrst_grid(const mxArray *G);
grid_t *mrst_grid_topology(const mxArray *G);
void    free_mrst_grid(grid_t *g);

int verify_mrst_grid(const mxArray *G);


#ifdef __cplusplus
}
#endif


#endif /* MRST_API_H_INCLUDED */
