/*===========================================================================
//
// File: preprocess.c
//
// Created: Fri Jun 19 08:42:39 2009
//
// Author: Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
//
// $Date: 2012-11-15 00:18:05 +0100 (Thu, 15 Nov 2012) $
//
// $Revision: 1010 $
//
//==========================================================================*/

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

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "preprocess.h"
#include "uniquepoints.h"
#include "facetopology.h"


#define MIN(i,j) ((i)<(j) ? (i) : (j))
#define MAX(i,j) ((i)>(j) ? (i) : (j))

static void
compute_cell_index(const int dims[3], int i, int j, int *neighbors, int len);

static int
checkmemory(int nz, struct processed_grid *out, int **intersections);

static void
process_vertical_faces(int direction,
                       int **intersections,
                       int *plist, int *work,
                       struct processed_grid *out);

static void
process_horizontal_faces(int **intersections,
                         int *plist,
                         struct processed_grid *out);

static int
linearindex(const int dims[3], int i, int j, int k)
{
    assert (0 <= i);
    assert (0 <= j);
    assert (0 <= k);

    assert (i < dims[0]);
    assert (j < dims[1]);
    assert (k < dims[2]);

    return i + dims[0]*(j + dims[1]*k);
}


/*-----------------------------------------------------------------
  Given a vector <field> with k index running faster than i running
  faster than j, and Cartesian dimensions <dims>, find pointers to the
  (i-1, j-1, 0), (i-1, j, 0), (i, j-1, 0) and (i, j, 0) elements of
  field.  */
static void
igetvectors(int dims[3], int i, int j, int *field, int *v[])
{
    int im = MAX(1,       i  ) - 1;
    int ip = MIN(dims[0], i+1) - 1;
    int jm = MAX(1,       j  ) - 1;
    int jp = MIN(dims[1], j+1) - 1;

    v[0] = field + dims[2]*(im + dims[0]* jm);
    v[1] = field + dims[2]*(im + dims[0]* jp);
    v[2] = field + dims[2]*(ip + dims[0]* jm);
    v[3] = field + dims[2]*(ip + dims[0]* jp);
}






/*-----------------------------------------------------------------
  Special purpose

  Convert from k-index to Cartesian index i+nx*(j+ny*k) for every
  other element in neighbors.

*/
static void
compute_cell_index(const int dims[3], int i, int j,
                   int *neighbors, int len)
{
    int k;

    if (((i < 0) || (i >= dims[0])) || /* 'i' outside [0, dims[0]) */
        ((j < 0) || (j >= dims[1]))) { /* 'j' outside [0, dims[1]) */

        for (k = 0; k < len; k += 2) {
            neighbors[k] = -1;  /* Neighbour is outside domain */
        }
    }
    else {
        for (k = 0; k < len; k += 2) {
            if (neighbors[k] != -1) {
                neighbors[k] = linearindex(dims, i, j, neighbors[k]);
            }
        }
    }
}


/*-----------------------------------------------------------------
  Ensure there's sufficient memory */
static int
checkmemory(int nz, struct processed_grid *out, int **intersections)
{
    int r, m, n, ok;

    /* Ensure there is enough space to manage the (pathological) case
     * of every single cell on one side of a fault connecting to all
     * cells on the other side of the fault (i.e., an all-to-all cell
     * connectivity pairing). */
    r = (2*nz + 2) * (2*nz + 2);
    m = out->m;
    n = out->n;

    if (out->number_of_faces +  r > m) {
        m += MAX(m / 2,  2 * r);
    }
    if (out->face_ptr[out->number_of_faces] + 6*r > n) {
        n += MAX(n / 2, 12 * r);
    }

    ok = m == out->m;
    if (! ok) {
        void *p1, *p2, *p3, *p4;

        p1 = realloc(*intersections     , 4*m   * sizeof **intersections);
        p2 = realloc(out->face_neighbors, 2*m   * sizeof *out->face_neighbors);
        p3 = realloc(out->face_ptr      , (m+1) * sizeof *out->face_ptr);
        p4 = realloc(out->face_tag      , 1*m   * sizeof *out->face_tag);

        if (p1 != NULL) { *intersections      = p1; }
        if (p2 != NULL) { out->face_neighbors = p2; }
        if (p3 != NULL) { out->face_ptr       = p3; }
        if (p4 != NULL) { out->face_tag       = p4; }

        ok = (p1 != NULL) && (p2 != NULL) && (p3 != NULL) && (p4 != NULL);

        if (ok) { out->m = m; }
    }

    if (ok && (n != out->n)) {
        void *p1;

        p1 = realloc(out->face_nodes, n * sizeof *out->face_nodes);

        ok = p1 != NULL;

        if (ok) {
            out->face_nodes = p1;
            out->n          = n;
        }
    }

    return ok;
}

/*-----------------------------------------------------------------
  For each vertical face (i.e. i or j constant),
  -find point numbers for the corners and
  -cell neighbors.
  -new points on faults defined by two intgersecting lines.

  direction == 0 : constant-i faces.
  direction == 1 : constant-j faces.
*/
static void
process_vertical_faces(int direction,
                       int **intersections,
                       int *plist, int *work,
                       struct processed_grid *out)
{
    int i,j;
    int *cornerpts[4];
    int d[3];
    int f;
    enum face_tag tag[] = { LEFT, BACK };
    int *tmp;
    int nx = out->dimensions[0];
    int ny = out->dimensions[1];
    int nz = out->dimensions[2];
    int startface;
    int num_intersections;
    int *ptr;
    int len;

    assert ((direction == 0) || (direction == 1));

    d[0] = 2 * (nx + 0);
    d[1] = 2 * (ny + 0);
    d[2] = 2 * (nz + 1);

    for (j = 0; j < ny + direction; ++j) {
        for (i = 0; i < nx + (1 - direction); ++i) {

            if (! checkmemory(nz, out, intersections)) {
                fprintf(stderr,
                        "Could not allocate enough space in "
                        "process_vertical_faces()\n");
                exit(1);
            }

            /* Vectors of point numbers */
            igetvectors(d, 2*i + direction, 2*j + (1 - direction),
                        plist, cornerpts);

            if (direction == 1) {
                /* 1   3       0   1    */
                /*       --->           */
                /* 0   2       2   3    */
                /* rotate clockwise     */
                tmp          = cornerpts[1];
                cornerpts[1] = cornerpts[0];
                cornerpts[0] = cornerpts[2];
                cornerpts[2] = cornerpts[3];
                cornerpts[3] = tmp;
            }

            /* int startface = ftab->position; */
            startface = out->number_of_faces;
            /* int num_intersections = *npoints - npillarpoints; */
            num_intersections = out->number_of_nodes -
                out->number_of_nodes_on_pillars;

            /* Establish new connections (faces) along pillar pair. */
            findconnections(2*nz + 2, cornerpts,
                            *intersections + 4*num_intersections,
                            work, out);

            /* Start of ->face_neighbors[] for this set of connections. */
            ptr = out->face_neighbors + 2*startface;

            /* Total number of cells (both sides) connected by this
             * set of connections (faces). */
            len = 2*out->number_of_faces - 2*startface;

            /* Derive inter-cell connectivity (i.e. ->face_neighbors)
             * of global (uncompressed) cells for this set of
             * connections (faces). */
            compute_cell_index(out->dimensions, i-1+direction, j-direction, ptr    , len);
            compute_cell_index(out->dimensions, i            , j          , ptr + 1, len);

            /* Tag the new faces */
            f = startface;
            for (; f < out->number_of_faces; ++f) {
                out->face_tag[f] = tag[direction];
            }
        }
    }
}


/*-----------------------------------------------------------------
  For each horizontal face (i.e. k constant),
  -find point numbers for the corners and
  -cell neighbors.

  Also define map from logically Cartesian
  cell index to local cell index 0, ..., #<active cells>.   Exclude
  cells that are have collapsed coordinates. (This includes cells with
  ACTNUM==0)

*/
static void
process_horizontal_faces(int **intersections,
                         int *plist,
                         struct processed_grid *out)
{
    int i,j,k;

    int nx = out->dimensions[0];
    int ny = out->dimensions[1];
    int nz = out->dimensions[2];

    int *cell  = out->local_cell_index;
    int cellno = 0;
    int *f, *n, *c[4];
    int prevcell, thiscell;
    int idx;

    /* dimensions of plist */
    int  d[3];
    d[0] = 2*nx;
    d[1] = 2*ny;
    d[2] = 2+2*nz;


    for(j=0; j<ny; ++j) {
        for (i=0; i<nx; ++i) {


            if (! checkmemory(nz, out, intersections)) {
                fprintf(stderr,
                        "Could not allocate enough space in "
                        "process_horizontal_faces()\n");
                exit(1);
            }


            f = out->face_nodes     + out->face_ptr[out->number_of_faces];
            n = out->face_neighbors + 2*out->number_of_faces;


            /* Vectors of point numbers */
            igetvectors(d, 2*i+1, 2*j+1, plist, c);

            prevcell = -1;


            for (k = 1; k<nz*2+1; ++k){

                /* Skip if space between face k and face k+1 is collapsed. */
                /* Note that inactive cells (with ACTNUM==0) have all been  */
                /* collapsed in finduniquepoints.                           */
                if (c[0][k] == c[0][k+1] && c[1][k] == c[1][k+1] &&
                    c[2][k] == c[2][k+1] && c[3][k] == c[3][k+1]){

                    /* If the pinch is a cell: */
                    if (k%2){
                        idx = linearindex(out->dimensions, i,j,(k-1)/2);
                        cell[idx] = -1;
                    }
                }
                else{

                    if (k%2){
                        /* Add face */
                        *f++ = c[0][k];
                        *f++ = c[2][k];
                        *f++ = c[3][k];
                        *f++ = c[1][k];

                        out->face_tag[  out->number_of_faces] = TOP;
                        out->face_ptr[++out->number_of_faces] = f - out->face_nodes;

                        thiscell = linearindex(out->dimensions, i,j,(k-1)/2);
                        *n++ = prevcell;
                        *n++ = prevcell = thiscell;

                        cell[thiscell] = cellno++;

                    }
                    else{
                        if (prevcell != -1){
                            /* Add face */
                            *f++ = c[0][k];
                            *f++ = c[2][k];
                            *f++ = c[3][k];
                            *f++ = c[1][k];

                            out->face_tag[  out->number_of_faces] = TOP;
                            out->face_ptr[++out->number_of_faces] = f - out->face_nodes;

                            *n++ = prevcell;
                            *n++ = prevcell = -1;
                        }
                    }
                }
            }
        }
    }
    out->number_of_cells = cellno;
}


/*-----------------------------------------------------------------
  On input,
  L points to 4 ints that indirectly refers to points in c.
  c points to array of coordinates [x0,y0,z0,x1,y1,z1,...,xn,yn,zn].
  pt points to array of 3 doubles.

  On output,
  pt holds coordinates to intersection between lines given by point
  numbers L[0]-L[1] and L[2]-L[3].
*/
static void approximate_intersection_pt(int *L, double *c, double *pt)
{
    double a;
    double z0, z1, z2, z3;
    double b1, b2;
    double x1, y1;
    double x2, y2;
    double z;

    /* no intersection on pillars expected here! */
    assert (L[0] != L[2]);
    assert (L[1] != L[3]);

    z0 = c[3*L[0] + 2];
    z1 = c[3*L[1] + 2];
    z2 = c[3*L[2] + 2];
    z3 = c[3*L[3] + 2];

    /* find parameter a where lines L0L1 and L2L3 have same
     * z-coordinate */
    if (fabs((z1 - z0) - (z3 - z2)) > 0.0) {

        a = (z2 - z0) / ((z1 - z0) - (z3 - z2));

    } else {

        a = 0;

    }

    /* the corresponding z-coordinate is */
    z =  z0*(1.0 - a) + z1*a;


    /* find point (x1, y1, z) on pillar 1 */
    b1 = (z2 - z) / (z2 - z0);
    b2 = (z - z0) / (z2 - z0);
    x1 = c[3*L[0] + 0]*b1 + c[3*L[2] + 0]*b2;
    y1 = c[3*L[0] + 1]*b1 + c[3*L[2] + 1]*b2;

    /* find point (x2, y2, z) on pillar 2 */
    b1 = (z - z3) / (z1 - z3);
    b2 = (z1 - z) / (z1 - z3);
    x2 = c[3*L[1] + 0]*b1 + c[3*L[3] + 0]*b2;
    y2 = c[3*L[1] + 1]*b1 + c[3*L[3] + 1]*b2;

    /* horizontal lines are by definition ON the bilinear surface
       spanned by L0, L1, L2 and L3.  find point (x, y, z) on
       horizontal line between point (x1, y1, z) and (x2, y2, z).*/
    pt[0] = x1*(1.0 - a) + x2*a;
    pt[1] = y1*(1.0 - a) + y2*a;
    pt[2] = z;
}

/*-----------------------------------------------------------------
  Compute x,y and z coordinates for points on each pillar.  Then,
  append x,y and z coordinates for extra points on faults.  */
static void
compute_intersection_coordinates(int                   *intersections,
                                 struct processed_grid *out)
{
    int n  = out->number_of_nodes;
    int np = out->number_of_nodes_on_pillars;
    int    k;
    double *pt;
    int    *itsct = intersections;
    /* Make sure the space allocated for nodes match the number of
     * node. */
    void *p = realloc (out->node_coordinates, 3*n*sizeof(double));
    if (p) {
        out->node_coordinates = p;
    }
    else {
        fprintf(stderr, "Could not allocate extra space for intersections\n");
    }


    /* Append intersections */
    pt    = out->node_coordinates + 3*np;

    for (k=np; k<n; ++k){
        approximate_intersection_pt(itsct, out->node_coordinates, pt);
        pt    += 3;
        itsct += 4;

    }
}


/* ------------------------------------------------------------------ */
static int*
copy_and_permute_actnum(int nx, int ny, int nz, const int *in, int *out)
/* ------------------------------------------------------------------ */
{
    int i,j,k;
    int *ptr = out;

    /* Permute actnum such that values of each vertical stack of cells
     * are adjacent in memory, i.e.,
     *
     *    out = [in(0,0,:), in(1,0,:),..., in(nx-1, ny-1,:)]
     *
     * in MATLAB pseudo-code.
     */
    if (in != NULL) {
        for (j = 0; j < ny; ++j) {
            for (i = 0; i < nx; ++i) {
                for (k = 0; k < nz; ++k) {
                    *ptr++ = in[i + nx*(j + ny*k)];
                }
            }
        }
    }
    else {
        /* No explicit ACTNUM.  Assume all cells active. */
        for (i = 0; i < nx * ny * nz; i++) {
            out[ i ] = 1;
        }
    }

    return out;
}

/* ------------------------------------------------------------------ */
static double*
copy_and_permute_zcorn(int nx, int ny, int nz, const double *in,
                       double sign, double *out)
/* ------------------------------------------------------------------ */
{
    int i,j,k;
    double *ptr = out;
    /* Permute zcorn such that values of each vertical stack of cells
     * are adjacent in memory, i.e.,

     out = [in(0,0,:), in(1,0,:),..., in(2*nx-1, 2*ny-1,:)]

     in Matlab pseudo-code.
    */
    for (j=0; j<2*ny; ++j){
        for (i=0; i<2*nx; ++i){
            for (k=0; k<2*nz; ++k){
                *ptr++ = sign * in[i+2*nx*(j+2*ny*k)];
            }
        }
    }
    return out;
}

/* ------------------------------------------------------------------ */
static int
get_zcorn_sign(int nx, int ny, int nz, const int *actnum,
               const double *zcorn, int *error)
/* ------------------------------------------------------------------ */
{
    /* Ensure that zcorn (i.e., depth) is strictly nondecreasing in
       the k-direction.  This is required by the processign algorithm.

       1) if  z(i,j,k) <= z(i,j,k+1) for all (i,j,k), return 1.0

       2) if -z(i,j,k) <=-z(i,j,k+1) for all (i,j,k), return -1.0

       3) if (1) and (2) fails, return -1.0, and set *error = 1.

    */
    int    sign;
    int    i, j, k;
    int    c1, c2;
    double z1, z2;

    for (sign = 1; sign>-2; sign = sign - 2)
    {
        *error = 0;

        for (j=0; j<2*ny; ++j){
            for (i=0; i<2*nx; ++i){
                for (k=0; k<2*nz-1; ++k){
                    z1 = sign*zcorn[i+2*nx*(j+2*ny*(k))];
                    z2 = sign*zcorn[i+2*nx*(j+2*ny*(k+1))];

                    c1 = i/2 + nx*(j/2 + ny*(k/2));
                    c2 = i/2 + nx*(j/2 + ny*((k+1)/2));

                    if (((actnum == NULL) ||
                         (actnum[c1] && actnum[c2]))
                        && (z2 < z1)) {

                        fprintf(stderr, "\nZCORN should be strictly "
                                "nondecreasing along pillars!\n");
                        *error = 1;
                        goto end;
                    }
                }
            }
        }

    end:
        if (!*error){
            break;
        }
    }

    if (*error){
        fprintf(stderr, "Attempt to reverse sign in ZCORN failed.\n"
                "Grid definition may be broken\n");
    }

    return sign;
}


static void
ind2sub(const size_t nx,
        const size_t ny,
        const size_t nz,
        size_t       c ,
        size_t *i, size_t *j, size_t *k)
{
    assert (c < (nx * ny * nz));

#if defined(NDEBUG)
    (void) nz;
#endif

    *i = c % nx;  c /= nx;
    *j = c % ny;
    *k = c / ny;
}


/* ---------------------------------------------------------------------- */
static double
vert_size(const struct grdecl *in,
          const size_t         c ,
          const size_t         off[8])
/* ---------------------------------------------------------------------- */
{
    size_t        i, j, k, nx, ny, start;
    double        dz;
    const double *zcorn;

    nx = in->dims[ 0 ];
    ny = in->dims[ 1 ];

    ind2sub(nx, ny, in->dims[ 2 ], c, &i, &j, &k);

    zcorn = in->zcorn;
    start = (2 * i) + (2 * nx)*((2 * j) + (2 * ny)*(2 * k));

    for (k = 0, dz = 0.0; (! (fabs(dz) > 0)) && (k < 4); k++) {
        dz = zcorn[start + off[k + 4]] - zcorn[start + off[k]];
    }

    return dz;
}


/* ---------------------------------------------------------------------- */
static int
is_lefthanded(const struct grdecl *in)
/* ---------------------------------------------------------------------- */
{
    int           active, searching;
    size_t        nx, ny, nz, c;
    size_t        origin, imax, jmax;
    size_t        off[8];
    double        dx[2], dy[2], dz, triple;
    const double *pt_coord;

    nx = in->dims[0];
    ny = in->dims[1];
    nz = in->dims[2];

    off[0] = 0;
    off[1] = off[0] + 1;
    off[2] = off[0] + (2 * nx);
    off[3] = off[2] + 1;
    off[4] = off[0] + ((2 * nx) * (2 * ny));
    off[5] = off[4] + 1;
    off[6] = off[4] + (2 * nx);
    off[7] = off[6] + 1;

    pt_coord = in->coord;

    origin = 0;
    imax   = (nx + 0) * 1        * (2 * 3);
    jmax   = (nx + 1) * (ny + 0) * (2 * 3);

    dx[0] = pt_coord[imax + 0] - pt_coord[origin + 0];
    dy[0] = pt_coord[imax + 1] - pt_coord[origin + 1];

    dx[1] = pt_coord[jmax + 0] - pt_coord[origin + 0];
    dy[1] = pt_coord[jmax + 1] - pt_coord[origin + 1];

    c = 0;  dz = 0.0;
    do {
        active = (in->actnum == NULL) || (in->actnum[c] != 0);

        if (active) {
            dz = vert_size(in, c, off);
        }

        searching = ! (active && (fabs(dz) > 0.0));

        c += 1;
    } while (searching && (c < (nx * ny * nz)));

    assert (! searching);       /* active && (fabs(dz) > 0) */

    /* Compute vector triple product to distinguish left-handed (<0)
     * from right-handed (>0) coordinate systems. */
    triple = dz * (dx[0]*dy[1] - dx[1]*dy[0]);

    assert (fabs(triple) > 0.0);

    return triple < 0.0;
}


/* ---------------------------------------------------------------------- */
static void
reverse_face_nodes(struct processed_grid *out)
/* ---------------------------------------------------------------------- */
{
    int f, t, *i, *j;

    for (f = 0; f < out->number_of_faces; f++) {
        i = out->face_nodes + (out->face_ptr[f + 0] + 0);
        j = out->face_nodes + (out->face_ptr[f + 1] - 1);

        assert (i <= j);

        while (i < j) {
            t  = *i;
            *i = *j;
            *j =  t;

            i += 1;
            j -= 1;
        }
    }
}


/*-----------------------------------------------------------------
  Public interface
*/
void process_grdecl(const struct grdecl   *in,
                    double                tolerance,
                    struct processed_grid *out)
{
    struct grdecl g;

    size_t i;
    int    sign, error, left_handed;
    int    cellnum;

    int    *actnum, *iptr;
    int    *global_cell_index;

    double *zcorn;

    const size_t BIGNUM = 64;
    const int    nx = in->dims[0];
    const int    ny = in->dims[1];
    const int    nz = in->dims[2];
    const size_t nc = ((size_t) nx) * ((size_t) ny) * ((size_t) nz);

    /* internal work arrays */
    int    *work;
    int    *plist;
    int    *intersections;




    /* -----------------------------------------------------------------*/
    /* Initialize output structure:
       1) allocate space for grid topology (which may need to be
          increased)
       2) set Cartesian imensions
    */
    out->m                = (int) (BIGNUM / 3);
    out->n                = (int) BIGNUM;

    out->face_neighbors   = malloc( BIGNUM      * sizeof *out->face_neighbors);
    out->face_nodes       = malloc( out->n      * sizeof *out->face_nodes);
    out->face_ptr         = malloc((out->m + 1) * sizeof *out->face_ptr);
    out->face_tag         = malloc( out->m      * sizeof *out->face_tag);
    out->face_ptr[0]      = 0;

    out->dimensions[0]    = in->dims[0];
    out->dimensions[1]    = in->dims[1];
    out->dimensions[2]    = in->dims[2];
    out->number_of_faces  = 0;
    out->number_of_nodes  = 0;
    out->number_of_cells  = 0;

    out->node_coordinates = NULL;
    out->local_cell_index = malloc(nc * sizeof *out->local_cell_index);



    /* Do actual work here:*/

    /* -----------------------------------------------------------------*/
    /* For each pillar, compare zcorn values for adjacent cells to
     * find the unique node z-coordinates specified by the input.
     * While here, enumerate unique points and assign point numbers
     * (in plist) for each cornerpoint cell. In other words, plist has
     * 8 node numbers for each cornerpoint cell.*/

    /* initialize grdecl structure "g" that will be processd by
     * "finduniquepoints" */
    g.dims[0] = in->dims[0];
    g.dims[1] = in->dims[1];
    g.dims[2] = in->dims[2];

    actnum    = malloc (nc *     sizeof *actnum);
    g.actnum  = copy_and_permute_actnum(nx, ny, nz, in->actnum, actnum);

    zcorn     = malloc (nc * 8 * sizeof *zcorn);
    sign      = get_zcorn_sign(nx, ny, nz, in->actnum, in->zcorn, &error);
    g.zcorn   = copy_and_permute_zcorn(nx, ny, nz, in->zcorn, sign, zcorn);

    g.coord   = in->coord;


    /* allocate space for cornerpoint numbers plus INT_MIN (INT_MAX)
     * padding */
    plist = malloc(8 * (nc + ((size_t)nx)*((size_t)ny)) * sizeof *plist);

    finduniquepoints(&g, plist, tolerance, out);

    free (zcorn);
    free (actnum);

    /* Determine if coordinate system is left handed or not. */
    left_handed = is_lefthanded(in);
    if (left_handed) {
        /* Reflect Y coordinates about XZ plane to create right-handed
         * coordinate system whilst processing intersections. */
        for (i = 1; i < ((size_t) 3) * out->number_of_nodes; i += 3) {
            out->node_coordinates[i] = -out->node_coordinates[i];
        }
    }

    /* -----------------------------------------------------------------*/
    /* Find face topology and face-to-cell connections */

    /* internal */
    work = malloc(2 * ((size_t) (2*nz + 2)) * sizeof *work);
    for(i = 0; i < ((size_t)4) * (nz + 1); ++i) { work[i] = -1; }

    /* internal array to store intersections */
    intersections = malloc(BIGNUM* sizeof(*intersections));



    process_vertical_faces   (0, &intersections, plist, work, out);
    process_vertical_faces   (1, &intersections, plist, work, out);
    process_horizontal_faces (   &intersections, plist,       out);

    free (plist);
    free (work);

    /* -----------------------------------------------------------------*/
    /* (re)allocate space for and compute coordinates of nodes that
     * arise from intersecting cells (faults) */
    compute_intersection_coordinates(intersections, out);

    free (intersections);

    /* -----------------------------------------------------------------*/
    /* Enumerate compressed cells:
       -make array [0...#cells-1] of global cell numbers
       -make [0...nx*ny*nz-1] array of local cell numbers,
       lexicographically ordered, used to remap out->face_neighbors
    */
    global_cell_index = malloc(nc * sizeof *global_cell_index);
    cellnum = 0;
    for (i = 0; i < nc; ++i) {
        if (out->local_cell_index[i] != -1) {
            global_cell_index[cellnum] = (int) i;
            out->local_cell_index[i]   = cellnum;
            cellnum++;
        }
    }

    /* Remap out->face_neighbors */
    iptr = out->face_neighbors;
    for (i = 0; i < ((size_t) 2) * out->number_of_faces; ++i, ++iptr) {
        if (*iptr != -1){
            *iptr = out->local_cell_index[*iptr];
        }
    }


    free(out->local_cell_index);
    out->local_cell_index = global_cell_index;

    /* Reflect Y coordinate back to original position if left-handed
     * coordinate system was detected and handled earlier. */
    if (left_handed) {
        for (i = 1; i < ((size_t) 3) * out->number_of_nodes; i += 3) {
            out->node_coordinates[i] = -out->node_coordinates[i];
        }
    }

    /* if sign==-1 in ZCORN preprocessing, the sign of the
     * z-coordinate need to change before we finish */
    if (sign == -1)
    {
        for (i = 2; i < ((size_t) 3) * out->number_of_nodes; i += 3)
            out->node_coordinates[i] *= sign;
    }

    /* If an odd number of coordinate reflections were applied, the
     * processing routines--especially facetopology()--will produce
     * node orderings that lead to normals pointing from 2 to 1.
     * Reverse nodes to reestablish expected normal direction (and
     * positive cell volumes). */
    if (left_handed ^ (sign == -1)) {
        reverse_face_nodes(out);
    }
}

/*-------------------------------------------------------*/
void free_processed_grid(struct processed_grid *g)
{
    if( g ){
        free ( g->face_nodes       );
        free ( g->face_ptr         );
        free ( g->face_tag         );
        free ( g->face_neighbors   );
        free ( g->node_coordinates );
        free ( g->local_cell_index );
    }
}

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
