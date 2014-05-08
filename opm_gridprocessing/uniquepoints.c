/*===========================================================================
//
// File: uniquepoints.c
//
// Created: Fri Jun 19 08:46:05 2009
//
// Author: Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
//
// $Date: 2012-06-19 15:22:05 +0200 (Tue, 19 Jun 2012) $
//
// $Revision: 983 $
//
//==========================================================================*/

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

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
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "preprocess.h"
#include "uniquepoints.h"

#define MIN(i,j) (((i) < (j)) ? (i) : (j))
#define MAX(i,j) (((i) > (j)) ? (i) : (j))

/*-----------------------------------------------------------------
  Compare function passed to qsortx  */
static int compare(const void *a, const void *b)
{
    const double a0 = *(const double*) a;
    const double b0 = *(const double*) b;

    /*                { -1, a < b
     * compare(a,b) = {  0, a = b
     *                {  1, a > b */
    return (a0 > b0) - (a0 < b0);
}

/*-----------------------------------------------------------------
  Creat sorted list of z-values in zcorn with actnum==1x */
static int createSortedList(double *list, int n, int m,
                            const double *z[], const int *a[])
{
    int i,j;
    double *ptr = list;
    for (i=0; i<n; ++i){
        for (j=0; j<m; ++j){
            if (a[j][i/2])  *ptr++ = z[j][i];
            /* else        fprintf(stderr, "skipping point in inactive cell\n"); */
        }
    }

    qsort(list, ptr-list, sizeof(double), compare);
    return ptr-list;
}


/*-----------------------------------------------------------------
  Remove points less than <tolerance> apart in <list> of increasing
  doubles.  */
static int uniquify(int n, double *list, double tolerance)
{
    int    i;
    int    pos;
    double val;

    assert (!(tolerance < 0.0));

    if (n<1) return 0;
    pos = 0;
    val = list[pos++];/* Keep first value */

    for (i=1; i<n; ++i){
        if (val + tolerance < list [i]){
            val         = list[i];
            list[pos++] = val;
        }
    }

    /*
      Preserve outer z-boundary.

      This operation is a no-op in the case

      list[pos-2] + tolerance < list[n-1].

      If, however, the second to last point is less than <tolerance>
      away from the last point (list[n-1]), we remove this
      second-to-last point as it cannot be distinguished from "final"
      point.
    */
    list[pos-1] = list[n-1];

    return pos;
}


/*-----------------------------------------------------------------
  Along single pillar: */
static int assignPointNumbers(int    begin,
                              int    end,
                              const double *zlist,
                              int    n,
                              const double *zcorn,
                              const int    *actnum,
                              int    *plist,
                              double tolerance)
{
    /* n     - number of cells */
    /* zlist - list of len unique z-values */
    /* start - number of unique z-values processed before. */

    int i, k;
    /* All points should now be within tolerance of a listed point. */


    const double *z = zcorn;
    const int    *a = actnum;
    int    *p = plist;

    k = begin;
    *p++ = INT_MIN; /* Padding to ease processing of faults */
    for (i=0; i<n; ++i){

        /* Skip inactive cells */
        if (!a[i/2]) {
            p[0] = p[-1];  /* Inactive cells are collapsed leaving
                            * void space.*/
            ++p;
            continue;
        }

        /* Find next k such that zlist[k] < z[i] < zlist[k+1] */
        while ((k < end) && (zlist[k] + tolerance < z[i])){
            k++;
        }

        /* assert (k < len && z[i] - zlist[k] <= tolerance) */
        if ((k == end) || ( zlist[k] + tolerance < z[i])){
            fprintf(stderr, "Cannot associate  zcorn values with given list\n");
            fprintf(stderr, "of z-coordinates to given tolerance\n");
            return 0;
        }

        *p++ = k;
    }
    *p++ = INT_MAX;/* Padding to ease processing of faults */


    return 1;
}


/* ---------------------------------------------------------------------- */
static void
vector_positions(const int dims[3] ,
                 const int i       ,
                 const int j       ,
                 size_t    start[4])
/* ---------------------------------------------------------------------- */
{
    size_t im, ip, jm, jp;

    im = MAX(1,       i  ) - 1;
    jm = MAX(1,       j  ) - 1;
    ip = MIN(dims[0], i+1) - 1;
    jp = MIN(dims[1], j+1) - 1;

    start[ 0 ] = dims[2] * (im + dims[0]*jm);
    start[ 1 ] = dims[2] * (im + dims[0]*jp);
    start[ 2 ] = dims[2] * (ip + dims[0]*jm);
    start[ 3 ] = dims[2] * (ip + dims[0]*jp);
}


/*-----------------------------------------------------------------
  Given a vector <field> with k index running faster than i running
  faster than j, and Cartesian dimensions <dims>, find pointers to the
  (i-1, j-1, 0), (i-1, j, 0), (i, j-1, 0) and (i, j, 0) elements of
  field.  */
static void igetvectors(const int dims[3], int i, int j,
                        const int *field, const int *v[])
{
    size_t p, start[4];

    vector_positions(dims, i, j, start);

    for (p = 0; p < 4; p++) {
        v[p] = field + start[p];
    }
}


/*-----------------------------------------------------------------
  Given a vector <field> with k index running faster than i running
  faster than j, and Cartesian dimensions <dims>, find pointers to the
  (i-1, j-1, 0), (i-1, j, 0), (i, j-1, 0) and (i, j, 0) elements of
  field.  */
static void dgetvectors(const int dims[3], int i, int j,
                        const double *field, const double *v[])
{
    size_t p, start[4];

    vector_positions(dims, i, j, start);

    for (p = 0; p < 4; p++) {
        v[p] = field + start[p];
    }
}


/*-----------------------------------------------------------------
  Given a z coordinate, find x and y coordinates on line defined by
  coord.  Coord points to a vector of 6 doubles [x0,y0,z0,x1,y1,z1].
  */
static void interpolate_pillar(const double *coord, double *pt)
{
    double a;

    if (fabs(coord[5] - coord[2]) > 0) {

        a = (pt[2] - coord[2]) / (coord[5] - coord[2]);

    } else {

        a = 0;

    }

#if 0
    pt[0]       = coord[0] + a*(coord[3]-coord[0]);
    pt[1]       = coord[1] + a*(coord[4]-coord[1]);
#else
    pt[0]       = (1.0 - a)*coord[0] + a*coord[3];
    pt[1]       = (1.0 - a)*coord[1] + a*coord[4];
#endif
}

/*-----------------------------------------------------------------
  Assign point numbers p such that "zlist(p)==zcorn".  Assume that
  coordinate number is arranged in a sequence such that the natural
  index is (k,i,j) */
int finduniquepoints(const struct grdecl *g,
                     /* return values: */
                     int           *plist, /* list of point numbers on
                                            * each pillar*/
                     double tolerance,
                     struct processed_grid *out)

{

    const int nx = out->dimensions[0];
    const int ny = out->dimensions[1];
    const int nz = out->dimensions[2];
    const int nc = g->dims[0]*g->dims[1]*g->dims[2];


    /* zlist may need extra space temporarily due to simple boundary
     * treatement  */
    int            npillarpoints = 8*(nx+1)*(ny+1)*nz;
    int            npillars      = (nx+1)*(ny+1);

    double *zlist = malloc(npillarpoints*sizeof *zlist);
    int     *zptr = malloc((npillars+1)*sizeof *zptr);




    int     i,j,k;

    int     d1[3];
    int     len    = 0;
    double  *zout  = zlist;
    int     pos    = 0;
    double *pt;
    const double *z[4];
    const int *a[4];
    int *p;
    int pix, cix;
    int zix;

    const double *coord = g->coord;

    d1[0] = 2*g->dims[0];
    d1[1] = 2*g->dims[1];
    d1[2] = 2*g->dims[2];

    out->node_coordinates = malloc (3*8*nc*sizeof(*out->node_coordinates));

    zptr[pos++] = zout - zlist;

    pt    = out->node_coordinates;

    /* Loop over pillars, find unique points on each pillar */
    for (j=0; j < g->dims[1]+1; ++j){
        for (i=0; i < g->dims[0]+1; ++i){

            /* Get positioned pointers for actnum and zcorn data */
            igetvectors(g->dims,   i,   j, g->actnum, a);
            dgetvectors(d1,      2*i, 2*j, g->zcorn,  z);

            len = createSortedList(     zout, d1[2], 4, z, a);
            len = uniquify        (len, zout, tolerance);

            /* Assign unique points */
            for (k=0; k<len; ++k){
                pt[2] = zout[k];
                interpolate_pillar(coord, pt);
                pt += 3;
            }

            /* Increment pointer to sparse table of unique zcorn
             * values */
            zout        = zout + len;
            zptr[pos++] = zout - zlist;

            coord += 6;
        }
    }
    out->number_of_nodes_on_pillars = zptr[pos-1];
    out->number_of_nodes            = zptr[pos-1];

    /* Loop over all vertical sets of zcorn values, assign point
     * numbers */
    p = plist;
    for (j=0; j < 2*g->dims[1]; ++j){
        for (i=0; i < 2*g->dims[0]; ++i){

            /* pillar index */
            pix = (i+1)/2 + (g->dims[0]+1)*((j+1)/2);

            /* cell column position */
            cix = g->dims[2]*((i/2) + (j/2)*g->dims[0]);

            /* zcorn column position */
            zix = 2*g->dims[2]*(i+2*g->dims[0]*j);

            if (!assignPointNumbers(zptr[pix], zptr[pix+1], zlist,
                                    2*g->dims[2],
                                    g->zcorn  + zix, g->actnum + cix,
                                    p, tolerance)){
                fprintf(stderr, "Something went wrong in assignPointNumbers");
                return 0;
            }

            p += 2 + 2*g->dims[2];
        }
    }

    free(zptr);
    free(zlist);

    return 1;
}

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
