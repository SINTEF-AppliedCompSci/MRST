#ifndef GRID2D_H
#define GRID2D_H

#include <stdlib.h>
#include <vector>
#include <omp.h>
#include <assert.h>
#include "configure.h"

/**
 *  Structure for containing cells. The exact same structure as used in matlab.
 */
template<class Real>
struct cells_t
{
    /// The number of cells.
    int     num;
    /// The number of faces for each cell.
    int     *numFaces;

    /// Index of each face connected to the cells.
    int     *faces;
    /// Height of the column.
    Real    *H;
    /// Index to column that the cell belongs to.
    int     *columnPos;
    /// Normals for each cell.
    Real  *normals;
    /// Volumes for each cell.
    Real  *volumes;
    /// Centroids for each cell in a (x,y) pair.
    Real  *centroids;
    /// z value for each cell.
    Real  *z;
};

/**
 *  Structure for containing faces. The exact same structure as used in matlab.
 */
struct faces_t
{
    /// Number of faces.
    int num;
    /// Indirection map of size num+1 into the nodes array.
    int *neighbors;
    /// Index into nodes for the faces.
    int *nodes;
};

/**
 *  Structure for containing nodes. As simple as it gets, only keeps
 *  coordinates. The exact same structure as used in matlab.
 */
template<class Real>
struct nodes_t
{
    /// Number of nodes.
    int     num;
    /// The (x,y) coordinates for the node.
    Real *coords;
};

/**
 *  Structure for containing columns. The exact same structure as used in
 *  matlab.
 */
template<class Real>
struct columns_t
{
    /// Number of cells.
    int     num;
    /// The index for the cells in each column.
    int     *cells;

    /// Change in z direction.
    Real  *dz;
    /// z value for the column.
    Real  *z;
};

/**
 *  Grid structure. Basicly a wrapper around the raw C structures used.
 */
template<class Real>
class Grid2D
{
    public:
        /**
         *  Make a new grid using the given structures.
         *  \param
         *      cells2d     - The cell structure to copy.
         *  \param
         *      faces2d     - The face structure to copy.
         *  \param
         *      nodes2d     - The node structure to copy.
         *  \param
         *      columns2d   - The column structure to copy.
         *  \param
         *      dim         - The dimensions in  x and y direction.
         *  \param
         *      topFaces    - The topFaces in the 3D grid. Defaults to NULL.
         */
        Grid2D(cells_t<double> *cells2d, faces_t *faces2d,
                nodes_t<double> *nodes2d, columns_t<double> *columns2d,
                int dim[2]);

        //constructor when topsurfacegrid is made in C++
        Grid2D(cells_t<double> *cells2d, faces_t *faces2d,
                nodes_t<double> *nodes2d, columns_t<double> *columns2d,
                int dim[2], bool cpp);

        /**
         *  Destructor.
         */
        ~Grid2D();

        /**
         *  Compute the mobility for a cell.
         *  \param
         *      cell - The index of the cell to compute for.
         *  \param
         *      height - The height of the current cell.
         *  \param
         *      perm - An array of the column permeabilities.
         *  \return
         *      The mobility for cell number cell.
         */
        /*
        Real cellMobility(const int& cell, const Real& height,
                const Real *perm) const;
        */
        void cellRelperm(Real& kr,Real& dkr,const int& cell, const Real& height,
                const Real *perm) const;
    public:
        /// Dimension of grid in (x,y) direction.
        int dims[2];
        /// Pointer to cell structure used.
        cells_t<Real>   *cells;
        /// Pointer to face structure used.
        faces_t         *faces;
        /// Pointer to node structure used.
        nodes_t<Real>   *nodes;
        /// Pointer to column structure used.
        columns_t<Real> *columns;
    protected:
        // True if grid is created with C++ (i.e. not Matlab)
        bool isCpp;
        /**
         *  Copy a cell structure to this structure.
         *  \param
         *      oldCells - The cell structure to copy.
         */
        void copyCells(cells_t<double> *oldCells);

        /**
         *  Copy a face structure to this structure.
         *  \param
         *      oldFaces - The face structure to copy.
         */
        void copyFaces(faces_t *oldFaces);

        /**
         *  Copy a node structure to this structure.
         *  \param
         *      oldNodes - The node structure to copy.
         */
        void copyNodes(nodes_t<double> *oldNodes);

        /**
         *  Copy a column structure to this structure.
         *  \param
         *      oldCols - The column structure to copy.
         */
        void copyColumns(columns_t<double> *oldCols);
};

#include "tmpl/Grid2D.tmpl"
#endif
