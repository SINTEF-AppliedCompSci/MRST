#ifndef STATE_H
#define STATE_H

#include <vector>

#include "configure.h"

/**
 *  A struct describing the current state of the solution.
 */
template<class Real>
class State
{
    public:
        /**
         *  Make a state given values.
         */
        State(int num_cells, int num_faces, double *f=NULL, double *h=NULL,
                double *max_h=NULL);

        /**
         *  Destroy the state.
         */
        ~State();

        void getHeight(double *h, double *max_h);
    public:
        /// The number of cells in the grid.
        int numCells;

        /// The number of faces in the grid.
        int numFaces;

        /// The current flux in the system.
        Real *flux;

        /// The current height of the state.
        Real *height;

        /// The historical maximum height
        Real *max_height;
};

#include "tmpl/State.tmpl"

#endif
