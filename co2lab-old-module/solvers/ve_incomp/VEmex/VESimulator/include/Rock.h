#ifndef ROCK_H
#define ROCK_H

#include <string.h>

#include <vector>
#include "configure.h"

/**
 *  A rock structure. Same as used in matlab.
 */
template<class Real>
class Rock
{
    public:
        /**
         *  Makes a new rock object.
         *  \param
         *      cells   - The number of cells in the grid.
         *  \param
         *      perm    - The permeability tensor for the rock.
         *  \param
         *      poros   - The porosity for the rock.
         */
        Rock(int cells, double *perm, double *poros);

        /**
         *  Destroy the rock object.
         */
        ~Rock();
    public:
        /// The number of cells in the rock grid.
        int numCells;

        /// The permeability of the rock.
        Real *permeability;

        /// The porosity of the rock.
        Real *porosity;
};

#include "tmpl/Rock.tmpl"

#endif
