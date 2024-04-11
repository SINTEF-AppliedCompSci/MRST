#ifndef WELL_H
#define WELL_H

/**
 *  A structure describing the properties of wells.
 */
template<class Real>
class Well
{
    public:
        /**
         *  Make a new well structure.
         *  \param
         *      num             - The number of wells.
         *  \param
         *      cell_idx        - The indexes of the cells containing a well,
         *      one pr. num.
         *  \param
         *      w_type          - The type of well. One for each well.
         *  \param
         *      rate            - The rate of injection/extraction for each
         *      well.
         *  \param
         *      rad             - The radius of the wells.
         *  \param
         *      direction       - The direction of each well.
         *  \param
         *      productivity    - Well productivity.
         *  \param
         *      displacement    - The displacement of each well in the
         *      z-direction.
         *  \param
         *      w_name          - The name of each well.
         *  \param
         *      composit        - Fluid composition, only used for injectors.
         *  \param
         *      depth           - The depth of each well.
         *  \param
         *      sig             - The sign of the wells, either positive or
         *      negative.
         *  \param
         *      h               - The height of each well.
         */
        Well(int num=0, int *cell_idx=NULL, double *rate=NULL, double *h=NULL);

        /*
         *  Copy a well.
         *  \param
         *      w - The new well.
         */
        Well(Well<double> *w);

        /**
         *  Destroy the well.
         */
        ~Well();

        /**
         *  Compute contribution from the wells to the q vector.
         *  \param
         *      q - The sources/sink vector.
         *  \param
         *      H - The height of the reservoir for each cell.
         */
        void contribute_source(Real *q, const Real *H) const;
    public:
        /// The number of cells that contains a well.
        int                 numCells;

        /// Indexes of well cells.
        int *cells;

        /// Value of injection rate. One for each well.
        Real *val;

        /// Height of the well.
        Real *height;
};

#include "tmpl/Well.tmpl"

#endif
