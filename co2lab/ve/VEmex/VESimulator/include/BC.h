#ifndef BC_H
#define BC_H

#include <vector>
#include <string>

#include "configure.h"

/**
 *  Structure for holding boundary conditions.
 */
template<class Real>
class BC
{
    public:
        /**
         *  Make a new set of boundary conditions.
         *  \param
         *      num_faces   - The number of boundary condition faces. Defaults
         *      to 0 (No boundary conditions).
         *  \param
         *      bc_face     - The index of the faces in the faces struct.
         *      Should be num_faces of them. Defaults to NULL.
         *  \param
         *      bc_values   - Values for the boundary conditions. Should be
         *      num_faces of them. Defaults to NULL.
         *  \param
         *      bc_height   - The height of the node in the boundary
         *      conditions. Should be num_faces of them. Default to NULL.
         *  \param
         *      bc_type     - The type of boundary condition we have. One for
         *      each face.
         */
        BC(int num_faces=0, int *bc_face=NULL, double *bc_values=NULL,
                double *bc_height=NULL, char **bc_type=NULL);

        /**
         *  Copy a boundary condition.
         *  \param
         *      bc - A boundary condition with double precision data.
         */
        BC(BC<double> *bc);

        /**
         *  Destructor.
         */
        ~BC();
    public:
        /// Number of boundary condition faces.
        int                 num;

        /// Face indexes for the boundaries.
        int *face;

        /// Value for each boundary face.
        Real *value;

        /// Height for each face.
        Real *height;

        /// Type of boundary condition for each face.
        std::string *type;
};

#include "tmpl/BC.tmpl"

#endif
