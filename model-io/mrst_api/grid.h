/*
  Copyright 2010, 2011, 2012 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_GRID_HEADER_INCLUDED
#define OPM_GRID_HEADER_INCLUDED


#ifdef __cplusplus
extern "C" {
#endif


struct UnstructuredGrid {
    int    dimensions;

    int    number_of_cells;
    int    number_of_faces;
    int    number_of_nodes;

    int    *face_nodes;
    int    *face_nodepos;
    int    *face_cells;

    int    *cell_faces;
    int    *cell_facepos;

    double *node_coordinates;

    double *face_centroids;
    double *face_areas;
    double *face_normals;

    double *cell_centroids;
    double *cell_volumes;


    int    *global_cell;

    int     cartdims[3];
    int    *cell_facetag;
};

#ifdef __cplusplus
}
#endif

#endif /* OPM_GRID_HEADER_INCLUDED */
