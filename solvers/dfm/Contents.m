% Discrete Fracture Matrix (DFM) module
%
% Routines supporting the DFM method. Examples in /examples
% Add module by typing MrstModule add dfm
% AUTHOR: tor.harald.sandve@iris.no
%
% new files
%   addhybrid                   - Add hybrid cells to the grid structure.
%   build_fractures_mod         - Build fracture network that serves as constraints in the triangulation.
%   computeHybridTrans          - Computes the hybrid-hybrid transmissibilities between a hybrid cell 1 and 2
%   nodeType                    - return node type
%   plotEdges                   - plot lines
%   testNormals                 - Tests if the normals points from neighbor 1 to 2 and returns the index of
%   triangulate                 - Make a Delaunay Triangulation and return a Mrst grid
%
%   /plotting/plotFractures - Plots 2d hybridcells.
%
% Files modified from  core MRST functions. Use these to
%   computeTrans_DFM            - Compute transmissibilities using a two-point scheme.
%   computeMultiPointTrans_DFM  - Compute multi-point transmissibilities.
%   twophaseJacobian_DFM        - Residual and Jacobian of single point upwind solver for two-phase flow.
%   explicitTransport_DFM       - Explicit single point upwind transport solver for two-phase flow.
%   implicitTransport_DFM       - Implicit single point upwind transport solver for two-phase flow.
%   removeInternalBoundary_DFM  - Remove internal boundary in grid by merging faces in face list N
%
%   /plotting/                 - Modification to plotting rutines to handle
%
% /private.
%  - Files needed since we can not access private folders from this location
%
% /examples:
%   add_point                   - Add a new point to the cloud.
%   computeHybridMPTrans        - Computes the hybrid-hybrid transmissibilities. Also update the grid and
%   computeTimeOfFlight_DFM     - Compute time of flight using finite-volume scheme.
%   distance_to_closest_line    - Copyright 2011-2012 University of Bergen
%   incompMPFA_DFM              - Solve incompressible flow problem (fluxes/pressures) using MPFA-O method.
%   incompTPFA_DFM              - Solve incompressible flow problem (fluxes/pressures) using TPFA method.
%   lines_intersect             - [pt] = lines_intersect (vertices, edge1, edge2)
%   partition_edges             - [vertices, new_edges, num_added] = ...
%   read_openoffice             - Read fractures drawn in Libreoffice, and represent them as points and
%   removeFractureIntersections - Partition intersecting edges by adding a point in the intersection.
%   remove_closepoints          - Remove all points that are too close to the lines, since these will only
%   snap_to_grid                - Move vertices to the closest (structured) grid point.
%   split_edge                  - Split an edge into two parts
