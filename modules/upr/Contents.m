% Unstructured PEBI-grids for Reservoir (upr) module
%
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  
%
% Routines for creating unstructured PEBI-grids adapting to faults and 
% wells
% 
% AUTHOR: Runar Lie Berge                              (runar.berge@uib.no)
%
%
% ./
%  COPYRIGHT.txt         - Copyright file
%  license.txt           - GNU License file
%  startup.m             - Creates all necessary mappings to run UPR.
%  README                - General readme file
%  [Unstructured 
%  PEBI-grids Adapting 
%  to Geological Features
%  in Subsurface 
%  Reservoir.pdf]        - Master thesis containing detailed description
%                          of UPR
%
% /voronoi2D
%   README2D            - README for files in voronoi2D folder
%   clippedPebi2D       - 2D PEBI-grid clipped against polygon boundary
%   compositePebiGrid2D   - Creates a conforming PEBI-grid embedded in a
%                         Cartesian background grid.
%   surfaceSites2D - Creates a set of Voronoi sites equidistant on
%                           both sides of a set of faults
%   lineSites2D - Creates a set of Voronoi sites tracing a set of
%                          wells
%   CPG2D               - Create a Centroidal Voronoi Diagram
%   surfaceSufCond2D        - Enforce suficient fault condition
%   pebiGrid2D            - Creates a fully unstructured PEBI-grid by using
%                         DistMesh
%   plotLinePath        - Plots a fault or a well path
%   removeConflictPoints - Remove sites that are too close to each other
%   splitAtInt2D          - Split a fault or well paths at the intersections.
%   /examples
%     intersectingWellAndFaults         - one well intersected of several
%                                         faults
%     reservoirWithComplexFaultNetwork  - A reservoir with many faults
%     showOptionValuesCompositePebiGrid - Show the different parameters
%                                         that you can pass to
%                                         compositePebiGrid2D
%     showOPtionValuesPebiGrid          - Show the different parameters you
%                                         can pass to pebiGrid2D
%     statisticalFractures              - Reservoir with many fractures
%     threeIntersectingFaults           - Three faults that intersects
%     wellBranching                     - Intersecting wells
%
% /voronoi3D
%   README3D            - README for files in the voronoi3D directory
%   ballInt             - Calculate the intersection of spheres
%   calcNormals         - Find the normals of each face in a triangulation
%   clipGrid            - Find the intersection of a Voronoi diagram and a
%                         surface in 3D
%   compositePebiGrid3D - Construct a conforming PEBI-grid embedded in a
%                         cartesian background grid
%   surfaceSites3D - Creates a set of points equidistant on both
%                             sides of a triangulation
%   lineSites3D - Creates a set of points tracing a well path
%   CPG3D               - create a 3D Centroidal Voronoi diagram
%   surfaceSufCond3D      - enforce sufficent and necessary fault condition
%   removeSurfaceConflictSites3D - Removes fault points at the intersection of
%                          faults
%   isColinear          - Tests if points lie on a straight line
%   mirroredPebi3D        - Creates a 3D PEBI-grid of a convext hull by
%                         mirroring all points around the boundary faces
%   voronoi2mrstGrid3D        - Transform the output of VORONOIN() to a MRST grid
%                         structure
%   lineSufCond3D       - Enforces the well condition
%   /examples
%       compositePebiGridWithTwoFaults - Create a grid using the composite
%                                        wrapping function of a reservoir
%                                        with two intersecting faults
%       CVDoptimization     - Create a CVD grid of a reservoir with a
%                             single fault
%       saltDome            - Create a salt dome
%       transformQhullGridStructure2MrstGridStructure - Show how to
%                             create a valid MRST grid from the returned
%                             values of voronoin()
%       wellBranching   - A multilateral well
%       
%
% /distMesh:            - Contains a modifed version of distMesh.
% /util:                - Contains some utility functions.
