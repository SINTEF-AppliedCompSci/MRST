% Unstructured PEBI-grids for Reservoir (upr) module
%
% Routines for creating unstructured PEBI-grids adapting to faults and wells
% Examples in /examples
% 
% AUTHOR: runarlb@stud.ntnu.no
%
%
% ./
%  startup.m             - Creates all necessary mappings to run UPR.
%  README                - General readme file
%  COPYRIGHT.txt         - Copyright file
%
% /voronoi2D
%   README              - README for files in voronoi2D folder
%   clippedPebi2D       - 2D PEBI-grid clipped against polygon boundary
%   createCVD           - Centroidal Voronoi Diagram
%   faultSufCond        - Enforce suficient fault condition
%   compositePebiGrid   - Creates a conforming PEBI-grid embedded in a
%                         Cartesian background grid.
%   pebiGrid            - Creates a fully unstructured PEBI-grid by using
%                         DistMesh
%   createFaultGridPoints - Creates a set of Voronoi sites equidistant on
%                           both sides of a fault
%   createWellGridPoints - Creates a set of Voronoi sites tracing a well
%   plotLinePath         - Plots a fault or a well path
%   README2d             - Readme file for voronoi2D
%   removeConflictPoints2 - Remove sites that are too close to each other
%   splitAtInt          - Split a fault or well paths at the intersections.
%   /examples
%     intersectingWellAndFaults         - one well intersected of several
%                                         faults
%     reservoirWithComplexFaultNetwork  - A reservoir with many faults
%     showOptionValuesCompositePebiGrid - Show the different parameters
%                                         that you can pass to
%                                         compositePebiGrid
%     showOPtionValuesPebiGrid          - Show the different parameters you
%                                         can pass to pebiGrid
%     statisticalFractures              - Reservoir with many fractures
%     threeIntersectingFaults           - Thre faults that intersects
%     wellBranching                     - Intersecting wells
%
%
% /distMesh:            - Contains a modifed version of distMesh.
% /util:                - Contains some utility functions.  
%
