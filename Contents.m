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
%   compositePebiGrid   - Creates a conforming PEBI-grid embedded in a
%                         Cartesian background grid.
%   pebiGrid            - Creates a fully unstructured PEBI-grid by using
%                         DistMesh
%   createFaultGridPoints - Creates a set of Voronoi sites equidistant on
%                           both sides of a fault
%   createWellGridPoints - Creates a set of Voronoi sites tracing a well
%   plotLinePath         - Plots a fault or a well path
%   README              - Readme file for voronoi2D
%   removeConflictPoints2 - Remove sites that are too close to each other
%   splitAtInt          - Split a fault or well paths at the intersections.
%   /examples
%
% /distMesh:            - Contains a modifed version of distMesh.
%  
%
