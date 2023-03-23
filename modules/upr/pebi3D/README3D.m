%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This matlab code was written as a part of my master thesis. It was
% written during the fall of 2015 and spring of 2016. My masters thesis is
% written in cooperation with SINTEF ICT.
% Runar Lie Berge                                    (runarlb@stud.ntnu.no)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% upr contains several functions for generating unstructured Voronoi
% grids (or Pebi grids) in MRST (http://www.sintef.no/MRST). The main focus 
% of this software is the construction of Pebi grids conforming to face
% constraints (e.g., faults or fractures), or cell constraints (e.g., wells)
% 
% upr supports two different types of structures.
% For the first structure one wants cell faces to follow the
% structure (a face constraint). For the second type, one wants cell centroids
% to follow the structure (a cell constraint).
% The pebi3D folder contains the following functions
%
% Functions:
%   ballInt
%   calcNormals
%   compositePebiGrid3D
%   CPG3D
%   ellipticSurface3D
%   lineGrid3D
%   lineSufCond3D
%   mirroredPebi3D
%   removeSurfaceConflictSites3D
%   surfaceGrid3D
%   surfaceIntersections3D
%   surfaceSites3D
%   surfaceSufCond3D
%   volumeGrid3D
%   voronoiGrid3D
%   voronoi2mrstGrid3D
%
%
% - ballInt calculates the two unique intersection points of three balls
% - 
% - calcNormals calculate the normal vectors of triangles.
%
% - compositePebiGrid3D is an interface function that can be called to 
%   create a valid MRST grid structure. It creates a Pebi grid conforming
%   to face constraints and cell constraints. The functions creates a
%   semi-structured grid with a structured background grid and local refinement
%   around the constraints
%
% - CPG3D creates a centroidal Pebi grid diagram by minimizing the CPG energy
%   function. It does this by using the L-BFGS optimization algorithm.
%
% - ellipticSurface3D approximates an ellipsoid by a polygon. Can be used
%   with surfaceIntersections3D
%
% - lineGrid3D create the 1D grids along all surface intersections as
%   returned from surfaceIntersections3D
%
% - lineSites3D places points along given well-lines.
%
% - lineSufCond3D enforces the condition on a set of sites such that the
%   cells of the resulting Pebi grid traces the cell constraint
%
% - mirroredPebi3D creates a PEBI-grid with a valid MRST grid structure
%   innside a convex domain.
%
% - removeSurfaceConflictSites3D: Given the importance of different surface
%   constraints, this function removes the conflict sites of surfaces with
%   lower importance if they ruin the conformity of surfaces with higher
%   importance.
%
% - surfaceGrid3D creates 2D PEBI-grids of a given set of surfaces, 1D
%   grids of surface intersections and desired mesh sizes.
%
% - surfaceIntersections3D calculate the intersections of surfaces given as
%   polygons.
%
% - surceSites3D places a set of sites on each side of the surface such
%   that the faces of the 3D grid conforms to the surface.
%
% - surfaceSufCond2D: given the cell centers and radii of each 2D vertex,
%   this function enforces the sufficient fault condition on a given set
%   of sites.
%
% - surfaceSufCondFromGrid3D: given the 2D grids of surfaces and a mesh
%   size, this function removes the 3D sites that violate the conformity
%   requirement of the surface.
%
% - volumeGrid3D creates a 3D grid that conform to a set of 2D surface grids.
%
% - voronoi2mrstGrid3D converts a Qhull grid structure into a MRST grid
%   structure.
