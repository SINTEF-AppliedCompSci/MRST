%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This matlab code was written as a part of my master thesis. It was
% written during the fall of 2015 and spring of 2016. My masters thesis is
% written in cooperation with SINTEF ICT.
% Runar Lie Berge                                    (runarlb@stud.ntnu.no)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Voronoi-3D contains several functions for generating unstructured Voronoi
% grids (or Pebi grids) in MRST (http://www.sintef.no/MRST). The main focus 
% of this software is the construction of Pebi grids conforming to faults, 
% fractures or other geological structures. 
% 
% Voronoi-3D supports two different types of structures.
% For the first structure one wants cell faces to follow the
% structure. For the second type, one wants cell centroids to follow
% the structure. In Voronoi3D the first type of structure is
% exclusively called for faults, and the second type for wells.
% Voronoi-3D contains the following functions
%
% Functions:
%   ballInt
%   calcNormals
%   clipGrid
%   clippedPebi3D
%   compositePebiGrid3D
%   createFaultGridPoints3D
%   createWellGridPoints3D
%   CVD3D
%   faultSufCond3D
%   fixFaultIntersection
%   isColinear
%   mirroredPebi
%   voronoi2mrst
%   wellSufCond3D
%
%
% - ballInt calculates the two unique intersection points of three balls
% - 
% - calcNormals calculate the normal vectors of triangles.
%
% - clipGrid calculates the intersection of a Voronoi diagram and a surface
%
% - clippedPebi3D generate a 3D Voronoi diagram that is restricted to a
%   bounded domain.
%
% - compositePebiGrid3D is an interface function that can be called to 
%   create a valid MRST grid structure. It creates a Pebi grid conforming
%   to faults and wells. The functions creates a semi-structured grid, by 
%   inserting voronoi seeds around wells and fractures.
%
% - createWellGridPoints3D places points along given well-lines.
%
% - CVD3D creates a centroidal Voronoi diagram by minimizing the CVD energy
%   function. It does this by using the L-BFGS optimization algorithm.
%
% - faultSufCond enforces the sufficient fault condition.
%
% - fixFaultIntersection removes points from faults that are intersecting
%
% - isColinear tests if a set of points lie on a straight line
%
% - mirroredPebi creates a PEBI-grid with a valid MRST grid structure
%   innside a convex domain.
%
% - voronoi2mrst converts a Qhull grid structure into a MRST grid
%   structure.
%
% - wellSufCond3D enforces the sufficient well condition.
%
% EXAMPLES:
% see subfolder voronoi3D/examples