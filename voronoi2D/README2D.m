%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This matlab code was written as a part of my master thesis. It was
% written during the fall of 2015 and spring of 2016. My masters thesis is
% written in cooperation with SINTEF ICT.
% Runar Lie Berge                                      (runar.berge@uib.no)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Voronoi-2D contains several functions for generating unstructured Voronoi
% grids (or Pebi grids) in MRST (http://www.sintef.no/MRST). The main focus
% of this software is the construction of Pebi grids conforming to faults, 
% fractures or other geological structures. 
% 
% Voronoi-2D supports two different types of structures.
% For the first structure one wants cell edges to follow the
% structure. For the second type, one wants cell centroids to follow
% the structure. In Voronoi2D the first type of structure is
% exclusively called for faults, and the second type for wells.
% 
% 
% Voronoi-2D contains the following functions
% Functions:
%     clippedPebi2D
%     compositePebiGrid
%     createFaultGridPoints
%     createWellGridPoints
%     CVD2D
%     faultSufCond
%     pebiGrid
%     plotLinePath
%     removeConflicPoints2
%     splitAtInt
% 
% - clippedPebiGrid2D creates a clipped Voronoi diagram from a point set
%   and a polygon boundary
%
% - compositePebiGrid is an interface function that can be called to create 
%   a valid MRST grid structure. It creates a Pebi grid conforming to 
%   faults and wells. The functions creates a semi-structured grid, by 
%   inserting voronoi seeds around wells and fractures.
% 
% - createFaultGridPoints creates points equiv-distant on each side of the
%   given faults. 
% 
% - createWellGridPoints places points along given well-lines.
%
% - CVD2D creates a centroidal Voronoi diagram by minimizing the CVD energy
%   function. It does this by using the L-BFGS optimization algorithm.
% 
% - faultSufCond enforces the sufficient fault condition.
%
% - pebiGrid is an interface function that can be called to create a valid 
%   MRST grid structure. It creates a fully unstructured PEBI-grid that 
%   conforms to faults and wells. It uses the software DistMesh to create 
%   the background grid. DistMesh is a software for creating unstructured 
%   delaunay triangulations: Per-Olof Persson and Gilbert Strang, "A Simple
%   Mesh Generator in MATLAB," SIAM Review Vol. 46 (2) 2004.
% 
% - plotLinePath plots the line paths that are stored in a cell array
%
% - removeConflicPoints2 removes any points from one set that are too close 
%   to points from another set. This function can be used to remove small 
%   orconstricted cells.
% 
% - splitAtInt is a function that splits a set of paths at each 
%   intersection. It can be used to split all faults and wells at their 
%   intersections.
% 
% 
% EXAMPLES: 
% For examples see the subfolder voronoi2d/examples