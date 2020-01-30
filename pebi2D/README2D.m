%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This matlab code was written as a part of my master thesis. It was
% written during the fall of 2015 and spring of 2016. My masters thesis is
% written in cooperation with SINTEF ICT.
% Runar Lie Berge                                      (runar.berge@uib.no)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Be aware that the code went through a major renaming in Januar 2020!
%
% upr contains several functions for generating unstructured Voronoi
% grids (or Pebi grids) in MRST (http://www.sintef.no/MRST). The main focus
% of this software is the construction of Pebi grids conforming to faults, 
% fractures or other geological structures. 
% 
% upr supports two different types of structures.
% For the first structure one wants cell edges to follow the
% structure. For the second type, one wants cell centroids to follow
% the structure. 
% 
% The subfolder pebi2D contains the following functions
% Functions:
%     clippedPebi2D
%     compositePebiGrid2D
%     surfaceSites2D
%     lineSites2D
%     CPG2D
%     surfaceSufCond2D
%     pebiGrid2D
%     plotLinePath
%     removeConflicPoints2
%     splitAtInt2D
% 
% - clippedPebiGrid2D creates a clipped Pebi grid from a point set
%   and a polygon boundary
%
% - compositePebiGrid2D is an interface function that can be called to create 
%   a valid MRST grid structure. It creates a Pebi grid conforming to 
%   face constraints and cell constraints. The functions creates a
%   semi-structured grid, by inserting voronoi seeds around wells and fractures.
% 
% - surfaceSites2D creates sites equiv-distant on each side of the
%   given surfaces. 
% 
% - lineSites2D places sites along given lines.
%
% - CPG2D creates a centroidal Pebi grid by minimizing the CPG energy
%   function. It does this by using the L-BFGS optimization algorithm.
% 
% - surfaceSufCond2D enforces the sufficient surface condition.
%
% - pebiGrid2D is an interface function that can be called to create a valid 
%   MRST grid structure. It creates a fully unstructured PEBI-grid that 
%   conforms to face constraints and cell constraints. It uses the software DistMesh to create 
%   the background grid. DistMesh is a software for creating unstructured 
%   delaunay triangulations: Per-Olof Persson and Gilbert Strang, "A Simple
%   Mesh Generator in MATLAB," SIAM Review Vol. 46 (2) 2004.
% 
% - plotLinePath plots the line paths that are stored in a cell array
% 
% - splitAtInt2D is a function that splits a set of paths at each 
%   intersection. It can be used to split all face constraints and cell
%   constraints at their intersections.
