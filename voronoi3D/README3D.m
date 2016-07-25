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
% 
% Voronoi-3D contains the following functions
% Functions:
%   calcNormals
%   clipGrid
%   clippedPebi3D
%   isColinear
%   voronoi2mrst
%
% 
% - calcNormals calculate the normal vectors of triangles.
%
% - clipGrid calculates the intersection of a Voronoi diagram and a surface
%
% - clippedPebi3D generate a 3D Voronoi diagram that is restricted to a
%   bounded domain.
% 
% - isColinear tests if a set of points lie on a straight line
%
% - voronoi2mrst converts a Qhull grid structure into a MRST grid
% structure.