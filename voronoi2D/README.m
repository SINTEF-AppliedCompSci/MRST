%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This matlab code was written as a part of my master thesis during the 
% fall of 2015 and spring of 2016. My masters thesis is written in 
% cooperation with SINTEF ICT.
% Runar Lie Berge                                      (runarlb@stud.ntnu.no)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Voronoi-2D contains several functions for generating unstructured Voronoi
% grids (or Pebi grids) in MRST (http://www.sintef.no/MRST). The main focus 
% of this software is the construction of Pebi grids conforming to faults, 
% fractures or other geological structures. 
% 
% Voronoi-2D supports two different types of structures.
% For the first structure type one wants cell edges to follow the 
% structure. For the second type, one wants cell centroids to follow 
% the structure. In Voronoi-2D is the first type of structure is exclusively 
% called for faults, and the second type for wells. 
% 
% 
% Voronoi-2D contains the following functions
% Functions:
%     compositePebiGrid
%     createFaultGridPoints
%     createWellGridPoints
%     pebiGrid
%     removeConflicPoints2
%     splitAtInt
% 
% - compositePebiGrid is an interface function that can be called to create 
% a valid MRST grid structure. It creates a Pebi grid conforming to faults 
% and wells. The functions creates a semi-structured grid, by inserting 
% voronoi seeds around wells and fractures.
% 
% - createFaultGridPoints creates points equiv-distant on each side of the
% given faults. 
% 
% - createWellGridPoints places points along given well-lines.
% 
% - pebiGrid is an interface function that can be called to create a valid 
% MRST grid structure. It creates a fully unstructured pebi grid that 
% conforms to faults and wells. It uses the software DistMesh to create the 
% background grid. DistMesh is a software for creating unstructured delaunay 
% triangulations: Per-Olof Persson and Gilbert Strang, "A Simple Mesh 
% Generator in MATLAB," SIAM Review Vol. 46 (2) 2004.
% 
% - removeConflicPoints2 removes any points from one set that are too close 
% to points from another set. This function can be used to remove small or
% constricted cells.
% 
% - splitAtInt is a function that splits a set of paths at each intersection.
% It can be used to split all faults and wells at their intersections.
% 
% 
% EXAMPLES: 
% For examples see the subfolder voronoi2d/examples
% 
% 
% Changes in distMesh:
% distmesh2d.m:
%   Removed plotting, added maximum number of iterations and fixed initial grid
%   size from relative to absolute. 
% 
%   Line 1:  function [p,t]= ... -> function [p,t, IC]= ...
%   Line 60: densityctrlfreq=30; -> densityctrlfreq=30; maxIt = 1000;
%   Line 70: p=p(rand(size(p,1),1)<r0./max(r0),:); -> p=p(rand(size(p,1),1)<r0,:);
%   Line 78: removed
%   Line 79: while 1 -> while count<maxIt
%   Line 91-92: removed
%   Line 124: added; if count == maxIt
%                        warning('DistMesh did not converge in maximum number of iterations.')
%                    end
%  Line 126: [p,t]=fixmesh(p,t); -> [p,t, ~, sorting]=fixmesh(p,t);
%  Line 127: removed
% 
% fixmesh:
%   Added a vector which keeps track of the how fixmesh changes the ordering
%   of p. 
%   line  1: function [p,t,pix]= ... -> function [p,t,pix,sorting]= ...
%   line  9: added; sorting = (1:size(p,1))';
%   line 13: added; sorting = sorting(ix,:);
%   line 19.5: added sorting = sorting(pix);
