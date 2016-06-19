%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This matlab code was written as a part of my master thesis during the 
fall of 2015 and spring of 2016. My masters thesis is written in 
cooperation with SINTEF ICT.
Runar Lie Berge                                      (runarlb@stud.ntnu.no)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

UPR contains two subfolders with functions. voronoi2D contains functions
for generate a grid adapting to surfaces in 2D, while voronoi3D contains 
functions for generate a grid adapting to surfaces in 3D. Please see the
respecite folders README file for more information on the specific
functions.


EXAMPLES: 
For examples see the subfolder examples/


We have included distMesh in this module. DistMesh is a software for creating unstructured delaunay  triangulations: Per-Olof Persson and Gilbert Strang, "A Simple Mesh Generator in MATLAB," SIAM Review Vol. 46 (2) 2004.
We have made some modifications to DistMesh:
distmesh2d.m:
  Removed plotting, added maximum number of iterations and fixed initial grid
  size from relative to absolute. 

  Line 1:  function [p,t]= ... -> function [p,t, IC]= ...
  Line 60: densityctrlfreq=30; -> densityctrlfreq=30; maxIt = 1000;
  Line 70: p=p(rand(size(p,1),1)<r0./max(r0),:); -> p=p(rand(size(p,1),1)<r0,:);
  Line 78: removed
  Line 79: while 1 -> while count<maxIt
  Line 91-92: removed
  Line 124: added; if count == maxIt
                       warning('DistMesh did not converge in maximum number of iterations.')
                   end
 Line 126: [p,t]=fixmesh(p,t); -> [p,t, ~, sorting]=fixmesh(p,t);
 Line 127: removed

fixmesh:
  Added a vector which keeps track of the how fixmesh changes the ordering
  of p. 
  line  1: function [p,t,pix]= ... -> function [p,t,pix,sorting]= ...
  line  9: added; sorting = (1:size(p,1))';
  line 13: added; sorting = sorting(ix,:);
  line 19.5: added sorting = sorting(pix);
