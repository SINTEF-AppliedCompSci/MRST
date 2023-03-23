%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This matlab code was written as a part of my master thesis during the 
% fall of 2015 and spring of 2016. My masters thesis is written in 
% cooperation with SINTEF ICT.
% Runar Lie Berge                                   (rube@norceresearch.no)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2015-2020   Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
%
% UPR contains four subfolders with functions. 
%
% pebi2D/ contains functions for generate a grid adapting to surfaces in
% 2D. Please see the README file in this folder for more information. 
%
% pebi3D contains functions for generating PEBI grids in 3D. Please
% see readme file in this subfolder for more information.
%
% The subfolder util/ contains some utility functions used by the functions
% in  pebi2D/ and pebi3D/
%
% distmesh/ contains a modified version of the software distmesh: 
% Per-Olof Persson and Gilbert Strang, "A Simple Mesh Generator in MATLAB,"
% SIAM Review Vol. 46 (2) 2004.
%
%
% If you use upr, please cite the paper:
% Berge, R.L., Klemetsdal, Ø.S. & Lie, KA. Comput Geosci (2019) 23: 169. https://doi.org/10.1007/s10596-018-9790-0