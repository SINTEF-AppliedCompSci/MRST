function [Gt rock2D, petrodata] = getFormationTopGrid(formation,coarsening_level)
% Load the formation grid at a given coarsening level
%
% SYNOPSIS:
%   function [Gt rock2D] = getUtsiraTopGrid(coarsening_level)
%
% DESCRIPTION:
%
% PARAMETERS:
%   coarsening_level - Coarsening factor.  If set to 1, a grid with
%                      approximately one cell per datapoint is produced.
%
% RETURNS:
%   Gt     - top surface grid
%   rock2D - rock structure
%
    moduleCheck('libgeometry');
    [grdecl dataset petroinfo] = ...
        getAtlasGrid(formation, 'coarsening', coarsening_level);

    % Computing the Utsira top-surface grid
    G = processGRDECL(grdecl{1});
    ncvec=nan(numel(G),1);
    for i=1:numel(ncvec)
       ncvec(i)=G(i).cells.num; 
    end
    [nc,j]=max(ncvec);
    G=G(j);
    G = mcomputeGeometry(G); 
    Gt = topSurfaceGrid(G);
    
    % Setting up the rock structure
    petrodata   = petroinfo{1};   
    rock2D.perm = repmat(petrodata.avgperm, Gt.cells.num, 1);
    rock2D.poro = repmat(petrodata.avgporo, Gt.cells.num, 1); 
end
