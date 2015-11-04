function [Gt rock2D] = getUtsiraTopGrid(coarsening_level, do_cut)
% Load the Utsira grid at a given coarsening level, and optional cutting-away
% of Eastern part. 
%
% SYNOPSIS:
%   function [Gt rock2D] = getUtsiraTopGrid(coarsening_level, do_cut)
%
% DESCRIPTION:
%
% PARAMETERS:
%   coarsening_level - Coarsening factor.  If set to 1, a grid with
%                      approximately one cell per datapoint is produced.
%   do_cut           - If 'true', the part of the grid east of coordinate
%                      6575000 will be cut away.
%
% RETURNS:
%   Gt     - top surface grid
%   rock2D - rock structure
%
    moduleCheck('libgeometry');
    [grdecl dataset petroinfo] = ...
        getAtlasGrid('Utsirafm', 'coarsening', coarsening_level);

    % Computing the Utsira top-surface grid
    G = processGRDECL(grdecl{1});  G = G(1);
    G = mcomputeGeometry(G);
    if do_cut
        % cut away eastern part of the grid
        east_bound = 6575000;
        rm_cells = G.cells.centroids(:,2) > east_bound;
        G = removeCells(G, rm_cells);
    end    
    Gt = topSurfaceGrid(mcomputeGeometry(G));
    
    % Setting up the rock structure
    petrodata   = petroinfo{1};   
    rock2D.perm = repmat(petrodata.avgperm, Gt.cells.num, 1);
    rock2D.poro = repmat(petrodata.avgporo, Gt.cells.num, 1); 
end
