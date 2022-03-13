function G = processAAIGrid(meta, data, topgrid, cstrids)
% Process aii grid meta data to a grid
% SYNOPSIS:
%   G = processAAIGrid(meta, data, topgrid, cstrids)
%
% PARAMETERS:
%   G       - Grid data structure.
%
%   meta    - meta data of the grid
%
%   data    - data defining height of surface
%
%   topgrid - if true make topsurface grid
%
%   cstrids - stride to make coarser representation
%
% RETURNS:
%   G - valid MRST grid. If 'topgrid' is true, G has the format of
%       a top-surface grid. If not, G is a 2D grid embedded in 3D

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    mlist = mrstModule();
    
    dims = [meta.ncols, meta.nrows];
    % We have a cell centered grid, so subtract by one
    gdims = dims - 1;

    % We have cell areas and equidistant grid, find h
    h = meta.cellsize;

    if(~isempty(cstrids))        
        data=data(1:cstrids(1):end,1:cstrids(2):end);
        dims=size(data);
        gdims=dims-1;
        G = cartGrid(gdims, gdims.*h.*cstrids);
    else
        G = cartGrid(gdims, gdims.*h);  
    end
    G.nodes.coords(:,3) = abs(data(:));

    G.nodes.coords(:,1) = G.nodes.coords(:,1) + meta.xllcorner;
    G.nodes.coords(:,2) = G.nodes.coords(:,2) + meta.yllcorner;
    
    % We exploit the fact that anything computed using a nan results in
    % another nan.
    s = warning('off', 'GridType:unsupported');
    G = computeGeometry(G);
    badcells = any(isnan(G.cells.centroids),2);
    
    z = G.cells.centroids(~badcells, 3);
    
    % Extract the subgrid of the actual cells, set gdims correctly
    active = find(~badcells);
    G = extractSubgrid(G, active);
    G.cartDims = gdims;
    G.cells.indexMap = active;
    
    % Make 2D plots much faster by explicitly storing cell->node map.
    G.cells.sortedCellNodes = getSortedCellNodes(G);
    
    if topgrid
        % Make the grid a top surface grid
        G.cells.z = z;
        G.nodes.z= G.nodes.coords(:,3);
        G.nodes.coords(:,3) = [];
        %{
        if(all(strid==strid(1))
            dx=h*strid(1);
            G.faces.areas     = dx  
            G.faces.normals   = [];
            G.faces.centroids = ;
            G.cells.volumes   = cellVolumes;
            G.cells.centroids = cellCentroids;
        %}
        G = computeGeometry(G);
    end
    G.cells.z = z;
    
    warning(s);
    mrstModule('reset',mlist{:});
end
