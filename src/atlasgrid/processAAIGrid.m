function G = processAAIGrid(meta, data, topgrid, cstrids)
% Process aii grid meta data to a grid
% SYNOPSIS:
%   G = processAAIGrid(meta, data, topgrid, cstrids)
%
% PARAMETERS:
%   G       - Grid data structure.
%
%   meta  - meta data of the grid
%
%   data  - data defining hight of surface
%
%   topgrid - if true make topsurface grid
%
%   cstrids - stride to make coarser representation
%
% RETURNS:
%   G - valid mrst grid, if topgrid is true it has the format of
%       at topsurface grid or else it is a 2D grid embedded in 3D
% 
% NOTES:
%
% EXAMPLE:
%
% SEE ALSO:
%     
%
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
