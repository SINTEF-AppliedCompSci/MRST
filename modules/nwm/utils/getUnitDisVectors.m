function [ux, uy, uz] = getUnitDisVectors(G, cfCentersAll, cells)
% Get unit distance vectors of cells in Corner-point or Cartesian grid
%
% SYNOPSIS:
%   [ux, uy, uz] = getUnitDisVectors(G, cfCentersAll, cells)
%
% PARAMETERS:
%  G             - Corner-point or Cartesian grid structure
%  cfCentersAll  - Cell face center of the G, corresponding to G.cells.faces
%                  can be obtained by: 'computeCpGeometry' or 
%                  G.faces.centroids(G.cells.faces(:,1), :)
%  cells         - Cells of G
%
% RETURNS:
%  ux: Unit distance vector directing from X- face center to X+ face center
%  uy:          ......                     Y-                Y+ 
%  uz:          ......                     Z-                Z+ 
% 
% EXAMPLE:
%   G = cartGrid([10, 10, 5], [100, 100, 10]);
%   G = computeGeometry(G);
%   cfCentersAll = G.faces.centroids(G.cells.faces(:,1), :);
%   [ux, uy, uz] = getUnitDisVectors(G, cfCentersAll, (1:G.cells.num)');
%
% SEE ALSO:
%   `computeGeometry`, `computeCpGeometry`

    % Cell face directions
    cellFacesDir = arrayfunUniOut(@(x)G.cells.faces(G.cells.facePos(x): ...
        G.cells.facePos(x+1)-1, 2), cells);
    % Cell face centers
    cfCenters = arrayfunUniOut(@(x)cfCentersAll(G.cells.facePos(x): ...
        G.cells.facePos(x+1)-1, :), cells);
    cfCentersDir = cell(6,1);
    for dir = 1 : 6
        tmp = cellfun(@(x,y)y(x==dir, :), cellFacesDir, cfCenters, ...
            'UniformOutput', false);
        tmp = cellfunUniOut(@(x)sum(x,1)/size(x,1), tmp);
        tmp = cell2mat(tmp);
        cfCentersDir{dir} = tmp;
    end

    % Unit distance vectors
    ux = cfCentersDir{2} - cfCentersDir{1};
    ux = bsxfun(@rdivide, ux, sqrt( sum(ux.^2, 2) ));
    uy = cfCentersDir{4} - cfCentersDir{3};
    uy = bsxfun(@rdivide, uy, sqrt( sum(uy.^2, 2) ));
    uz = cfCentersDir{6} - cfCentersDir{5};
    uz = bsxfun(@rdivide, uz, sqrt( sum(uz.^2, 2) ));
end