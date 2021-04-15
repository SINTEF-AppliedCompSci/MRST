function [layerBoundaryFaces,layerFaceMapG] = getSleipnerLayerBoundaryFaces(G,varargin)
% Get faces in the 2019 Sleipner benchmark model which represent the top of
% each layer.
%
% SYNOPSIS:
%   [layerBoundaryFaces,layerFaceMapG] = getSleipnerLayerBoundaryFaces(G,varargin)
%
% REQUIRED PARAMETERS:
%   G - Full grid for 2019 Sleipner benchmark model created using
%       getMultilayerSleipnerGrid. Default assumption is that caprock
%       cells have been removed from the grid.
%
% OPTIONAL PARAMETERS:
%   includeCaprock - flag to indicate that grid G includes caprock
%                  cells. Default is false.
%
% RETURNS:
%   layerBoundaryFaces  - Faces in G corresponding to the top of each
%                           layer. Used when upscaling grid to layered VE.
%   layerFaceMapG       - Array indicating which layer each boundary face 
%                           belongs to. 0 - not a boundary face. 1 - 10
%                           face is at the top of that layer. Layer 1 is 
%                           the deepest layer. Layer 9 indicates the thick shale layer. 
%
%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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


   opt = struct('includeCaprock',false);
    opt = merge_options(opt, varargin{:});


fileID = fopen('SleipnerRefModel_GridLayers.txt');
C = textscan(fileID,'%s\t%f');
fclose(fileID);

    if opt.includeCaprock
        layers = C{2}(1:end);
        cumLayers = cumsum(layers);
        layerNum = cumLayers([1 2 3 5:2:17]);
    else
        layers = C{2}(2:end);
        cumLayers = cumsum(layers);
        layerNum = cumLayers([1 2 4:2:16]);
    end



layerSize = G.cartDims(1)* G.cartDims(2);
% 
%% Base Layer faces

boundaryFaces = zeros(G.faces.num,1);

layerFaceMapG = zeros(G.faces.num,1);

for k = 1:numel(layerNum)+1
    if k == 1
        cellInds = (1:1:layerSize);
    else
        cellInds = (1:1:layerSize) + layerNum(k-1).*layerSize;
    end
    [Gs,~,gf] = extractSubgrid(G,cellInds);
    fInds = gf(boundaryFaceIndices(Gs, 'Top', [],[]));
    boundaryFaces(fInds) = 1; 
    layerFaceMapG(fInds) = 11-k;

end

layerBoundaryFaces = find(boundaryFaces);


    
    
    
    
    
