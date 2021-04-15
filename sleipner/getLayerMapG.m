function layerMapG = getLayerMapG(G,varargin)
% Function to assign layer numbers to each cell in the full Sleipner 
% benchmark grid.
%
% SYNOPSIS:
%   layerMapG = getLayerMapG(G,varargin)
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
%   layerMapG  - Array indicating which layer each cell in the grid belongs
%                 to.Layer 1 is the deepest layer. Layer 9 indicates the 
%                 thick shale layer. 
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
    
    
    layerSize = G.cartDims(1)*G.cartDims(2);

    layerMapG = ones(G.cells.num,1); 
    for k = 1:numel(layerNum)
 
        if k == 1
            cellInds = 1:(layerSize*layerNum(k));    
            
        else
            cellInds = 1:(layerSize*(layerNum(k)-layerNum(k-1)));
            cellInds = cellInds + layerSize*(layerNum(k-1));
        end
        layerMapG(cellInds) = numel(layerNum)+2-k;
        
    end
    
end

