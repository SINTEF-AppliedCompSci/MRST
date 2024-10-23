function [layerMapG,layerMapShaleSand] = getLayerMapG(G,varargin)
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
%   layerMapShaleSand  - Array indicating which layer each cell in the grid belongs
%               to with separate shale and sand layers.
%               Layer 1 is the deepest sand layer. Shale layers are
%               listed from 101 - 110 (deepest = 101, caprock = 110).
%               If caprock is not included 109 is the shallowest shale layer
%               corresponding, to the thick shale layer.
%
%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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
    
    fileID = fopen(fullfile(fileparts(mfilename("fullpath")),'SleipnerRefModel_GridLayers.txt'));
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

    % Calculate layer map where each shale layer is given

    if nargout > 0
        [tempMap,layerMapShaleSand] = deal(zeros(G.cells.num,1));
        
        startIx = 0;

        for i = 1:numel(layers)

                ix =[1:(layerSize.*layers(i))] + startIx;                
                tempMap(ix) = i;                       
            startIx = ix(end);
        end


        if opt.includeCaprock
            shale_layers = [1:2:18];
            sand_layers = [2:2:18];
        else
            shale_layers = [2:2:17];
            sand_layers = [1:2:17];
        end
        
        for i = 1:numel(shale_layers)
            
            layerMapShaleSand(tempMap == shale_layers(i)) = 100 + numel(shale_layers) + 1 - i;
            
        end
        for i = 1:numel(sand_layers)
            
            layerMapShaleSand(tempMap == sand_layers(i)) = numel(sand_layers) + 1 - i;
            
        end
    end
    
end

