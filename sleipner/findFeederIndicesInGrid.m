function [feederCells,feederFaces,feederCellsMap] = findFeederIndicesInGrid(G,feeders,layerMapG,varargin)
% Get indices of grid cells in finescale Sleipner grid corresponding to
% feeder chimney locations from the Sleipner 2019 Benchmark dataset.
% See Sleipner 2019 Reference Model - Polygons_feeders_chimney.pptx from
% the Sleipner 2019 dataset for more details.
%
% SYNOPSIS:
%   [feederCells] = findFeederIndicesInGrid(G,feeders,layerMapG);
%
%
% PARAMETERS:
%  G        - Finescale Grid of Multilayered Sleipner Benchmark model without
%              caprock cells. Output from getMultiLayerSleipnerGrid.
%
%  feeders   - Output structure from getSleipnerFeederOutlines2019.
%
%  layerMapG - array indicating which layer each cell in G belongs to.
%               Output from getSleipnerLayerBoundaryFaces.
%
% OPTIONAL PARAMETERS:
%  feederId - ID of which feeder chimneys to include. 
%                 1 - main feeder chimney
%                 2 - NE feeder
%                 3 - SW feeder
%
% RETURNS:
%  feederCells - array of cell indices of feeder cells in fine scale
%                   Sleipner grid.
%
%  feederFaces - array of face indices of feeder cells in fine scale
%                Sleipner grid
% 
% feederCellsMap - array of size G.cells.num indicating which feeder (if
%                  any) a cell belongs to.
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

opt = struct('feederId',  [1 2 3],...
             'includeCaprock',false);
opt = merge_options(opt, varargin{:});


fIx = opt.feederId;
cartDims = G.cartDims;

feederCells = [];
feederFaces = [];
feederCellsMap = ones(G.cells.num,1);


for k = 1:numel(fIx)
    
    [c,on] = inpolygon(G.cells.centroids(:,1),G.cells.centroids(:,2),...
        feeders{fIx(k)}.outline{1}(:,1),feeders{fIx(k)}.outline{1}(:,2));
    c = find(c-on);
    [gi,gj,~] = gridLogicalIndices(G,c);
    
    feederCellsTemp = [];    
    
    % Add extra layer if layer 9 is included (layerMapG labels the thick
    % shale layer as 9 and L9 (sand layer above it) as 10.
    if ismember(9, feeders{fIx(k)}.layernos)
        layers = [feeders{fIx(k)}.layernos 10];
    else
        layers = feeders{fIx(k)}.layernos;
    end

    if fIx(k)~=1
        try
        % Pick a single cell from feeder polygon for the two less well
        % constrained chimneys
        
        % Cells are chosen to be same as those in "Sleipner 2019 Reference
        % Model - Polygons_feeders_chimney.pptx" Slide 3.
        fcIx = [2 3; 3 3];
        
        cellIxi = unique(gi);
        cellIxj = unique(gj);
        iIx =  cellIxi(fcIx(fIx(k)-1,1));
        jIx =  cellIxj(fcIx(fIx(k)-1,2));
        
        ijIx = logical((gi==iIx).*(gj==jIx));
        gi = gi(ijIx);
        gj = gj(ijIx);
        catch
            warning('Cannot add feeders 2 and 3. Probably using a coarsened mesh')
            break
        end
    end
    
    for j = 1:numel(layers)
        
        layerNo = layers(j);
        if opt.includeCaprock
            gridInds = zeros(cartDims(1),cartDims(2),cartDims(3));
        else
            gridInds = zeros(cartDims(1),cartDims(2),cartDims(3)-10);
        end
        
        gridInds(gi,gj,:) = 1;
        gridInds = reshape(gridInds,[],1);
        gridInds = gridInds.*(layerMapG == layerNo);
        feederCellsTemp = [feederCellsTemp; find(gridInds)];


 
        
    end

    feederCells = [feederCells; feederCellsTemp]; 
    feederCellsMap(feederCellsTemp) = k+1;    
    feederFaces = [feederFaces; gridCellFaces(G, feederCellsTemp)];
end
