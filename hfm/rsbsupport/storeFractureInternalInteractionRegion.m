function CGf = storeFractureInternalInteractionRegion(CGf, G, A, varargin)
% storeFractureInternalInteractionRegion can be used to compute the support
% region that lies inside a fracture for every fracture coarse block.
%
% SYNOPSIS:
%   CGf = storeFractureInternalInteractionRegion(CGf, G, A)
%   CGf = storeFractureInternalInteractionRegion(CGf, G, A, 'pn1', 'pv1', ...)
%
% REQUIRED PARAMETERS:
%
%   CGf  - Fracture coarse grid (supplied by 'generateCoarseGrid') with
%          geometry information (computed through 'coarsenGeometry').
%
% OPTIONAL PARAMETERS:
%    simpleInteractionRegion - connects neighboring coarse nodes, in a 2D
%                              domain where fractures are represented as 1D
%                              grids to compute the interaction region.
%
%
%    callStoreInteractionRegion - calls the function storeInteractionRegion
%                                 to compute fracture internal interaction
%                                 regions. Applicable for 3D domains.
%
% RETURNS:
%   CGf - Fracture coarse grid containing internal support regions in the
%         list CGf.cells.interaction.
%
%
% SEE ALSO:
%   storeInteractionRegionFrac

%{
Copyright 2009-2015: TU Delft and SINTEF ICT, Applied Mathematics.

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

opt = struct(...
    'simpleInteractionRegion'    , true  , ...
    'callStoreInteractionRegion' , false );
    
opt = merge_options(opt, varargin{:});

if opt.callStoreInteractionRegion
    CGf = storeInteractionRegion(CGf,...
    'adjustCenters'              , false  , ...
    'skipSingleCellBlocks'       , false  , ... 
    'simpleCellGrouping'         , true   , ...
    'ensureConnected'            , true   , ...
    'largeBasis'                 , true   , ...  % true is better
    'localTriangulation'         , false  , ...  % local triangulation will not work for fractures
    'useMultipoint'              , false  );  
elseif opt.simpleInteractionRegion
    CGf = addCoarseCenterPointsFrac(CGf,'option','useCoarseFaceCentroids','meantype','arithmetic'); % IR centers
    Af = A(G.Matrix.cells.num+1:end,G.Matrix.cells.num+1:end);
    Af = Af - diag(diag(Af));
    p = CGf.partition;
    for i = 1:max(p)
        rmcenters = [CGf.cells.centers(1:i-1);CGf.cells.centers(i+1:end)];
        fcg_fcells = find(p==i);
        irf = fcg_fcells;
        start = irf;
        count = 0;
        while ~isempty(start)
            start2 = [];
            for j = 1:numel(start)
                add = find(Af(start(j),:));
                if isempty(rmcenters)
                    add = setdiff(add,irf);
                else
                    add = setdiff(add,[irf;rmcenters]);
                end
                start2 = [start2;add']; %#ok
            end
            start = start2;
            irf = [irf;start]; %#ok
            count = count + 1;
            if count==20
                break;
            end
        end
        CGf.cells.interaction{i,1} = unique(irf);
    end
end
return