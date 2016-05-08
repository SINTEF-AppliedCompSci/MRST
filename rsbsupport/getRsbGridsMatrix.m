function CGm = getRsbGridsMatrix(G, pm, varargin)
% Generates the matrix coarse grid and computes support regions for the
% matrix coarse blocks. See generateCoarseGrid and storeInteractionRegion.

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


if isstruct(varargin{1})
    opts = varargin{1};
    coarseCenterOpts = {'edgeBoundaryCenters', true};
else
    opts = struct(  'fullyCoupled'        , true , ...
                    'Wells'               , []   );
    [opts, coarseCenterOpts] = merge_options(opts, varargin{:});
end

Gm = G.Matrix; nt = G.cells.num;
CGm = generateCoarseGrid(Gm, pm);
CGm = coarsenGeometry(CGm);
if isfield(Gm,'cartDims') && ~strcmp(G.type{1,1},'processGRDECL')
    CGm = addCoarseCenterPoints(CGm, coarseCenterOpts{:});
    if ~isempty(opts.Wells)
        CGm = changeCoarseCenterToWellLocation(CGm, opts.Wells);
    end
    CGm = storeInteractionRegionCart(CGm);
else
    CGm = addCoarseCenterPoints(CGm, coarseCenterOpts{:});
    if ~isempty(opts.Wells)
        CGm = changeCoarseCenterToWellLocation(CGm, opts.Wells);
    end
    CGm = storeInteractionRegion(CGm);
end
if opts.fullyCoupled
    % Add fracture cells within region
    fcent = G.cells.centroids(Gm.cells.num+1:nt,:);
    for i = 1:CGm.cells.num
        tempIR = CGm.cells.interaction{i,1};
        nn = unique(gridCellNodes(Gm,tempIR));
        tri = delaunayTriangulation(G.nodes.coords(nn,:));
        in = pointLocation(tri,fcent);
        add = Gm.cells.num + find(~isnan(in));
        CGm.cells.interaction{i,1} = [tempIR;add];
    end
end
return