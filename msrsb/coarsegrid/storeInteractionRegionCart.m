function CG = storeInteractionRegionCart(CG, varargin)
% Store MPFA-like interaction region for coarse grid made from logical partitioning.
%
% SYNOPSIS:
%   CG = storeInteractionRegionCart(CG)
%
% DESCRIPTION:
%   This function quickly computes the interaction regions based on logical
%   indices for a coarse grid. The interaction region for a given coarse
%   block I is defined as any cells within the convex hull of the coarse
%   centroids from the nodal neighbors of I. This version relies on logical
%   indices for speed, but it is not applicable for general grids.
%
% REQUIRED PARAMETERS:
%   CG - Coarse grid.
%
% RETURNS:
%   CG - Coarse grid with added field "CG.cells.interaction", which is a
%        cell array of length equal to CG.cells.num. Each entry j contains
%        a list of fine cells designated as interacting with coarse block
%        j.
%
% SEE ALSO:
%   `storeInteractionRegion`

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

    opt = struct('adjustCenters', true, ...
                 'centerOverride',    [], ...
                 'edgeBoundaryCenters', true);
    opt = merge_options(opt, varargin{:});

    % For when grid was created by partitionUI
    G = CG.parent;
    
    if ~isfield(CG.faces, 'centers') || ~isfield(CG.cells, 'centers')
        CG = addCoarseCenterPoints(CG, ...
                                   'adjustCenters', opt.adjustCenters, ...
                                   'centerOverride', opt.centerOverride, ...
                                   'edgeBoundaryCenters', opt.edgeBoundaryCenters);
    end
    ijk = gridLogicalIndices(G);
    ijk = [ijk{:}];
    
    interaction = cell(CG.cells.num, 1);
    for i = 1:CG.cells.num
        dispif(mrstVerbose, 'Interaction stored for block %d of %d\n', i, CG.cells.num);
        n = coarseNeighbors(CG, i, false);
        if isempty(n)
            n = i;
        end
        ijk_near = ijk(CG.cells.centers(n), :);
        
        ijkloc = ijk(CG.partition == i, :);
        
        M_ix = max(ijk_near, [], 1);
        m_ix = min(ijk_near, [], 1);
        
        M_ix = max(M_ix, max(ijkloc + 1, [], 1));
        m_ix = min(m_ix, min(ijkloc - 1, [], 1));
        
        ok = true(G.cells.num, 1);
        for j = 1:G.griddim
            ok = ok & ijk(:, j) < M_ix(j);
            ok = ok & ijk(:, j) >  m_ix(j);
        end
        interaction{i} = find(ok);
    end
    CG.cells.interaction = interaction;
end
