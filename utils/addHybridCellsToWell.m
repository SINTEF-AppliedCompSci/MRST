function W = addHybridCellsToWell(W, G, rock, varargin)
%Add hybrid cells adjacent to a set of wells. For use with discrete fracture models (DFM)
%
% SYNOPSIS:
%   W = addHybridCellsToWell(W, G, rock)
%   setup = spe10_wo(W, G, rock, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   For each well in W, this function finds all hybrid cells adjacent ot
%   the cells of the well and adds them to the well. Hybrid are grid faces
%   that are assigned a volume in the computational grid. This function is
%   useful if the wells are added using a function that does not recognize
%   hybrid cells, e.g., addWellFromTrajectory.
%
% REQUIRED PARAMETERS:
%   W    - MRST well structure
%
%   G    - MRST grid structure, with hybrid cells
%
%   rock - MRST rock structre
%
%   fluid  - Fluid object as defined by function 'initSimpleFluid'.
%
% OPTIONAL PARAMETERS:
%   aperture - Hybrid cell aperture. Used for computing well indices.
%
%   any      - All other optional input parameters are passed on to
%             computeWellIndex

% SEE ALSO:
%   addHybrid

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

    % Optional input arguments
    opt          = struct('aperture', []);
    [opt, extra] = merge_options(opt, varargin{:});

    % Topology (+1 to avoid zero index)
    N = G.faces.neighbors + 1;
    % Hybrid cell mask (pad with leading false to comply above trick)
    hc  = G.cells.hybrid > 0;
    hc  = [false; hc];
    % Number of matrix cells
    ncm = G.cells.num - nnz(G.cells.hybrid);
    for i = 1:numel(W)
        % Well cell mask
        wc = false(G.cells.num+1, 1); wc(W(i).cells+1) = true;
        % Find well cells that are adjecent to hybrid cells
        ix = any(hc(N), 2) & any(wc(N), 2);
        if any(ix)
            cells  = reshape(N(ix,:), [], 1)-1;
            hcells = unique(cells(cells > ncm));
            W(i)   = addWellCells(W(i), G, rock, hcells, ...
                opt.aperture, extra{:});
        end
    end
    
end

%-------------------------------------------------------------------------%
function W = addWellCells(W, G, rock, cells, aperture, varargin)
% Add cells

    nc0     = numel(W.cells);
    W.cells = [W.cells; cells];
    nc      = numel(W.cells);
    % Safeguard against negative values by enlarging cellDims if necessary
    [dx, dy, dz] = cellDims(G, cells);
    if ~isempty(aperture)
        dz = aperture;
        if numel(dz) == G.cells.num, dz = dz(cells); end
    end
    re   = 2*0.14*sqrt(dx.^2 + dy.^2)/2;
    C    = max(1.1*W.r(1)./re,1);
    dx   =  [dx.*C, dy.*C, dz];
    % Compute well index (should probably be calculated more rigourously)
    WI = computeWellIndex(G, rock, W.r(1), cells, ...
        'cellDims', dx, varargin{:});
    W.WI = [W.WI; WI];
    % Add delta z field
    dZ   = getDeltaZ(G, cells, W.refDepth);
    W.dZ = [W.dZ; dZ];
    [~, ix] = sort(W.dZ);
    
    W.cells = W.cells(ix);
    W.WI = W.WI(ix);
    W.dZ = W.dZ(ix);
    
    % Add remaining fields by repeating W.(fn)(1)
    for fn = fieldnames(W)'
        v = W.(fn{1});
        if size(v,1) ~= nc0, continue; end
        W.(fn{1}) = repmat(v(1,:), nc, 1);
    end
    
end

%-------------------------------------------------------------------------%
function dZ = getDeltaZ(G, cells, refDepth)
% Compute distance from ref height to perforation (Copied from addWell)

    direction = gravity();
    dims      = G.griddim;
    if norm(direction(1:dims)) > 0
       direction = direction ./ norm(direction(1:dims));
    else
       direction = zeros(1, dims);
       if dims > 2
          direction(end) = 1;
       end
    end
    xyz = G.cells.centroids(cells, :);
    xyz(isnan(xyz)) = 0;
    Z = xyz * direction(1:dims).';
    dZ = Z - refDepth;
    
end