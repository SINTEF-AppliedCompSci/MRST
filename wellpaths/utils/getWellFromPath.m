function [W, segInd, param] = getWellFromPath(W0, G, rock, wellpath, varargin)
% Convert well path to MRST well.
%
% SYNOPSIS:
%   W = getWellFromPath(W, G, rock, wellpath);
%
% DESCRIPTION:
%   This routine converts a well path (representing curves and points) into
%   a well (represented by cells and connectivity). 
%
% REQUIRED PARAMETERS:
%   W0       - Well array to be extended with the new well.
%
%   G        - The grid the well is to be placed in.
%
%   rock     - Rock structure which defines permeability and porosity.
%
%   wellpath - Well path as procued by makeSingleWellPath or
%              combineWellPaths.
%
%  OPTIONAL PARAMETERS:
%   
%   This function calls addWell. Any keyword arguments are
%   passed onto addWell.
%
% RETURNS:
%   W  - Updated wells
%
% NOTE:
%   Currently no special effort is made to ensure correct well indices for
%   the well. 
%
%
% SEE ALSO:
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

[cells, segInd, param, ptInd] = findWellPathCells(G, wellpath);
W = addWell([], G, rock, cells, varargin{:});


nperf = numel(cells);

topo = nan(nperf, 2);

% Treat connections in each segment first, as these are straightforward
offsets = cumsum([0; accumarray(segInd, 1)]);
nseg = max(segInd);
for i = 1:nseg
    nc = sum(segInd == i);

    of = offsets(i);

    topo((1:nc) + of, :) = getSegmentNeighborship(nc, of);
end

% Then deal with non-neighboring connections gluing the different
% segments together

nnc_topo = nan(nseg-1, 2);
for i = 2:nseg
    % Segment connection
    seg = wellpath.connectivity(i, 1);
    % Position in that array
    pos = wellpath.connectivity(i, 2);

    mask = false(nperf, 1);
    mask(segInd == seg) = true;

    connStart = find(ptInd == pos & mask, 1, 'last');
    connStop = find(segInd == i, 1, 'first');
    nnc_topo(i-1, :) = [connStart, connStop];
end

topo = [topo; nnc_topo];

% Remove zero entries, except for the first column which indicates
% where the well starts
mask = all(topo ~= 0, 2);
mask(1) = true;

W.topo = topo(mask, :);
W.status = true;
W.cstatus = true(nperf, 1);

W = [W0; W];
end

function topo = getSegmentNeighborship(nperf, offset)
    topo = [(0:(nperf-1))', (1:nperf)'];
    topo = topo + offset;
    topo(1, 1) = 0;
end

