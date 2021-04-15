function wellpath = makeSingleWellpath(pts, conn, active)
%Create well path from points (and optional connectivity and active flag)
%
% SYNOPSIS:
%   wellpath = makeSingleWellpath(pts);
%   wellpath = makeSingleWellpath(pts, conn)
%   wellpath = makeSingleWellpath(pts, conn, active)
%
% DESCRIPTION:
%   Create a well path from lists of points. The resulting structure will
%   be a single struct with fields:
%
%   - points:
%             Cell array, each consisting of N x Dim arrays of n points.
%             Each entry contains the points for one segment that are
%             assumed to be connected as a line according to their
%             ordering. The first entry is assumed to be closest to the
%             starting point of the well (closest here means along the well
%             bore).
%
%   - connectivity:
%            Array of size M x 2 where M is the number of entries in the
%            .points cell array. If entry number 5 of connectivity is [2 7]
%            it means that segment number 5 is connected to segment 2, at
%            point number 7 of segment 2's internal ordering. In effect,
%            segment 5 branches off from segment 2 from the coordinate
%            .points{2}(7, :).
%
%   - active: 
%            Cell array, containing active flags for the segments between
%            points. If points{i} contains N x dim entries, active{i}
%            should contain (N-1) x 1 entries, indicating if the
%            subsegments are active. 
%
%            If points is of size 6 x 3 and active looks like this:
%            [1; % 1 -> 2
%             1; % 2 -> 3
%             0; % 3 -> 4
%             1; % 4 -> 5
%             1] % 5 -> 6
%            it means that the part of the segment will be disabled from
%            point 3 to point 4.
%
% REQUIRED PARAMETERS:
%   pts     - Maps directly to points (see above). If pts is a numeric
%             array, it will be interpreted as a cell array with a single
%             entry.
%
%  connectivity - Maps directly into connectivity.
%
%  active       - Maps directly into active.
%
% RETURNS:
%   wellpath - Wellpath suitable for plotting or producing well
%              completions.
%
% SEE ALSO:
%   `plotWellPath`, `getWellFromPath`

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

if ~iscell(pts)
    pts = {pts};
end
npts = numel(pts);
if nargin < 2
    % Topology is trivial, all points are connected on a line
    assert(npts == 1);
    conn = [0 0];
end

if nargin < 3
    % Default: All segments are active
    active = cellfun(@(x) true(size(x, 1) - 1, 1), pts, 'UniformOutput', false);
end

wellpath = struct('points', [], ...
                  'connectivity', [], ...
                  'active', []);
wellpath.points = pts;
wellpath.connectivity = conn;
wellpath.active = active;
end