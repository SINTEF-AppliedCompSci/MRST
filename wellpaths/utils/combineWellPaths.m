function wellpath = combineWellPaths(wellpaths)
%Combine multiple simple paths into a full tree
%
% SYNOPSIS:
%   wellpath = combineWellPaths({wp1, wp2, wp3});
%
% DESCRIPTION:
%   Given multiple simple well paths this function assembles a single well
%   tree from the inputs. For this to work we assume that:
%       - The paths are ordered by depth in the tree. This is not the
%       vertical depth, but rather that a path will always be connected to
%       one of the curves preceding it in the list.
%       - Paths (aside from the first one) always start with a point that
%       also exists in one of the preceding paths. This is used to connect
%       the paths.
%
% REQUIRED PARAMETERS:
%
%   wellpaths - Cell array of simple well paths to be assembled together. A
%               simple well path is assumed to contain a single list of
%               points (i.e. it will only represent a line segment).
%
% RETURNS:
%   wellpath  - Composite wellpath made from the simple wellpaths. The
%               topology will be tree-like in nature.
%
%
% SEE ALSO:
%   `makeSingleWellpath`

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

% Check number of segments in each provided well path.
pathlens = cellfun(@(x) numel(x.points), wellpaths);
assert(all(pathlens == 1), ...
    'One or more wellpaths contained multiple distinct curves');

% Take points and active flags into cell arrays.
points = cellfun(@(x) x.points{1}, wellpaths, 'UniformOutput', false);
active = cellfun(@(x) x.active{1}, wellpaths, 'UniformOutput', false);

% First point of each path - this are the points that we will have to
% glue/match with parts of other curves..
pts = getPathsStartPoints(points);

npth = numel(wellpaths);
% First row is index of parent path, second index is the cell in that
% local list we are connected to. So if conn(7, :) = [3, 5] this
% indicates that path 7 starts from point number 5 of path number 3.
conn = zeros(npth, 2);
% Start at number 2 - first entry is always assumed to be the beginning
% of the well.
for i = 2:npth
    % Match the starting point against all points we have seen so far
    % in the tree (i.e. already processed segments). We do not check
    % against points that are not processed, as we assume that wells
    % are trees (no cycles in the topology).
    start = pts(i, :);
    attached = cellfun(@(x) getMatch(x, start), points(1:i-1));

    % Take first match as parent. Matches beyond the first represent
    % other curves that are attached to the same point.
    parent = find(attached ~= 0, 1);
    if isempty(parent)
        msg = ['Unable to parse topology. First path of path ', num2str(i),...
               ' did not connect to any points on preceding path. Aborting.'];
        error(msg)
    end
    % Parent and then include the index of the point we matched in the
    % parent's indexing system.
    conn(i, :) = [parent, attached(parent)];
end
% Finally make the trajectory.
wellpath = makeSingleWellpath(points, conn, active);
end

function pts = getPathsStartPoints(points)
    tmp = cellfun(@(x) x(1, :), points, 'uniformoutput', false)';
    
    pts = vertcat(tmp{:});
end

function sub = getMatch(pts, pt)
    sub = find(all(bsxfun(@eq, pts, pt), 2));
    if isempty(sub)
        sub = 0;
    end
end