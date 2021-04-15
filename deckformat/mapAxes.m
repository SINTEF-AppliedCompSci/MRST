function mpos = mapAxes(pos, MA)
% Transform lateral position from map to local coordinate system
%
% SYNOPSIS:
%   p = mapAxis(pos, ma)
%
% PARAMETERS:
%   pos - nx2 vector of map coordinates that are to be transformed
%
%   ma  - 6x1 vector of points: x1 y1 x2 y2 x3 y3
%         Here: (x1,y1) describes a point on the y-axis, (x2,y2) is the
%         grid origin, and (x3,y3) is a point on the x-axis
%
% RETURNS:
%   p   - nx2 vector of coordinates transformed to local coordinate system
%
% EXAMPLE:
%   p = mapAxes([x y], [100 100 0 100 100 0])

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


e1 = MA(5:6) - MA(3:4); e1 = e1/norm(e1);
e2 = MA(1:2) - MA(3:4); e2 = e2/norm(e2);
msign = dot([0 0 1],cross([e1,0],[e2,0]));
if msign < 0
    warning('mapAxes:signChange','mapAxes may change signature of coordinate system.')
end
% tranform coordinate system may change the sign of the system
mpos = bsxfun(@times, pos(:,1), e1) + bsxfun(@times, pos(:,2), e2);
mpos = bsxfun(@plus, mpos, MA(3:4));
