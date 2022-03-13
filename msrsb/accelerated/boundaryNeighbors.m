function [c_self, c_other] = boundaryNeighbors(g, subs)
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
    fa = boundaryFaces(g, subs);
    
    % Extract neighborship
    N = g.faces.neighbors(fa, :);
    % Remove faces on boundary
    bad = any(N == 0, 2);
    fa = fa(~bad);
    N = N(~bad, :);
    
    ok = false(g.cells.num, 1);
    ok(subs) = true;
    ok = [false; ok];
    
    isInt = ok(N(:, 1) + 1);

    [c_self, c_other] = deal(zeros(numel(fa), 1));
    
    c_self( isInt) = N( isInt, 1);
    c_self(~isInt) = N(~isInt, 2);
    
    c_other( isInt) = N( isInt, 2);
    c_other(~isInt) = N(~isInt, 1);
    c_self = unique(c_self);
    c_other = unique(c_other);

end