function [A, b] = reorderForILU(A, b, nc)
%Attempt to reorder a set of equations so that the diagonal is non-zero
%
% SYNOPSIS:
%   [A, b] = reorderForILU(A, b, nc);
%
% DESCRIPTION:
%   Reorder a set of equations to ensure non-zero diagonal. This is useful
%   when building ILU-based solvers. Notably, this utility is useful
%   whenever well equations are added, that may not have derivatives with
%   respect to all well controls.
%
% REQUIRED PARAMETERS:
%   A  - Linear system to be reordered.
%   b  - Right hand side of the system
%   nc - (OPTIONAL) The routine will only look at equations from nc+1 and
%        onwards to numel(b). 
% RETURNS:
%   A  - Reordered linear system
%   b  - Reordered right hand side

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
    if nargin == 2
        % Include all equations for consideration
        nc = 0;
    end
    d = diag(A);
    d_sub = d((nc+1):end);
    % Find bad entries
    bad = find(d_sub == 0);
    if isempty(bad)
        return
    end
    
    N = size(A, 1);
    
    A_sub = A((nc+1):end, (nc+1):end);
    
    [ii, jj, vv] = find(A_sub);
    
    nbad = numel(bad);
    newInd = zeros(nbad, 1);
    for i = 1:nbad
        % Switch equations around to avoid zero diagonal
        possible = jj(ii == bad(i));
        candidates = ii(jj == bad(i));
        
        new = intersect(possible, candidates);
        if isempty(new)
            warning('Unable to reorder equations, zeros on the diagonal still present');
            continue
        end
        newInd(i) = new(1);
        jj(ii == new(1)) = 0;
    end
    
    renum = (1:N);
    renum(nc + newInd) = nc + bad;
    renum(nc + bad) = nc + newInd;
    
    A = A(renum, :);
    b = b(renum, :);
end

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
