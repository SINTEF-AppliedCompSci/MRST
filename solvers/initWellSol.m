function wellSol = initWellSol(W, p0)
%Initialize well solution data structure.
%
% SYNOPSIS:
%   wellSol = initWellSol(W, p0)
%
% DESCRIPTION:
%   Initialize the well solution structure to uniform well bottom-hole
%   pressures and all-zero well rates.
%
% PARAMETERS:
%   W  - Well data structure as defined by addWell &c.
%   p0 - Initial uniform well bottom-hole pressure (scalar).
%
% RETURNS:
%   wellSol - Initialized reservoir solution structure having fields
%               - flux     -- Well rates in all perforations (== 0).
%               - pressure -- Well bottom-hole pressure (== p0).
%
% SEE ALSO:
%   `initResSol`, `solveIncompFlow`.

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


nW = numel(W);
wellSol = repmat(struct('flux', [], 'pressure', []), [1, nW]);

for w = 1 : nW,
   wellSol(w).flux     = zeros([numel(W(w).cells), 1]);
   wellSol(w).pressure = p0;%repmat(p0, [numel(W(w).cells), 1]);
end
end
