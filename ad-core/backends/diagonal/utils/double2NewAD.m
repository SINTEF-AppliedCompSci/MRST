function u = double2NewAD(u, sample)
% Convert a double to NewAD variable, using a sample NewAD variable for dimensions
%
% SYNOPSIS:
%   u = double2NewAD(u, adivar)
%
% REQUIRED PARAMETERS:
%   u      - Double to be converted to NewAD.
%
%   sample - Sample variable of the type NewAD to be used.
% RETURNS:
%
%   u      - Variable with same type as sample and same value as u
%            initially had. If u is a NewAD class instance, u will have zero
%            jacobians with the same number of primary variables as the
%            jacobians of sample.
%
% SEE ALSO:
%   `NewAD`

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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


assert(isa(u,'double'));
if isnumeric(sample)
    % The dummy variable is also a double, we simply return without
    % changing it.
    return
end
nval  = numel(u);
jac  = cellfun(@(j) makeZero(j, nval), ...
               sample.jac, 'UniformOutput', false);

u = NewAD(u,jac);
u.numVars = sample.numVars;
u.offsets = sample.offsets;
end

function x = makeZero(x, nval)
    if issparse(x)
        x = sparse([], [], [], nval, size(x, 2));
    else
        x = x.toZero(nval);
    end
end