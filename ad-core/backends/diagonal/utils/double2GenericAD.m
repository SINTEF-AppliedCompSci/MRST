function u = double2GenericAD(v, sample)
% Convert a double to GenericAD variable, using a sample GenericAD variable for dimensions
%
% SYNOPSIS:
%   u = double2GenericAD(u, adivar)
%
% REQUIRED PARAMETERS:
%   u      - Double to be converted to GenericAD.
%
%   sample - Sample variable of the type GenericAD to be used.
% RETURNS:
%
%   u      - Variable with same type as sample and same value as u
%            initially had. If u is a GenericAD class instance, u will have zero
%            jacobians with the same number of primary variables as the
%            jacobians of sample.
%
% SEE ALSO:
%   `GenericAD`

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

    assert(isa(v,'double'));
    if isnumeric(sample)
        % The dummy variable is also a double, we simply return without
        % changing it.
        u = v;
        return
    end
    u = sample.convertDouble(v);
end
