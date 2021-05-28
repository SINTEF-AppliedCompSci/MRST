function x = assignValue(x, v, inx)
% Assign values to ADI object by way of indices, without changing jacobians
%
% SYNOPSIS:
%   x = assignValue(x, v, inx)
% DESCRIPTION:
%   Replace the numerical values of a ADI or double, without changing the
%   Jacobians. This can lead to inconsistent Jacobians and variables so it
%   should only be used if you really know what you are doing!
%
% REQUIRED PARAMETERS:
%   x - ADI or double where values are to be replaced.
%
%   v - Values that will replace some subset of x.
%
% inx - The indices into x that v will replace. That is, after the call, 
%       double(x(ix)) == v
%
% RETURNS:
%   x - Modified version of input with same class.
%
% SEE ALSO:
%   ADI
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
    if isa(x, 'ADI')
        x.val(inx) = v;
    else
        x(inx) = v;
    end
end
