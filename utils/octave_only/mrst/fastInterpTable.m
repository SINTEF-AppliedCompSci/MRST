function yi = fastInterpTable(X, Y, xi)
% Fast interpolation of table, using interp1
%
% SYNOPSIS:
%   yi = fastInterpTable(X, Y, xi)
%
% DESCRIPTION:
%   A simple wrapper for griddedInterpolant for fast interpolation of
%   simple data in the AD-framework. Always defaults to linear
%   interpolation with linear extrapolation.
%
% REQUIRED PARAMETERS:
%   x  - Sample X-coordinates
%
%   y  - Sample function values
%
%   xi - The X-coordinates at which the linear interpolant is to be
%        evaluated.
%
% RETURNS:
%   yi - Linear function interpolating (x, y) evaluated at xi.

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

   if isempty(xi)
       yi = [];
       return
   end
    yi = interp1(X, Y, xi, 'linear','extrap');
end

