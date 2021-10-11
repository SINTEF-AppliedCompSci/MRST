function nk = polyDim(k, dim)
%Computes the dimension of the space of polynomials of degree k or less in
%R^dim.
%
%   SYNOPSIS:
%       polyDim(k, dim)
%
%   REQUIRED PARAMETERS:
%       k   - Polynomial order.
%
%       dim - Dimension of domain on which the polynomails are defined.
%
%   RETURNS:
%       nk  - dimension of polynomial function space.
%

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

%   Written by Øystein Strengehagen Klemetsdal, SINTEF/NTNU, 2016.

if k < 0
    nk = 0;
else
    nk = nchoosek(k+dim,k);
end
    
end