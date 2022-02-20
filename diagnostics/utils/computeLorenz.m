function Lc = computeLorenz(F,Phi)
%Compute the Lorenz coefficient
%
% SYNOPSIS:
%   Lc = computeLorenz(F, Phi)
%
% PARAMETERS:
%   F   - flow capacity
%   Phi - storage capacity
%
% RETURNS:
%   Lc  - the Lorenz coefficient, a popular measure of heterogeneity. It is
%         equal to twice the area under the curve and above the F=Phi line.
%         It varies between 0 (homogeneous displacement) to 1 (infinitely
%         heterogeneous displacement).
%
% SEE ALSO:
%   `computeFandPhi`, `computeSweep`

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


v      = diff(Phi,1);
Lc     = 2*(sum((F(1:end-1,:)+F(2:end,:))/2.*v) - .5);
end

