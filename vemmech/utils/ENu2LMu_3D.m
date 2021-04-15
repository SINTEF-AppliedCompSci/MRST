function [lambda, mu] = ENu2LMu_3D(E, nu)
%
%
% SYNOPSIS:
%   function [lambda, mu] = ENu2LMu_3D(E, nu)
%
% DESCRIPTION: Compute the Lamé parameters lambda and mu in the 2D case
% from the Young's modulus E and Poisson's ratio Nu. In the orthogonal
% direction (with respect to the 2D plane), it is assumed zero
% displacement. In particular, it implies that Nu = 0.5 corresponds as in
% the 3D case to an incompressible material. An other assumption, which is
% not used here, could be that there is no stress in the orthogonal
% direction.
%
% PARAMETERS:
%   E  - Young's modulus
%   nu - Poisson's ratio
%
% RETURNS:
%   lambda - First Lamé coefficient
%   mu     - Second Lamé coefficient
%
% EXAMPLE:
%
% SEE ALSO:
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

   lambda = (E .* nu) ./ ( (1 + nu) .* (1 - 2 * nu));
   mu = E ./ (2 * (1 + nu));

end
