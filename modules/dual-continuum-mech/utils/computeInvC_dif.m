function invC_dif = computeInvC_dif(G, C_m, C_f)
%
% SYNOPSIS:
%   function invC_dif = computeInvC_dif(G, C_m, C_f)
%
% DESCRIPTION: 
%   Function used to calcualte the quantity (C_m - C_f)^-1 (see Ashworth and
%   Doster 2020 for further details). 
%
% PARAMETERS:
%   G   - Grid struc
%   C_m - matrix stiffness tensor
%   C_f - fracture stiffness tensor
%
% RETURNS:
%   invC_dif 
%
% EXAMPLE:
%
% SEE ALSO:
%
%{
Copyright 2009-2020 SINTEF ICT, Applied Mathematics.

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
if G.griddim == 3
    nlin = 6;
else
    nlin = 3;
end

% First C_dif
C_dif = C_m - C_f;
C_dif = repmat(reshape(permute(reshape(C_dif, [], nlin), [2,1]), [], nlin),...
             1,G.cells.num);
blockID = kron(eye(G.cells.num),ones(nlin));
C_dif = C_dif.*blockID;

% Perform inversion and reshape back to original form
invC_dif = (C_dif\eye(G.cells.num*nlin))';
k = find(blockID');
invC_dif = reshape(invC_dif(k), nlin^2,[])';
end