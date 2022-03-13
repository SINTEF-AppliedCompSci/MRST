function [CG, CS] = speedUpMS(G, CG, CS, type)
%Precompute MsMFE basis reduction matrices in order speed up assembly
%
% SYNOPSIS:
%   [CG, CS] = speedUpMS(G, CG, CS, solver)
%
% PARAMETERS:
%   G, CG  - Grid and coarse grid data structures respectively.
%
%   CS     - Coarse-scale linear system data structure as defined by
%
%   solver - Particular solver mode for which to precompute the basis
%            reduction matrices.  This mode defines the kind of system of
%            linear equations that will later be assembled and solved in
%            the 'solveIncompFlowMSSpeedUp' function.
%
%            String.  Must be one of 'hybrid' or 'string'.
%
% RETURNS:
%   CG, CS - Updated coarse-scale grid and linear system data structures
%            with additional fields designed to lower the cost of forming
%            the coarse-scale system of linear equations from which to
%            derive block fluxes and block pressures.
%
% NOTE:
%   If you wish to change the solver mode in function
%   'solveIncompFlowMSSpeedup', you must call 'speedUpMS' with a different
%   'solver' parameter.
%
% SEE ALSO:
%   `generateCoarseGrid`, `generateCoarseSystem`, `solveIncompFlowMSSpeedUp`.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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


% preassemble basis
if strcmpi(type, 'hybrid'),
   [CS.Bv, CS.Phi] =  basisMatrixHybrid(G, CG, CS);
else
   [CS.Bv, CS.Phi] = basisMatrixMixed (g, cg, cs);
end

% compute some indices (which fine faces belongs to which coarse faces)
[nsub, sub] = subFaces(G, CG);
sub_ix      = cumsum([0; nsub]);
sub         = sub(mcolon(sub_ix(CS.activeFaces) + 1, ...
sub_ix(CS.activeFaces + 1)));
nsub        = nsub(CS.activeFaces);
CG.nsub = nsub;
CG.sub = sub;


end
