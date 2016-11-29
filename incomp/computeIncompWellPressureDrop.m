function dp = computeIncompWellPressureDrop(W, mob, rho, g)
% Compute incompressible connection pressure drop for a single well
%
% SYNOPSIS:
%       dp = computeIncompWellPressureDrop(W, mob, rho, g)
%
% PARAMETERS:
%       W   - Well struct (see addWell)
%
%       mob - G.cells.num x n_ph array of cell mobilities
%
%       rho - n_ph array of densities
%
%       g   - Norm of gravity along z-axis
%
% RETURNS:
%       dp  - Column vector with numel(W.cells). Contains pressure drop
%       from the bottom hole pressures, under certain assumptions for
%       the incompressible, linear pressure equation.
%
% NOTES:
%   In order to avoid nonlinear behavior for wells, this function assumes
%   that the well bore density is constant for all perforations and is
%   equal to the injection composition if the well is an injector,
%   either through positive rates or W.sign > 0. If the well is a producer,
%   the mixture is assumed to be taken from the perforated cells, with flux
%   proportional to the well indices.
%
% SEE ALSO:
%   incompTPFA, incompMimetic, incompMPFA

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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

% A well is a injector if the sign is set or the control is set to a
% positive rate. Otherwise, we assume that it is a producer.
isInj = W.sign > 0 || (strcmpi(W.type, 'rate') && W.val > 0);
if isInj
    % We compute a simple mixture of the injected fluids
    rhoMix = rho*W.compi';
else
    % We compute approximate phase fluxes by taking the fractional
    % mobility in each cell, weighted by the densities. The flux is
    % estimated roughly as proportional to the well index.
    wc = W.cells;
    mobw = mob(wc, :);
    f = bsxfun(@rdivide, mobw, sum(mobw, 2));
    rhoMix = sum(W.WI.*(f*rho'))./sum(W.WI);
end
dp = g * W.dZ*rhoMix;
