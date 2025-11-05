function [model, poro0, perm0] = setupBioCloggingModel(model, nbact0, nc, cp)
% setupBioCloggingModel -- Add bio-clogging effects to a compositional model
%
% SYNOPSIS:
%   [model, poro0, perm0] = setupBioCloggingModel(model, nbact0, nc, cp)
%
% DESCRIPTION:
%   Modifies the input MRST compositional model by including porosity
%   and permeability reduction due to bacterial growth (bio-clogging).
%   The reduction is parameterized by a characteristic bacterial
%   concentration and a clogging strength coefficient.
%
% PARAMETERS:
%   model  - MRST compositional model struct with fields `rock` and `fluid`.
%   nbact0 - Initial bacterial concentration (scalar, dimensionless/normalized).
%   nc     - Characteristic bacterial concentration at which clogging
%            effects become significant (scalar).
%   cp     - Dimensionless clogging strength coefficient (scalar).
%
% RETURNS:
%   model  - Modified model with porosity and permeability given as
%            function handles depending on bacterial concentration.
%   poro0  - Original porosity vector (before clogging modification).
%   perm0  - Original permeability vector (before clogging modification).
%
% SEE ALSO:
%   initCompositionalStateBacteria, TableCompositionalMixture

% Store original rock properties
poro0 = model.rock.poro;
perm0 = model.rock.perm(:, 1);

% Define porosity multiplier (function of bacterial concentration)
scale = 1 + cp.*(nbact0 / nc).^2;
pvMult_nbact = @(nbact) 1 ./ (1 + scale.*(nbact./nc).^2);

% Assign porosity update function to fluid/rock
model.fluid.pvMultR = @(p, nbact) pvMult_nbact(nbact);
poroFun = @(p, nbact) poro0 .* pvMult_nbact(nbact);
model.rock.poro = poroFun;

% Define permeability update function via Kozeny–Carman-like relation
tauFun = @(p, nbact) ((1 - poro0) ./ (1 - poroFun(p, nbact))).^2 .* ...
    (poroFun(p, nbact) ./ poro0).^3;
permFun = @(p, nbact) perm0 .* tauFun(p, nbact);
model.rock.perm = permFun;
end

%{
Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify it under the terms of the GNU 
General Public License as published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

MRST is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with MRST. If not, 
see <http://www.gnu.org/licenses/>.
%}