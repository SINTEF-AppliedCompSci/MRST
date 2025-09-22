function fluid = modifyRelPermForResidualSaturations(fluid, swc, swmax, sgr, sor, sgmax)
% modifyRelPermForResidualSaturations -- Shift relative permeability functions
%                                        to include residual saturations
%
% SYNOPSIS:
%   fluid = modifyRelPermForResidualSaturations(fluid, swc, swmax, ...
%                                               sgr, sor, sgmax)
%
% DESCRIPTION:
%   Adjusts the relative permeability functions of the input fluid model
%   to account for specified residual saturations of water, oil, and gas.
%   The original relperm curves are shifted such that:
%       - Water flow starts only after connate water saturation (swc).
%       - Gas flow starts only after residual gas saturation (sgr).
%       - Oil relative permeability is reduced according to sor and sgr.
%
%   The method interpolates the original relperm tables and constructs
%   new function handles with the shifted saturation ranges.
%
% PARAMETERS:
%   fluid  - Fluid struct with fields `krW`, `krG`, `krOW`, `krOG`.
%   swc    - Connate (irreducible) water saturation.
%   swmax  - Maximum water saturation (upper endpoint).
%   sgr    - Residual gas saturation.
%   sor    - Residual oil saturation.
%   sgmax  - Maximum gas saturation (upper endpoint).
%
% RETURNS:
%   fluid  - Modified fluid struct with updated relative permeability
%            function handles and endpoints stored in `fluid.krPts`.
%
% SEE ALSO:
%   setupBioCloggingModel, TableCompositionalMixture

% Define synthetic saturation ranges for interpolation
SW = linspace(0, 1, 100);
SO = linspace(0, 1, 100);
SG = linspace(0, 1, 100);

% Evaluate original relative permeability functions
krW_original  = arrayfun(fluid.krW,  SW);
krOW_original = arrayfun(fluid.krOW, SO);
krG_original  = arrayfun(fluid.krG,  SG);
krOG_original = arrayfun(fluid.krOG, SO);

% --- Water ---
SW_shifted = max(SW - swc, 0);
krW_new    = interp1(SW, krW_original, SW_shifted, 'linear', 0);
fluid.krW  = @(sw) interp1(SW, krW_new, value(max(sw - swc, 0)));

% --- Gas ---
SG_shifted = max(SG - sgr, 0);
krG_new    = interp1(SG, krG_original, SG_shifted, 'linear', 0);
fluid.krG  = @(sg) interp1(SG, krG_new, value(max(sg - sgr, 0)));

% --- Oil ---
SO_shifted  = max(SO - sor, 0);
krOW_new    = interp1(SO, krOW_original, SO_shifted, 'linear', 0);
krOG_new    = interp1(SO, krOG_original, SO_shifted, 'linear', 0);
fluid.krOW  = @(so) interp1(SO, krOW_new, value(max(so - sor, 0)));
fluid.krOG  = @(so) interp1(SO, krOG_new, value(max(so - sor, 0)));

% Store updated relative permeability endpoints
fluid.krPts.w  = [swc,  swmax];
fluid.krPts.g  = [sgr,  sgmax];
fluid.krPts.ow = [0,    1 - swc];
fluid.krPts.og = [0,    1 - sgr];
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