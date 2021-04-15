function fluid = getSPE10_model_1_fluid()
%Construct ADI Fluid Object for Model 1 of Tenth SPE CSP
%
% SYNOPSIS:
%   fluid = getSPE10_model_1_fluid
%
% PARAMETERS:
%   None.
%
% RETURNS:
%   fluid - Simplified ADI fluid object, suitable for use with simulator
%           function simulateScheduleAD, that represents the petrochemical
%           properties of the fluids in Model 1 of the Tenth SPE CSP.  The
%           fluids are immiscible and incompressible with constant
%           viscosities, but there is a non-trivial relative permeability
%           relation between the oil and gas phases.
%
% SEE ALSO:
%   `getSPE10_model_1_relperm`, `initSimpleADIFluid`, `initDeckADIFluid`.

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

   if ~ constructors_loaded(),
      error('Environment:Incomplete', ...
            'Function %s depends on module ''ad-props'' being active.', ...
            mfilename);
   end

   % The fluids in this simulation model are incompressible and immiscible
   % with constant viscosities.  This means we can use MRST's special
   % purpose fluid constructor |initSimpleADIFluid| to create the fluid
   % object.  We will however need to use sampled relative permeability
   % curves so we do not enter any relative permeability data in this call.

   fluid = initSimpleADIFluid('mu'    , [  1, 0.01]*centi*poise     , ...
                              'rho'   , [700, 1   ]*kilogram/meter^3, ...
                              'phases', 'OG');

   % Replace the synthetic relative permeability curves created through
   % function |initSimpleADIFluid| with the real benchmark values.

   fluid_kr = initDeckADIFluid(getSPE10_model_1_relperm());

   fluid.krG = fluid_kr.krG;
   fluid.krO = fluid_kr.krOG;
end

%--------------------------------------------------------------------------

function tf = constructors_loaded()
   tf = all(cellfun(@(fn) exist(fn, 'file') == 2, ...
                    { 'initSimpleADIFluid', 'initDeckADIFluid' }));
end
