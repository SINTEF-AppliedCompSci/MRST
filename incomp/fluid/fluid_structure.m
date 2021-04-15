%Fluid structure used in MATLAB Reservoir Simulation Toolbox.
%
% MRST's fluid representation is a structure of function handles which
% compute fluid properties (e.g., density or viscosity), fluid phase
% saturations, and phase relative permeability curves.  This representation
% supports generic implementation of derived quantities such as mobilities
% or flux functions for a wide range of fluid models.
%
% Specifically, the 'fluid' structure contains the following fields:
%
%   - properties -- Function handle for evaluating fluid properties such
%                   as density or viscosity.  Must support the syntaxen
%
%                        mu              = properties()
%                       [mu, rho]        = properties()
%                       [...]            = properties(state)
%                       [mu, rho, extra] = properties(state)
%
%                   The first call computes the phase viscosity, the
%                   second adds phase density, and the final computes
%                   additional parameters which may be needed in complex
%                   fluid models (e.g., Black Oil).
%
%                   If there are multiple PVT/property regions in the
%                   reservoir model--for instance defined by means of the
%                   ECLIPSE keyword 'PVTNUM'--then 'properties' must be a
%                   cell array of function handles with one element for
%                   each property region.  In this case each of the
%                   functions in 'properties' must support a call syntax of
%                   the form
%
%                       [output] = properties(state, i)
%
%                   where [output] represents one of the three output forms
%                   specified above and 'i' is an indicator (a LOGICAL
%                   array).  The call 'properties(state,i)' must compute
%                   property values only for those elements where 'i' is
%                   TRUE.
%
%   - saturation -- Function handle for evaluating fluid phase saturations.
%                   Must support the syntaxen
%
%                        s = saturation(state)
%                        s = saturation(state, extra)
%
%                   The first call computes fluid phase saturations based
%                   only on the current reservoir (and well) state, while
%                   the second call may utilize additional parameters and
%                   properties--typically derived from the third syntax of
%                   function 'properties'.
%
%                   If there are multiple PVT/property regions in the
%                   reservoir model--for instance defined by means of the
%                   ECLIPSE keyword 'PVTNUM'--then 'saturations' must be a
%                   cell array of function handles with one element for
%                   each property region.  In this case each of the
%                   functions in 'saturation' must support a call syntax of
%                   the form
%
%                       s = saturation(state, i)
%                       s = saturation(state, i, extra)
%
%                   where 'i' is an indicator (a LOGICAL array).  The
%                   'saturation' function must compute values only for
%                   those elements where 'i' is TRUE.
%
%   - relperm    -- Function handle, or a cell array of function handles
%                   (one element for each relative permeability region,
%                   i.e., each saturation function region) for evaluating
%                   relative permeability curves for each fluid phase.
%                   Must support the syntaxen
%
%                        kr             = relperm(s, state)
%                       [kr, dkr]       = relperm(s, state)
%                       [kr, dkr, d2kr] = relperm(s, state) [optional]
%
%                   The first call computes the relative permeability
%                   curves at the current saturations, 's'.  The second
%                   call additionally computes the first partial
%                   derivatives of these curves with respect to each phase
%                   saturation.  Finally, the third syntax additionally
%                   computes second order derivatives.  This syntax is
%                   optional and only needed in a few specific instances
%                   (e.g., when using the second order adjoint method).
%
% SEE ALSO:
%   In case the description scrolled off screen too quickly, you may access
%   the information at your own pace using the command
%
%               more on, help fluid_structure, more off

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

