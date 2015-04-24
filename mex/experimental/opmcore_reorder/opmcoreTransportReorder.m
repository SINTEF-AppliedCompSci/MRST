function state = opmcoreTransportReorder(state, G, rock, src, dt)
%Solve incompressible two-phase transport with reordering.
%
% SYNOPSIS:
%   state = opmcoreTransportReorder(state, G, rock, src, dt)
%
% DESCRIPTION:
%   This function calls C++ code from the opm-core library to perform
%   incompressible two-phase transport.
%
% REQUIRED PARAMETERS:
%   state  - Reservoir and well solution structure either properly
%            initialized from function 'initState', or the results from a
%            previous call to function 'solveIncompFlow' and, possibly, a
%            transport solver such as function 'explicitTransport'.
%            Assumed to contain the fields 'flux' for signed face fluxes
%            and 's' for saturations.
%
%   G      - Grid.
%
%   rock   - Rock data structure.
%
%   src    - Vector of source contributions, one for each cell. Positive
%            entries are taken to mean water inflow rates, while negative
%            entries are taken to mean total outflow rates.
%
%   dt     - Timestep.
%
% RETURNS:
%   state - Update reservoir solution structure with new values
%           for the field:
%              - s        -- Saturation values for all cells in the
%                            discretised reservoir model, 'G'.
%
% NOTE:
%   There is at the moment no way to provide fluid information to the
%   solver, which runs with a default, primitive fluid model (linear
%   relperm, default values for viscosities and densities). Future
%   revisions to this function will remedy this.
%
% TROUBLESHOOTING:
%   To compile the mex-function you must have the opm-core library
%   installed in a standard location (such as /usr/local). You can
%   find installation instructions for opm-core at
%
%     https://public.ict.sintef.no/openrs/wiki/Installing_opm-core
%
%   When following the instructions above, please note that if you
%   use the --prefix mechanism to choose a non-default installation
%   location, the library will probably not be found when trying to
%   compile the mex-function.
%
%   Sometimes a problem may occur with the C++ run-time library,
%   libstdc++, that is provided by Matlab for linking the
%   mex-file. This may be due to having opm-core (or one of its
%   dependent boost libraries) linked against a newer libstdc++ than
%   the one provided by Matlab. A workaround may be to replace
%   Matlab's libstdc++ with the system-provided one, but this may
%   run the risk of wrecking your Matlab installation, do this only
%   as a last resort, and do make a backup of the original
%   libstdc++ first.
%
% SEE ALSO:
%   solveIncompFlow.

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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


   pv = poreVolume(G, rock);
   % Taking the transpose in the next two lines since
   % the saturation matrix orders are different in
   % Matlab and the opm code.
   s0 = state.s';
   state.s = mex_opmcoreTransportReorder(G, pv, state.flux, src, dt, s0)';
end
