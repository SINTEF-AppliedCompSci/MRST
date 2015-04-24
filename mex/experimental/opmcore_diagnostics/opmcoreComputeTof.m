function tof = opmcoreComputeTof(state, G, rock, src)
% Solve single-phase time-of-flight equation with reordering.
%
% SYNOPSIS:
%   state = opmcoreComputeTof(state, G, rock, src, dt)
%
% DESCRIPTION:
%   This function calls C++ code from the opm-core library to perform
%   time-of-flight computations.
%
% REQUIRED PARAMETERS:
%   state  - Reservoir and well solution structure either properly
%            initialized from function 'initState', or the results from a
%            previous call to function 'solveIncompFlow' and, possibly, a
%            transport solver such as function 'explicitTransport'.
%            Assumed to contain the fields 'flux' for signed face fluxes.
%
%   G      - Grid.
%
%   rock   - Rock data structure.
%
%   src    - Vector of source contributions, one for each cell. Positive
%            entries are taken to mean water inflow rates, while negative
%            entries are taken to mean total outflow rates.
%
%
% RETURNS:
%   tof    - Time-of-flight values for all cells in the grid.
%
% TROUBLESHOOTING:
%   To compile the mex-function you must have the opm-core library
%   installed in a standard location (such as /usr/local). You can
%   find opm-core at
%
%     https://github.com/OPM/opm-core
%
%   Follow the build instructions in the README, and you should
%   have a working library. The mex-function used here makes
%   assumptions as to what options are required for headers and
%   linkage, you may have to edit mex_opmcoreTransportReorder.m to
%   fit your particular system.
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
   % the matrix orders are different in Matlab and the
   % opm code. It is most convenient to let Matlab do
   % the transformation here.
   % Also note that while the mex function takes optional
   % arguments to control discretization method (FV, DG0, DG1)
   % and velocity interpolation method (Constant, ECVI), this
   % has not been enabled since it is untested (only use FV).
   % Part of the problem is that for the DG implementation,
   % LAPACK is used, and it is hard to avoid conflicts with
   % Matlab's own versions of BLAS and LAPACK.
   tof = mex_opmcoreComputeTof(G, pv, state.flux, src)';
end
