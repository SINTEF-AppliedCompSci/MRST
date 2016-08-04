function command = simcommand(dir, sim)

% Construct a simulator command depending on user and context.
%
% SYNOPSIS:
%   command = simcommand(context, sim, user)
%
% PARAMETERS:
%   context - can be either 'opm' or 'urc'
%   sim     - the filename of the simulator you want to run,
%             including necessary paths from the root, for
%             example 'examples/sim_2p_comp_reorder'
%   user    - will be passed to opm_dir() or urc_dir()
%             depending on context
%
% RETURNS:
%   A string with the full path to the requested simulator program.

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

% This seems necessary on some Linux variants, and does no harm on Mac OS.
ldfix='export LD_LIBRARY_PATH=/lib/x86_64-linux-gnu/:/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH; ';

command = [ldfix, fullfile(dir, sim)];
