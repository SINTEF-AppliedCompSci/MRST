function [ws, states, reports] = simulateDataFileJutul(path)
% Simulate a .DATA file with JutulDarcy
%
% SYNOPSIS:
%   [ws, states, reports] = simulateDataFileJutul('/path/to/file.data')
%
% REQUIRED PARAMETERS:
%   path - Path to .DATA file that is to be simulated.
%
% RETURNS:
%   pth - The path of the output from the simulation.
%
%
% SEE ALSO:
%   readJutulOutput, runJutulOnDaemon, simulateScheduleJutul

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

    fprintf('To simulate .DATA file, call:\n\n\tusing JutulDarcy\n\tsimulate_mrst_case(raw"%s", write_mrst = true);\n\n', path)

    disp('Pausing. Hit any key to continue once simulation has been run.')
    pause()
    [fldr, fname, ext] = fileparts(path); %#ok
    jpth = fullfile(fldr, fname);
    [ws, states] = readJutulOutput(jpth, 'error', true);
    disp('Read simulated case.')
    % Not parsed at the moment.
    reports = struct();
end
