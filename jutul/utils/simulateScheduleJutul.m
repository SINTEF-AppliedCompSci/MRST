function [ws, states, reports] = simulateScheduleJutul(state0, model, schedule, varargin)
% Simulate a MRST case with Jutul, in the same manner as simulateScheduleAD
%
% SYNOPSIS:
%   [ws, states, reports] = simulateScheduleJutul(state0, model, schedule)
%
% REQUIRED PARAMETERS:
%   state0, model, schedule - Initial state, model and schedule to write.
%                             Same inputs as simulateScheduleAD.
%
%
% OPTIONAL PARAMETERS:
%   daemon - Use daemon mode to run the simulation, otherwise a case will
%            be written and a manual run must be performed. Defaults to
%            true.
%
%   pause  - If daemon = false, pause after the case has been written to
%            make manual simulation easier.
%
%   path - Valid folder path that will be used to write the output mat
%          file. By default, the Jutul output will be written to this
%          folder as well.
%
%   printcmd - Generate and print the command required to run the case in
%              Jutul and produce MRST-compatible output files. Default:
%              true.
%
%   extra - If provided, whatever data is given will be written to the
%           Jutul input .mat file. Useful if you want to pass along more
%           data for a workflow inside Jutul.
%
%   name  - Name that will be used for output file.
%
%   other - Other arguments will be passed onto either runJutulOnDaemon or
%           writeJutulInput.
%
% RETURNS:
%   pth - The path of the output file. Will be a normal .mat file.
%
% EXAMPLE:
% mrstModule add test-suite jutul
% setup = qfs_wo();
% writeJutulInput(setup.state0, setup.model, setup.schedule)
%
% SEE ALSO:
%   readJutulOutput, runJutulOnDaemon, simulateScheduleJutul

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

    opt = struct('name', class(model), 'daemon', true, 'pause', true, 'printcmd', true);
    [opt, extra] = merge_options(opt, varargin{:});
    if opt.daemon
        [ws, states] = runJutulOnDaemon(state0, model, schedule, 'name', opt.name, extra{:});
    else
        jpth = writeJutulInput(state0, model, schedule, opt.name, 'printcmd', opt.printcmd);
        if opt.pause
            disp('Pausing. Hit any key to continue once simulation has been run.')
            pause()
        end
        [ws, states] = readJutulOutput(jpth, 'error', true);
    end
    % Not parsed at the moment.
    reports = struct();
end
