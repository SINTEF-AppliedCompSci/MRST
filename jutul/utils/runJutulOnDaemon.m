function [ws, states] = runJutulOnDaemon(state0, model, schedule, varargin)
% Simulate a MRST case with Jutul using the DaemonMode module.
%
% DESCRIPTION:
%   This function calls a background Julia session to simulate cases,
%   before the output is retrieved back via a .mat file. The benefit of
%   this approach is that Jutul can be used as a direct accelerator for
%   MRST, with no manual management of the Julia session once the Daemon
%   has been started.
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
%   name  - Name for the case (used to create output folders only)
%
%   project - Path to the Julia enviroment where DaemonMode is installed.
%             Not required if DaemonMode is present in the base ENV. See
%             backgroundJutulExample for more details.
%
%   path - Path to folder where input and output files are written.
%
%   julia_path - Path to julia executable. If not provided, MRST will
%                assume that Julia is available directly on path in the
%                system() call.
%
% NOTE:
%    Additional inputs are passed as keyword arguments to Julia simulator.
%    First simulation can take some time, as this triggers compilation.
%    DaemonMode appears to be somewhat slower than a direct run, but still
%    fast relative to MRST for most cases. There is significant I/O
%    overhead for very small cases.
%
% RETURNS:
%   ws, states - Reformatted output that appears similar to MRST's.
%
% EXAMPLE:
%    See backgroundJutulExample for more details.
%
% SEE ALSO:
%   readJutulOutput, simulateScheduleJutul

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
    opt = struct('name', 'jutul_case', ...
                 'project', '', ...
                 'julia_path', 'julia', ...
                 'path', tempdir());
    [opt, extra] = merge_options(opt, varargin{:});
    v = mrstVerbose();
    if ~isempty(opt.project)
        opt.project = sprintf('--project="%s"', opt.project);
    end
    pth = writeJutulInput(state0, model, schedule, opt.name, 'printcmd', false, 'path', opt.path);
    % Create a julia script that runs the file
    cmd_pth = fullfile(opt.path, sprintf('run_%s.jl', opt.name));
    dispif(v, 'Creating Julia runscript at %s... ', cmd_pth) 
    f = fopen(cmd_pth, 'w');
    if f == 0
        error(ferror(f));
    end
    if v
        info = 1;
    else
        info = -1;
    end
    extra_str = '';
    ne = numel(extra);
    assert(mod(ne, 2) == 0, 'Additional inputs must be key/value pairs')
    info_found = false;
    for i = 1:2:(ne-1)
        k = extra{i};
        if strcmpi(k, 'info_level')
            info_found = true;
        end
        val = extra{i+1};
        if isnumeric(val)
            val = num2str(val);
        elseif islogical(val)
            if val
                val = 'true';
            else
                val = 'false';
            end
        end
        extra_str = sprintf('%s, %s=%s', extra_str, k, val);
    end
    if ~info_found
        extra_str = sprintf('%s, info_level=%d', extra_str, info);
    end
    dispif(v, 'ok.\n');
    dispif(v, 'Running Jutul simulation...\n');
    fprintf(f, 'using JutulDarcy\nsimulate_mrst_case(\"%s\", write_mrst = true%s, ascii_terminal = true)\n', ...
        pth, extra_str);
    fclose(f);
    % Finally put together the command to invoke the daemon in client mode
    % and run the case.
    cmd = sprintf('%s %s --startup-file=no --color=no -e "using DaemonMode; runargs()" %s', ...
        opt.julia_path, opt.project, cmd_pth);
    id = system(cmd);
    if id == 1
        error('Julia simulation was unable to complete successfully.')
    end
    dispif(v, 'Julia simulation complete.\n');
    dispif(v, 'Reading Julia output... ');
    [ws, states] = readJutulOutput(pth);
    dispif(v, 'ok.\n');
end
