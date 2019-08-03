function h = monitorBackgroundSimulations(problems, varargin)
% Monitor simulations running in the background or another session
%
% SYNOPSIS:
%   monitorBackgroundSimulations(problems)
%   monitorBackgroundSimulations(problems, 'useFigure', false)
%
% REQUIRED PARAMETERS:
%   problems  - A cell array of problems to monitor
%
% OPTIONAL PARAMETERS:
%   useFigure - Use a figure for showing progress. Default: true if GUI is
%               enabled.
%
% RETURNS:
%   Nothing


%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('pause', 0.25, 'useFigure', [], ...
                 'handle', struct('figure', [], 'text', []),...
                 'singleUpdate', false);
    opt = merge_options(opt, varargin{:});
    if ~iscell(problems)
        problems = {problems};
    end
    np = numel(problems);
    simulating = true;
    basenames = unique(cellfun(@(x) x.BaseName, problems, 'UniformOutput', false));
    bn = join(basenames, ',');
    fn = sprintf('Simulating: %s', bn{1});
    h = opt.handle;
    if isempty(opt.useFigure) || opt.useFigure
        if isempty(h.figure)
            h.figure = figure('Name', fn, 'ToolBar', 'none', 'NumberTitle', 'off', 'MenuBar', 'none');
        end
        opt.useFigure = isgraphics(h.figure);
    end
    dispstr = '';
    if ~opt.useFigure
        fprintf('\n');
        if ~isempty(h.text)
            dispstr = [h.text, '\n'];
        end
    end
    
    
    height = 0.9/np;
    descriptions = cellfun(@getDescription, problems, 'UniformOutput', false);
    while simulating
        n_prev = numel(dispstr);
        dispstr = '';
        if opt.useFigure
            clf(h.figure);
        end
        active = false(np, 1);
        for i = 1:np
            p = problems{i};
            [n, num, active(i)] = getStatus(p);
            descr = descriptions{i};
            if opt.useFigure
                simpleUIbar(h.figure, num/n, (i-1)*height, height, descr)
            else
                nextstr = sprintf('%d) %s\n-> %d of %d steps simulated (%2.2f%% done).', ...
                                  i, descr, num, n, 100*num/n);
               if i == 1
                   dispstr = nextstr;
               else
                   dispstr = [dispstr, newline, nextstr];
               end
            end
        end
        simulating = any(active);
        if ~opt.useFigure
            remstr = sprintf(repmat('\b', 1, n_prev+1));
            h.text = dispstr;
            disp([remstr, dispstr])
        end
        if opt.singleUpdate
            break
        else
            pause(opt.pause);
        end
    end
end

function [n, num, active] = getStatus(p)
    n = numel(p.SimulatorSetup.schedule.step.val);
    num = p.OutputHandlers.states.numelData();
    active = num < n;
end

function descr = getDescription(p)
    if isempty(p.Description)
        descr = p.Name;
    else
        descr = [p.Name, ': ', p.Description];
    end
end