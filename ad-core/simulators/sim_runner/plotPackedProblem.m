function [wh, gh] = plotPackedProblem(problem, varargin)
% Plot simulation results from a packed problem
%
% SYNOPSIS:
%   plotPackedProblem(problem)
%   plotPackedProblem(problem)
%
% REQUIRED PARAMETERS:
%   problem - packedSimulatorProblem
%
% OPTIONAL PARAMETERS:
%   plotWellSols - Boolean indicating if well sols are to be plotted
%
%   plotStates   - Boolean indicating if states are to be plotted
%
%   Other  - Various parameters controlling plotting
%
% RETURNS:
%   wh     - Handle to plotWellSols figure
%
%   gh     - Handle to plotToolbar patch
% EXAMPLE:
%   demoPackedProblems, simulateSPE1, simulateSPE9
%
% SEE ALSO:
%   packSimulationProblem, initEclipsePackedProblemAD

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
    mrstModule add mrst-gui
    warg = {'color',    [1.0, 0.2, 0.2], ...
            'color2',   [0.2, 0.2, 1.0], ...
            'fontsize', 12};
    if iscell(problem)
        p = problem{1};
    else
        p = problem;
    end
    opt = struct('readFromDisk', false,  ...
                 'plotStates',   true,   ...
                 'view',         [],     ...
                 'name',         p.BaseName, ...
                 'wellArg',      {warg}, ...
                 'gridArg',      {{}},   ...
                 'wsArg',        {{}},   ...
                 'time',         [],     ...
                 'plotWell',     isfield(p.SimulatorSetup.model.G, 'nodes'),   ...
                 'plotWellSols', true);
    opt = merge_options(opt, varargin{:});
    if iscell(problem)
        if opt.plotWellSols
            wh = plotWellSols(problem, opt.wsArg{:});
        else
            wh = nan;
        end
        if opt.plotStates
            np = numel(problem);
            gh = nan(np, 1);
            for i = 1:np
                p = problem{i};
                nstr = sprintf('%s (%s)', p.BaseName, p.Name);
                gh(i) = plotPackedProblem(p, varargin{:}, 'plotWellSols', false, 'name', nstr);
            end
        else
            gh = nan;
        end
    else
        if isempty(opt.time)
            opt.time = cumsum(problem.SimulatorSetup.schedule.step.val);
        end
        if opt.plotStates
            gh = plotReservoir(problem, opt);
        else
            gh = nan;
        end

        if opt.plotWellSols
            ws = problem.OutputHandlers.wellSols;
            wsd = reshape(ws(:), [], 1);
            wh = plotWellSols(wsd, opt.time(1:numel(wsd)), opt.wsArg{:});
            set(gcf, 'Name', opt.name);
        else
            wh = nan;
        end
    end
end

function h = plotReservoir(problem, opt)
    setup = problem.SimulatorSetup;
    G = setup.model.G;
    states = problem.OutputHandlers.states;
    if opt.readFromDisk
        states = states(:);
    end
    figure;
    h = plotToolbar(G, states,...
        'dynamicTitle', true, ...
        'time', opt.time, ...
        'title', opt.name, ...
        opt.gridArg{:});
    if isempty(opt.view)
        if G.griddim == 3
            opt.view = [45, 45];
        else
            opt.view = [0, 90];
        end
    end
    view(opt.view(1), opt.view(2));
    axis tight
    ctrl = setup.schedule.control(1);
    if isfield(ctrl, 'W') && opt.plotWell
        plotWell(G, ctrl.W, opt.wellArg{:});
    end
end