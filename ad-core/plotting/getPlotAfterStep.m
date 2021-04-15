function fn = getPlotAfterStep(state0, model, schedule, varargin)
% Get a function that allows for dynamic plotting in `simulateScheduleAD`.
%
% SYNOPSIS:
%   fn = getPlotAfterStep(state0, model, schedule, 'plotWell', true);
%
% DESCRIPTION:
%   The `simulateScheduleAD` function has a optional input argument
%   `afterStepFn` that allows for dynamic plotting after each step in the
%   simulation, for instance to show how the well curves progress during the
%   simulation, or to print out extra information to the command window. This
%   function is an implementation of one such function, that can add both a
%   panel showing the simulation progress, as well as interactive plots for
%   well and reservoir quantities.
%
% PARAMETERS:
%   state0 -   Initial state for simulateScheduleAD
%
%   model -    Simulation model which will be passed to simulateScheduleAD.
%
%   schedule - The simulation schedule containing wells, driving forces
%              and time-steps that will be passed to simulateScheduleAD.
%
% KEYWORD ARGUMENTS:
%
%  'plotWell' -      Launch interactive plotting for well quantities
%                    using `plotWellSols`
%
%  'plotReservoir' - Add an interactive plotting window for reservoir
%                    quantities during the simulation. Note that, due to
%                    limitations in the implementation, this window will
%                    only be truly interactive after the simulation
%                    finishes. You can, however, set the options (field for
%                    plotting, locked color axis and so on) before
%                    initiating the simulation itself.
%
%  'view' -          View angle for the reservoir plotting. See Matlab
%                    builtin `view` for more information. Defaults to
%                    empty for no modification to the default.
%
%  'wells' -         Wells for the reservoir plotting (using `plotWell`) 
%
% RETURNS:
%  fn -              Function handle suitable for the `afterStepFn`
%                    input in `simulateScheduleAD`. 
%
% EXAMPLE: 
%  fn = getPlotAfterStep(state0, model, schedule, 'plotWell', true);
%  simulateScheduleAD(state0, model, schedule, 'afterStepFn', fn);
%
% SEE ALSO:
%   `simulateScheduleAD`, `plotWellSols`

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
    opt = struct('plotWell',        true, ...
                 'plotReservoir',   true, ...
                 'view',            [],   ...
                 'wells',           [],   ...
                 'useTimesteps',    false, ...
                 'plotWellSolArgs', {{}});
    [opt, extra] = merge_options(opt, varargin{:});
    
    G = model.G;
    
    W = schedule.control(1).W;
    
    if ~isempty(W) && opt.plotWell
        ws = initWellSolAD(schedule.control(1).W, model, state0);
        nc = numel(vertcat(W.cells));
        d = ones(nc, 1);
        sources.phaseVolume = {d, d, d};
        sources.phaseVolume = sources.phaseVolume(model.getActivePhases());
        
        model.FacilityModel = FacilityModel(model);
        model.FacilityModel = model.FacilityModel.setupWells(W);
        ws = model.FacilityModel.setWellSolStatistics(ws, sources);

        ws0 = {ws; ws};
        if opt.useTimesteps
            [hwell, injectWell] = plotWellSols(ws0, schedule.step.val([1,1]), opt.plotWellSolArgs{:});
        else
            [hwell, injectWell] = plotWellSols(ws0, opt.plotWellSolArgs{:});
        end
    else
        [hwell, injectWell] = deal(nan);
    end
    
    if opt.plotReservoir
        mrstModule add mrst-gui
        hdata = figure;
        [~, injData] = plotToolbar(G, {state0; state0}, extra{:});
        axis tight
        if ~isempty(opt.view)
            view(opt.view)
        end
        if ~isempty(opt.wells)
           plotWell(G, opt.wells)
        end
    else
        [hdata, injData] = deal(nan);
    end
    hdata = double(hdata);
    hwell = double(hwell);
    fn = @(model, states, reports, solver, schedule, simtime) afterStepFunction(model, states, reports, solver, schedule, simtime, injData, injectWell, hdata, hwell);
end

function [model, states, reports, solver, ok] = afterStepFunction(model, states, reports, solver, schedule, simtime, injData, injectWell, hdata, hwell)
    computed = cellfun(@(x) ~isempty(x), states);
    ctrl_computed = cellfun(@(x) ~isempty(x), reports);
    
    current = find(computed, 1, 'last');
    
    st = states(computed);   
    rep = reports(ctrl_computed);
    simtime = simtime(ctrl_computed);
    if ~isnan(hdata) && ishandle(hdata)
        injData(model.G, st, current);
    end
    
    if ~isnan(hwell) && ishandle(hwell)
        set(0, 'CurrentFigure', hwell);
        ws = cellfun(@(x) x.wellSol, st, 'uniformoutput', false);
        % Note: We are not actually considering the case where ministeps
        % are being inputed, which would result in an error
        T = getTimesteps(rep);
        injectWell({ws}, T)
    end
    
    ok = true;
    if 1
        ok = ok & simulationRuntimePanel(model, states, reports, solver, schedule, simtime);
    end
end

function T = getTimesteps(reports)
    T = [];
    for i = 1:numel(reports)
        for j = 1:numel(reports{i}.StepReports)
            r = reports{i}.StepReports{j};
            if r.Converged
                T = [T; r.Timestep]; %#ok<AGROW>
            end
        end
    end
    T = cumsum(T);
end