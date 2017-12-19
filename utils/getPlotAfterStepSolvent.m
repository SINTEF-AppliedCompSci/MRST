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
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

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
    
    G = model.G;
    
    df = get(0, 'defaultFigurePosition');
    resfig = figure('position', [df(1:2), 800, 400]);
    subplot(1,2,1);
    ssh = plotCellData(G, state0.s(:,4), 'edgeColor', 'none');
    caxis([0,1]); colorbar; axis equal tight
    subplot(1,2,2);
    soh = plotCellData(G, state0.s(:,2), 'edgeColor', 'none');
    caxis([0,1]); colorbar; axis equal tight;
    colormap(jet);
    
    fn = @(model, states, reports, solver, schedule, simtime) afterStepFunction(model, states, reports, solver, schedule, simtime, resfig, ssh, soh);
end

function [model, states, reports, solver, ok] = afterStepFunction(model, states, reports, solver, schedule, simtime, resfig, ssh, soh)
    computed = cellfun(@(x) ~isempty(x), states);
    ctrl_computed = cellfun(@(x) ~isempty(x), reports);
    
    current = find(computed, 1, 'last');
    
    st = states(computed);   
    rep = reports(ctrl_computed);
    
    current = find(computed, 1, 'last');
    
    simtime = simtime(ctrl_computed);
    [~, bc] = boundaryFaces(model.G);
    
    set(0, 'CurrentFigure', resfig);
    ssh.CData = st{current}.s(bc,4);
    soh.CData = st{current}.s(bc,2);
    
    
    ok = true;
    if 1
        ok = ok & simulationRuntimePanel(model, states, reports, solver, schedule, simtime);
    end
end