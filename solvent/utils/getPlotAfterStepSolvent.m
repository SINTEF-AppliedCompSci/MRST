function fn = getPlotAfterStepSolvent(state0, model, varargin)
% Get a function that allows for dynamic plotting in `simulateScheduleAD`
% for the black-oil solvent model.
%
% SEE ALSO:
%   `getPlotAfterStep`

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
    
    st = states(computed);
    
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