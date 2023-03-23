function fn = getPlotAfterStepCO2(state0, model, az, el)
% Get a function that allows for dynamic plotting in `simulateScheduleAD`
% for the assesment of CO2 in the ad-micp module. The parameters az and 
% el are the azimuth and elevation angles for view of the current axes.
% 
% This function is modified from a file in The MATLAB Reservoir Simulation
% Toolbox (MRST), see
%   mrst/modules/solvent/utils/getPlotAfterStepSolvent.m
%
% We refer to that function for a complete commented version of the file. 
% In this file we comment on some of the lines. 

%{
Partial copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
Partial copyright 2021 NORCE Norwegian Research Centre AS, Computational 
Geosciences and Modeling.

This file is part of the ad-micp module.

ad-micp is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ad-micp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%}
    
    G = model.G; 
    df = get(0, 'defaultFigurePosition');
    resfig = figure('position', [df(1 : 2), 800, 1200]);
    c = flipud(jet);
    sz = size(c,1);
    sco2 = plotCellData(G, state0.s(:, 2), 'edgeColor', 'none');
    title('CO2 [-]', 'Interpreter','latex');
    colorbar; view(az, el); 
    colormap(c((round( 70 * sz / 256)):(round(100 * sz / 256)), :));
    caxis([0 1]);
    
    fn = @(model, states, reports, solver, schedule, simtime) ... 
            afterStepFunction(model, states, reports, solver, schedule, ...
                                  simtime, resfig, sco2);
end

function [model, states, reports, solver, ok] = afterStepFunction( ...
           model, states, reports, solver, schedule, simtime, resfig, sco2)
    computed = cellfun(@(x) ~isempty(x), states);
    ctrl_computed = cellfun(@(x) ~isempty(x), reports);
    
    st = states(computed);
    
    current = find(computed, 1, 'last');
    
    simtime = simtime(ctrl_computed);
    [~, bc] = boundaryFaces(model.G);
    
    set(0, 'CurrentFigure', resfig);
    sco2.CData = st{current}.s(bc, 2);
    
    ok = true;
    if 1
        ok = ok & simulationRuntimePanel(model, states, reports, ...
                                                solver, schedule, simtime);
    end
end