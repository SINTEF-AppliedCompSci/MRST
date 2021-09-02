function fn = getPlotAfterStepMICP(state0, model, az, el)
% Get a function that allows for dynamic plotting in `simulateScheduleAD`
% for the micp model. The parameters az and el are the azimuth and 
% elevation angles for view of the current axes.
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
    sz = size(c, 1);
    n1 = subplot(3, 2, 1);
    sm = plotCellData(G, state0.m, 'edgeColor', 'none');
    title('Microbes [kg/m$^3$]', 'Interpreter', 'latex');
    colorbar; view(az, el);
    colormap(n1, c((round( 70 * sz / 256)) : (round(120 * sz / 256)), :));
    n2 = subplot(3, 2, 2);
    so = plotCellData(G, state0.o, 'edgeColor', 'none');
    title('Oxygen [kg/m$^3$]', 'Interpreter', 'latex');
    colorbar; view(az, el);
    colormap(n2, c((round( 70 * sz / 256)) : (round(180 * sz / 256)), :));
    caxis([0 model.fluid.Comax]);
    n3 = subplot(3, 2, 3);
    sb = plotCellData(G, state0.b, 'edgeColor', 'none');
    title('Biofilm [-]', 'Interpreter', 'latex');
    colorbar; view(az, el);
    colormap(n3, c((round( 70 * sz / 256)) : (round(120 * sz / 256)), :));
    n4 = subplot(3, 2, 4);
    su = plotCellData(G, state0.u, 'edgeColor', 'none');
    title('Urea [kg/m$^3$]', 'Interpreter', 'latex');
    colorbar; view(az, el);
    colormap(n4, c((round(70 * sz / 256)) : (round(100 * sz / 256)), :));
    caxis([0 model.fluid.Cumax]);
    n5 = subplot(3, 2, 5); 
    sc = plotCellData(G, state0.c, 'edgeColor', 'none');
    title('Calcite [-]', 'Interpreter', 'latex');
    colorbar; view(az, el);  
    colormap(n5, c((round(70 * sz / 256)) : end, :));
    caxis([0 model.fluid.crit]);
    n6 = subplot(3, 2, 6);
    sk = plotCellData(G, 100 * (1 - model.fluid.K(model.rock.poro - ...
            state0.b - state0.c) ./ model.rock.perm), 'edgeColor', 'none');
    title('Permeability [\%]', 'Interpreter', 'latex');
    colorbar; view(az, el);
    colormap(n6, c((round(70 * sz / 256)) : end, :));
    caxis([0 100]);
    
    fn = @(model, states, reports, solver, schedule, simtime) ... 
            afterStepFunction(model, states, reports, solver, schedule, ...
                                  simtime, resfig, sm, so, sb, su, sc, sk);
end

function [model, states, reports, solver, ok] = afterStepFunction( ...
 model, states, reports, solver, schedule, simtime, resfig, sm, so, sb, ...
                                                                su, sc, sk)
    computed = cellfun(@(x) ~isempty(x), states);
    ctrl_computed = cellfun(@(x) ~isempty(x), reports);
    
    st = states(computed);
    
    current = find(computed, 1, 'last');
    
    simtime = simtime(ctrl_computed);
    [~, bc] = boundaryFaces(model.G);
    
    set(0, 'CurrentFigure', resfig);
    sm.CData = st{current}.m(bc, 1);
    so.CData = st{current}.o(bc, 1);
    sb.CData = st{current}.b(bc, 1);
    su.CData = st{current}.u(bc, 1);
    sc.CData = st{current}.c(bc, 1);
    st{current}.k = 100 * (1 - model.fluid.K(model.rock.poro - ...
                        st{current}.b - st{current}.c) ./ model.rock.perm);
    sk.CData = st{current}.k(bc, 1);
    
    ok = true;
    if 1
        ok = ok & simulationRuntimePanel(model, states, reports, ...
                                                solver, schedule, simtime);
    end
    
    % The simulator stops if clogging has been reached in any of the cells 
    fn = checkCloggingMICP(ok);
    [model, states, reports, solver, ok] = fn(model, states, reports, ...
                                                solver, schedule, simtime);
end