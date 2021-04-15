%% A sloped example with hybrid VE
% In this example we build a sloping grid with internal sealing faces which
% are distributed throughout the model. We compare results for finescale,
% hybrid VE, coarse scale and simple (non-hybrid) VE simulations.
%
% For more information see:
% Møyner, O., & Nilsen, H. M. (2019). Multiresolution coupled vertical
% equilibrium model for fast flexible simulation of CO2 storage.
% Computational Geosciences, 23(1), 1-20.

%% Add modules
gravity reset on;
mrstModule add ad-core ad-blackoil ad-props co2lab matlab_bgl coarsegrid;
mrstModule add mrst-gui

%% Setup case
% First we set up the finescale simulation. This is moved to a separate
% function here to make the code less cluttered. 

% The function returns nearWell containing the cell index of cells near the
% well which will be kept as fine scale cells in the VE simulation. 

% The setToZero variable contains the face indices of sealing faces in the 
% model where transmissiblity is set to zero.

[state0, model, schedule, nearWell, setToZero] = setupSloped();

% Define some parameters for plotting
G = model.G;
W = schedule.control(1).W;
bc = schedule.control(1).bc;
T = model.operators.T_all;
da = 0.25;
xl = [-50, 1050];
yl = [-50, 110];

% Plot fine scale grid
figure();
plotOutlinedGrid(G, W, bc, T);
plotGrid(G, 'edgealpha', 0.2, 'facecolor', 'none')
view(0, 0)
axis equal tight
daspect([1, 1, da])
xlim(xl)
zlim(yl)
xlabel('Vertical position [m]');
ylabel('Depth [m]');
title('Fine scale grid with sealing layers')


% Setup nonlinear solver and pack simulation problem
nls = NonLinearSolver('maxIterations', 50);
problem = packSimulationProblem(state0, model, schedule, ...
    'sloped_finescale', 'NonLinearSolver', nls);
 
% Simulate and get the output
simulatePackedProblem(problem);
[ws, states, report] = getPackedSimulatorOutput(problem);

%% Convert to VE model

[ii, jj, kk] = gridLogicalIndices(G);
isFine = nearWell;

% Set up model
[model_ve, model_c] = convertToMultiVEModel(model, isFine, 'sealingFaces', find(model.operators.T_all == 0), 'sumTrans', true);
schedule_ve = upscaleSchedule(model_ve, schedule);
state0_ve = upscaleState(model_ve, model, state0);

% Plot grid 
figure();
plotOutlinedGrid(G, W, bc, T);
plotGrid(model_ve.G, 'edgealpha', 0.2, 'facecolor', 'none')
view(0, 0)
axis equal tight
daspect([1, 1, da])
xlim(xl)
zlim(yl)
title('Grid for hybrid VE model')

% Plot VE partitions on grid
% Here we can see how VE cells exist between sealing layers and can be
% stacked on top of each other.
figure();
plotOutlinedGrid(G, W, bc, T);
plotCellData(G, model_ve.G.partition,'edgealpha', 0.2)
view(0, 0)
axis equal tight
daspect([1, 1, da])
xlim(xl)
zlim(yl)
xlabel('Vertical position [m]');
ylabel('Depth [m]');
title('Partition of hybrid VE model')

nls = NonLinearSolver('maxIterations', 50);
nls.useRelaxation = true;
problem_ve = packSimulationProblem(state0_ve, model_ve, schedule_ve, ...
    'sloped_ve', 'NonLinearSolver', nls);

%% Simulate model
simulatePackedProblem(problem_ve);
%% Get the output
[ws_ve, states_ve, report_ve] = getPackedSimulatorOutput(problem_ve);
states_f = convertMultiVEStates(model_ve, states_ve);
%% Simple VE.
% In this model we have a simple VE simulation with no sealing faces.
model_ves = convertToMultiVEModel(model);
schedule_ves = upscaleSchedule(model_ves, schedule);
state0_ves = upscaleState(model_ves, model, state0);

problem_ves = packSimulationProblem(state0_ves, model_ves, schedule_ves, ...
    'sloped_ves', 'NonLinearSolver', nls);

[ok, status] = simulatePackedProblem(problem_ves);

[ws_ves, states_ves, report_ves] = getPackedSimulatorOutput(problem_ves);
states_fs = convertMultiVEStates(model_ves, states_ves);

%% Plot results for the finescale, hybrid VE and coarse simulations
% We loop through results, evenly spaced through the first 101 states, and
% plot the CO2 saturation for teh finescale, hybrid VE and coarse 
% simulations.

close all
fafa = find(setToZero);

parg = {3, 1};
for i = 1:10:101
    figure(1); clf
    subplot(parg{:}, 1)
    plotCellData(model.G, states{i}.s(:, 2), 'edgec', 'none');
    plotFaces(G, fafa, 'facec', 'w', 'linewidth', 2)
    view(0, 0);
    axis tight off
    title('Fine-scale saturation')
    
    subplot(parg{:}, 2)
    plotCellData(model.G, states_f{i}.s(:, 2), 'edgec', 'none');
    plotFaces(G, fafa, 'facec', 'w', 'linewidth', 2)
    view(0, 0);
    axis tight off
    title('VE reconstructed saturation')
    
    subplot(parg{:}, 3)
    plotCellData(model_ve.G, states_ve{i}.s(:, 2), 'edgec', 'none')
    plotGrid(model_ve.G, 'facec', 'none', 'edgec', 'w', 'edgea', .3)
    view(0, 0);
    axis tight off
    title(['Coarse saturation', num2str(i)])
    drawnow
end

%% Plot final results 
% CO2 saturation after injection and migration for all four 
% simulations
close all
nstep = numel(schedule.step.val);
end_inj = find(schedule.step.control == 1, 1, 'last');
end_mig = nstep;

substeps = [end_inj, end_mig];
names = {'injected', 'migrated'};
c1 = [48, 37, 255]/255;
c2 = [0, 255, 0]/255;
cc = interp1([0; 1], [c1; c2], (0:0.01:1)');

for i = 1:numel(substeps)
    ss = substeps(i);
    figure;
    for j = 1:4
        if j == 1
            g = G;
            s = states{ss}.s(:, 2);
            nm = 'Fine-scale';
        elseif j == 2
            g = G;
            s = states_f{ss}.s(:, 2);
            nm = 'Hybrid VE';
        elseif j == 3
            g = model_ve.G;
            s = states_ve{ss}.s(:, 2);
            nm = 'Coarse';
        else
            g = G;
            s = states_fs{ss}.s(:, 2);
            nm = 'VE (no layers)';
        end
        subplot(2, 2, j);
        plotCellData(g, s, 'EdgeColor', 'none')
        plotOutlinedGrid(G, W, bc, T);

        view(0, 0)
        axis equal tight
        daspect([1, 1, da])
        xlim(xl)
        zlim(yl)
        colormap(cc)
        title(nm);
    end
    set(gcf, 'Name', names{i});
end

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
