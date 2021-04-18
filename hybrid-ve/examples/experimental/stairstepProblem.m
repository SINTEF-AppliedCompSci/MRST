%% A problem where VE fails: Impermeable shales
gravity reset on;
mrstModule add ad-core;
mrstModule add ad-blackoil ad-props co2lab;
mrstModule add matlab_bgl coarsegrid;
% Parameters
nx = 50;
ny = 1;
nz = 100;
whdist = 5;
wvdist = 13;

L = 1000;
G = cartGrid([nx, ny, nz], [L, 1, 100]);
G = computeGeometry(G);

rock = makeRock(G, 300*milli*darcy, 0.3);

[ii, jj, kk] = gridLogicalIndices(G);

figure(1);
clf
plotGrid(G, 'facec', 'none');
view(0, 0)

K = G.cartDims(3);
% Define "shale" layers that impede flow in vertical direction and violate
% the VE assumption
setToZero = false(G.faces.num, 1);
ranges = [400, 600, ceil(7/8*K); ...
          200, 500, ceil(3*K/4); ...
          550, 800, ceil(6.5*K/8); ...
          700, 900, ceil(4*K/8); ...
          400, 850, ceil(5*K/8); ...
          0, 225, ceil(4*K/8); ...
          150, 670, ceil(3*K/8); ...
          600, 1000, ceil(2*K/8); ...
          110, 300, ceil(1*K/8); ...
          ];
for i = 1:size(ranges, 1)
    faces = addSealingFaces(G, 'x_range', ranges(i, 1:2), 'k_range', [ranges(i, 3)-1, ranges(i, 3)]);
    setToZero(faces) = true;
end
G.nodes.coords(:, 3) = G.nodes.coords(:, 3) + 2*sin(4*pi*G.nodes.coords(:, 1)/L);
G = computeGeometry(G);


% time = 10*year;
time = 5*year;

pv = poreVolume(G, rock);
inj_rate = sum(pv)/time;
nt = 100;

% Specify well information
W = [];
W = verticalWell(W, G, rock, ceil(nx/2), 1, nz, ...
            'type', 'rate', ...  % inject at constant rate
            'val', inj_rate, ... % volumetric injection rate
            'comp_i', [0 1]);    % inject CO2, not water
        
W = verticalWell(W, G, rock, 1, 1, 1, ...
            'type', 'bhp', ...  
            'val', 100*barsa, ... 
            'comp_i', [0 1]);

W = verticalWell(W, G, rock, nx, 1, 1, ...
            'type', 'bhp', ...  
            'val', 100*barsa, ... 
            'comp_i', [0 1]);
        
nearWell = false(G.cells.num, 1);
for i = 1:numel(W)
    c = W(i).cells(1);
    hdist = abs(ii - ii(c));
    vdist = abs(kk - kk(c));
    
    nearWell(hdist < whdist & vdist < wvdist) = true;
end
figure(1);
clf
plotGrid(G, 'facec', 'none');
plotFaces(G, find(setToZero), 'facec', 'r', 'linewidth', 2)
view(0, 0)
plotWell(G, W);
plotGrid(G, nearWell)

%% Define flow problem
r = 0.2;
fluid = initSimpleADIFluid('rho', [1000, 100], 'n', [1, 1], 'phases', 'WG', 'mu', [1, r]*centi*poise);
schedule = simpleSchedule(rampupTimesteps(time, time/nt), 'W', W);
% Set impermeable layers
T = getFaceTransmissibility(G, rock);
T(setToZero) = 0;
model = twoPhaseGasWaterModel(G, rock, fluid, 1, 1, 'useCNVConvergence', true);
model.operators.T = T(model.operators.internalConn);
model.operators.T_all(model.operators.internalConn) = model.operators.T;
model.extraStateOutput = true;
state0 = initResSol(G, 0, [1, 0]);

%% Simulate
nls = NonLinearSolver('maxIterations', 50);
nls.useRelaxation = true;
[ws, states] = simulateScheduleAD(state0, model, schedule, 'nonlinearsolver', nls);

%% Convert to VE model
% Set cells close to well as fine-scale
isFine = nearWell;
% Convert model
[model_ve, model_c] = convertToMultiVEModel(model, isFine);
% Just a regular coarse model with some extra information - we use the
% standard tools to upscale initial state & schedule
schedule_ve = upscaleSchedule(model_ve, schedule);
state0_ve = upscaleState(model_ve, model, state0);
%% Simulate and reconstruct
[ws_ve, states_ve] = simulateScheduleAD(state0_ve, model_ve, schedule_ve, 'nonlinearsolver', nls);
% Convert back to fine-scale by performing pressure and saturation
% reconstruction
states_f = convertMultiVEStates(model_ve, states_ve);
%% Plot saturations
fafa = find(setToZero);
for i = 1:numel(states)
    figure(1); clf
    subplot(3, 1, 1)
    plotCellData(model.G, states{i}.s(:, 2), 'edgec', 'none');
    plotFaces(G, fafa, 'facec', 'w', 'linewidth', 2)
    view(0, 0);
    axis tight off
    title('Fine-scale saturation')
    
    subplot(3, 1, 2)
    plotCellData(model.G, states_f{i}.s(:, 2), 'edgec', 'none');
    plotFaces(G, fafa, 'facec', 'w', 'linewidth', 2)
    view(0, 0);
    axis tight off
    title('VE reconstructed saturation')
    
    subplot(3, 1, 3)
    plotCellData(model_ve.G, states_ve{i}.s(:, 2), 'edgec', 'none')
    plotGrid(model_ve.G, 'facec', 'none', 'edgec', 'w', 'edgea', .3)
    view(0, 0);
    axis tight off
    title('Coarse saturation')
    drawnow
    colormap parula
    drawnow
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
