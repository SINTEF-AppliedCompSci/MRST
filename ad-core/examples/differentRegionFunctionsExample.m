%% The use of regions: Different functions in different parts of the domain
% In many applications there is a need for varying functions in different
% parts of the domain. A typical example of this is the problem of varying
% rock types, where multiphase flow happens in a medium where the surface
% properties and porous structure varies significantly in different regions
% of the domain. In this example we demonstrate how to set up multiple rock
% types for relative permeability. The purpose of this example is to show
% the function interfaces used to have varying functions. The problem is
% not intended to be realistic.
mrstModule add ad-core ad-blackoil ad-props mrst-gui
%% Set up parameters
% We create a small grid, with initially uniform fluid parameters.
dims = [20, 20];
G = cartGrid(dims, [1, 1]);
G = computeGeometry(G);
fluid = initSimpleADIFluid('phases', 'wo', 'rho', [100, 100],...
                            'cR', 1e-10/barsa, ...
                                           'mu', [1, 1]);
%% Add multiple relative permeability functions
% We create a quadratic and a linear relative-permeability relation as
% function handles of the saturation. We define the water and oil relative
% permeability to be equal as cell arrays. In doing so, the first entry
% corresponds to the first region and the second function to the second
% region and so on.
% 
% We also create a cell-wise region indicator to map the two functions we
% have defined onto the grid. The exact mechanism for how we set up the
% region in the model can vary.
kr_1 = @(S) S.^2;
kr_2 = @(S) S;
fluid.krW = {kr_1, kr_2};
fluid.krO = {kr_1, kr_2};

rock = makeRock(G, 0.1*darcy, 0.3);
reg = 1 + double(G.cells.centroids(:, 2) > 0.5);
%% Plot the fluid regions and the corresponding relative permeability curves
% We have a top and bottom part of the domain
figure;
subplot(1, 2, 1)
plotCellData(G, reg);
axis equal tight
colorbar('horiz')
colormap(lines(2));

subplot(1, 2, 2); hold on;
s = (0:0.05:1)';
nkr = numel(fluid.krW);
l = cell(1, nkr);
for i = 1:nkr
    k = fluid.krW{i};
    plot(s, k(s), 'linewidth', 2);
    l{i} = sprintf('Region %d', i);
end
legend(l, 'location', 'northwest');
xlabel('Saturation');
title('Relative permeability functions')
%% Set up simulation scenario
% We inject 1 pore-volume over 100 days, from the left part of the domain
% to the right.
pv = poreVolume(G, rock);
time = 100*day;
irate = sum(pv)/time;
bc = [];
bc = fluxside(bc, G, 'Xmin', irate, 'sat', [1, 0]);
bc = pside(bc, G, 'XMax', 100*barsa, 'sat', [1, 0]);

state0 = initResSol(G, 100*barsa, [0, 1]);
schedule = simpleSchedule(rampupTimesteps(time, time/10), 'bc', bc);
%% Approach 1: Use different rock-types
% We can set the rock.regions.saturation field to set up the different
% regions for the model. This is the default place the relative
% permeability and capillary pressure implementations will look for a
% region.
rock_reg = rock;
rock_reg.regions = struct('saturation', reg);
model = GenericBlackOilModel(G, rock_reg, fluid, 'gas', false);
[~, states] = simulateScheduleAD(state0, model, schedule);

figure;
plotToolbar(G, states, 'startPlayBack', true, 'field', 's:1')
title('Different regions (via rock)');
%% Approach 2: Directly modify the state functions
% For a more fine-grained approach, we can manually update the
% RelativePermeability state function to have a new region. If a state
% function is using a function handle from the fluid in it's implementation
% and it has a "regions" field set up, it will use that to evaluate the
% function. Since it is possible to swap out the StateFunctions to change
% the implementation of a single function, there are many ways to
% incorporate spatially varying functions.
model = GenericBlackOilModel(G, rock, fluid, 'gas', false);
model = model.validateModel();
model.FlowPropertyFunctions.RelativePermeability.regions = reg;

[~, states2] = simulateScheduleAD(state0, model, schedule);
figure;
plotToolbar(G, states2, 'startPlayBack', true, 'field', 's:1')
title('Different regions (via rock)');
%% Approach 3: Use input files
% We do not demonstrate this in this example, but if you build your model
% using inputs from the deckformat module, regions are automatically set up
% if present.

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
