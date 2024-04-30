%% Load modules, turn on gravity and set utility functions
mrstModule add co2lab-common co2lab-ve ad-core ad-props ad-blackoil
gravity on;

view_side = @() view(0,0);
view_oblique = @() view(-13, 2);

%%
% In this example, we compare VE and 3D simulations of injected CO2 migrating up
% a (synthetic) sloping aquifer, while looking at the effect of residual
% trapping and capillary effects.


%% Residual trapping, sharp interface
% In the first example, we include residual trapping but no capillary
% pressure in the simulation.  The result is a sharp-interface between CO2
% and water, and a trail of trapped CO2 in the wake of the moving plume.

srw = 0.25; % residual water satuation
srg = 0.25; % residual gas saturation

% Simulate cross-sections
[G, states_res, G_VE, states_VE_res, sat_VE_converted_3D_res, timing_res, wellcell] = ...
    sloping_aquifer('cross_sectional', true, 'zres', 40, 'srw', srw, 'srg', srg);

% plot grid with well cell
fig = figure;
plotGrid(G, 'facecolor', 'yellow', 'facealpha', 0.6, 'edgealpha', 0.3);
plotGrid(G, wellcell, 'facecolor', 'red');
set_fig_standardformat(fig, "Cross-sectional simulation grid, with cell containing injection point");
view_side()

% Compare 3D and reconstructed VE plots
fig = figure;
subplot(2,1,1); plotCellData(G, states_res{20}.s(:,2), 'edgealpha', 0.1); view_side(); 
subplot(2,1,2); plotCellData(G, sat_VE_converted_3D_res{20}, 'edgealpha', 0.1); view_side();
set_fig_standardformat(fig, {'CO2 saturation, timestep 20, 3D', ...
                       'CO2 saturation, timestep 20, VE'});

fig = figure;
subplot(2,1,1); plotCellData(G, states_res{end}.s(:,2), 'edgealpha', 0.1); view_side(); 
subplot(2,1,2); plotCellData(G, sat_VE_converted_3D_res{end}, 'edgealpha', 0.1); view_side();
set_fig_standardformat(fig, {'CO2 saturation, last timestep, 3D', ...
                       'CO2 saturation, last timestep, VE'});

% Note that the main difference is the region with residual trapping around
% the wellbore.  In this region, with significant vertical flow during
% injection, the VE model does not capture the effect of CO2 rising from the 
% injection point at the bottom.   For the rest of the model, the results are
% very similar.

%% Capillary fringe, no residual trapping

% In this second example, we ignore residual trapping, but impose a capillary
% pressure that will introduce a capillary transition zone between the CO2
% and brine region ("capillary fringe" rather than "sharp interface").  

% Simulate cross-sections
[G, states_fringe, G_VE, states_VE_fringe, sat_VE_converted_3D_fringe, timing_fringe, wellcell] = ...
    sloping_aquifer('cross_sectional', true, 'zres', 40, 'cap_fringe', true);

% Compare 3D and reconstructed VE plots
fig = figure;
subplot(2,1,1); plotCellData(G, states_fringe{20}.s(:,2), 'edgealpha', 0.1); view_side(); 
subplot(2,1,2); plotCellData(G, sat_VE_converted_3D_fringe{20}, 'edgealpha', 0.1); view_side();
set_fig_standardformat(fig, {'CO2 saturation, timestep 20, 3D', ...
                       'CO2 saturation, timestep 20, VE'});

fig = figure;
subplot(2,1,1); plotCellData(G, states_fringe{end}.s(:,2), 'edgealpha', 0.1); view_side(); 
subplot(2,1,2); plotCellData(G, sat_VE_converted_3D_fringe{end}, 'edgealpha', 0.1); view_side();
set_fig_standardformat(fig, {'CO2 saturation, last timestep, 3D', ...
                       'CO2 saturation, last timestep, VE'});


%% Both capillary fringe and residual trapping

% In the last example, we combine capillary pressure and residual trapping
% and once again compare the simulated migration from the VE and 3D models.
% In this case, hysteresis play an important role.  The VE model comes with
% an implementation of an endpoint-scaling model to represent hysteresis.
% The 3D model used here has no equivalent mechanism for handling hysteresis,
% leading to a significantly higher amount of residual trapping and limited
% migration distance.

% Simulate cross-sections
[G, states_fringe, G_VE, states_VE_fringe, sat_VE_converted_3D_fringe, timing_fringe, wellcell] = ...
    sloping_aquifer('cross_sectional', true, 'zres', 40, 'cap_fringe', true, ...
                    'srw', srw, 'srg', srg);

% Compare 3D and reconstructed VE plots
fig = figure;
subplot(2,1,1); plotCellData(G, states_fringe{20}.s(:,2), 'edgealpha', 0.1); view_side(); 
subplot(2,1,2); plotCellData(G, sat_VE_converted_3D_fringe{20}, 'edgealpha', 0.1); view_side();
set_fig_standardformat(fig, {'CO2 saturation, timestep 20, 3D', ...
                       'CO2 saturation, timestep 20, VE'});

fig = figure;
subplot(2,1,1); plotCellData(G, states_fringe{end}.s(:,2), 'edgealpha', 0.1); view_side(); 
subplot(2,1,2); plotCellData(G, sat_VE_converted_3D_fringe{end}, 'edgealpha', 0.1); view_side();
set_fig_standardformat(fig, {'CO2 saturation, last timestep, 3D', ...
                       'CO2 saturation, last timestep, VE'});

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.
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
