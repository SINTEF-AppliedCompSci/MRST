%% Load modules, turn on gravity and set utility functions
mrstModule add co2lab-common co2lab-ve ad-core ad-props ad-blackoil
gravity on;

view_side = @() view(0,0);
view_oblique = @() view(-13, 2);

%% Compare VE solution with 3D solutions with increasing vertical resolution

% In this example, we demonstrate how the solution from the 3D model
% converges towards the VE solution when vertical grid resolution is
% increased.  Although the result is the same, the VE simulation runs at a
% fraction of the time used to simulate the 3D model.  We simulate on a 
% section of a synthetic, sloping aquifer model.

%% Compare cross-sections, with low vertical resolution (5)
[G5, states5, G5_VE, states5_VE, sat5_VE_converted_3D, timing5, wellcell5] = ...
    sloping_aquifer('cross_sectional', true, 'zres', 5);

% plot grid
fig = figure;
plotGrid(G5, 'edgealpha', 0.3);
set_fig_standardformat(fig, 'Cross-sectional simulation grid'); view_side(); 

% plot grid with well cell
fig = figure;
plotGrid(G5, 'facecolor', 'yellow', 'facealpha', 0.6, 'edgealpha', 1);
plotGrid(G5, wellcell5, 'facecolor', 'red');
set_fig_standardformat(fig, "Cross-sectional simulation grid, with cell containing injection point");
view_side()

% Compare 3D and reconstructed VE plots
fig = figure;
subplot(2,1,1); plotCellData(G5, states5{end}.s(:,2)); view_side(); 
subplot(2,1,2); plotCellData(G5, sat5_VE_converted_3D{end}); view_side();
set_fig_standardformat(fig, {'CO2 saturation, last timestep, 3D', ...
                       'CO2 saturation, last timestep, VE'});

%% compare cross-sections, using higher vertical resolution (15)
[G15, states15, G15_VE, states15_VE, sat15_VE_converted_3D, timing15, wellcell15] = ...
    sloping_aquifer('cross_sectional', true, 'zres', 15);

% Compare 3D and reconstructed VE plots
fig = figure;
subplot(2,1,1); plotCellData(G15, states15{end}.s(:,2), 'edgealpha', 0.2); view_side(); 
subplot(2,1,2); plotCellData(G15, sat15_VE_converted_3D{end}, 'edgealpha', 0.2); view_side();
set_fig_standardformat(fig, {'CO2 saturation, last timestep, 3D', ...
                       'CO2 saturation, last timestep, VE'});

% compare cross-sections, using high vertical resolution (40)
[G40, states40, G40_VE, states40_VE, sat40_VE_converted_3D, timing40, wellcell40] = ...
    sloping_aquifer('cross_sectional', true, 'zres', 40);

% Compare 3D and reconstructed VE plots
fig = figure;
subplot(2,1,1); plotCellData(G40, states40{end}.s(:,2), 'edgealpha', 0.2); view_side(); 
subplot(2,1,2); plotCellData(G40, sat40_VE_converted_3D{end}, 'edgealpha', 0.2); view_side();
set_fig_standardformat(fig, {'CO2 saturation, last timestep, 3D', ...
                       'CO2 saturation, last timestep, VE'});

%% compare cross-sections, for extremely high vertical resolution (140)
[G140, states140, G140_VE, states140_VE, sat140_VE_converted_3D, timing140, wellcell140] = ...
    sloping_aquifer('cross_sectional', true, 'zres', 140);

% Compare 3D and reconstructed VE plots
fig = figure;
subplot(2,1,1); plotCellData(G140, states140{end}.s(:,2), 'edgealpha', 0.2); view_side(); 
subplot(2,1,2); plotCellData(G140, sat140_VE_converted_3D{end}, 'edgealpha', 0.2); view_side();
set_fig_standardformat(fig, {'CO2 saturation, last timestep, 3D', ...
                       'CO2 saturation, last timestep, VE'});

% Compare runtimes
fprintf('Runtime for the 3D simulation: %4.2f seconds.\n', timing140.sim3D);
fprintf('Runtime for the VE simulation: %4.2f seconds.\n', timing140.simVE);

% compare curves
states3D = {states5, states15, states140};
grids3D = {G5, G15, G140};
markers = {'--', '-o', '-x'};

abscissa = G5.cells.centroids(1:G5.cartDims(1), 1); % x-coordinates of cell
                                                    % centroids of first row of
                                                    % cells along x-axis in grid
selected_tsteps = [10, 20, 31];
aquifer_thickness = 20;

fig = figure;

for tix = 1:numel(selected_tsteps)
    tstep = selected_tsteps(tix);
    subplot(3, 1, tix);
    hold on;
    for i = 1:numel(states3D)
        averaged_sat = sum(reshape(states3D{i}{tstep}.s(:,2), grids3D{i}.cartDims), 3) / grids3D{i}.cartDims(end);
        plot(abscissa, averaged_sat * aquifer_thickness, markers{i});
    end
    plot(abscissa, states5_VE{tstep}.s(:,2) * aquifer_thickness, 'k-');
    days = floor(states3D{1}{tstep}.time/day);
    title(sprintf('%d days', days))
    xlabel('x (meter)'); ylabel('Thickness (m)');
    set(gca, 'fontsize', 14);
    if tix == numel(selected_tsteps)
        % only show legend on one subplot
        legend('z resolution = 5', 'z resolution = 15', 'z resolution = 140', ...
               'vertical equilibrium');
    end
end

set(fig, 'color', 'white')
set(fig, 'position', [1000 350 1200 900]);

%% Compare in larger domain (not cross section)

% Compare cross-sections, with low vertical resolutions (5 and 15)
[G5_full, states5_full, G5_VE_full, states5_VE_full, sat5_VE_converted_3D_full, timing5_full, wellcell5_full] = ...
    sloping_aquifer('cross_sectional', false, 'zres', 5);

[G15_full, states15_full, G15_VE_full, states15_VE_full, sat15_VE_converted_3D_full, timing15_full, wellcell15_full] = ...
    sloping_aquifer('cross_sectional', false, 'zres', 15);


[G40_full, states40_full, G40_VE_full, states40_VE_full, sat40_VE_converted_3D_full, timing40_full, wellcell40_full] = ...
    sloping_aquifer('cross_sectional', false, 'zres', 40);

% plot 3D grid with well cell
fig = figure;
plotGrid(G5_full, 'facecolor', 'yellow', 'facealpha', 0.7, 'edgealpha', 0.4);
plotGrid(G5_full, wellcell5_full, 'facecolor', 'red');
set_fig_standardformat(fig, "3D simulation grid, with cell containing injection point");
view_oblique()

% Compare 3D and reconstructed VE plots
fig = figure;
subplot(2,2,1); plotCellData(G5_full, states5_full{end}.s(:,2), 'edgealpha', 0.2); view_oblique();
subplot(2,2,2); plotCellData(G15_full, states15_full{end}.s(:,2), 'edgealpha', 0.2); view_oblique(); 
subplot(2,2,3); plotCellData(G40_full, states40_full{end}.s(:,2), 'edgealpha', 0.2); view_oblique(); 
subplot(2,2,4); plotCellData(G40_full, sat40_VE_converted_3D_full{end}, 'edgealpha', 0.2); view_oblique();
set_fig_standardformat(fig, {'CO2 saturation, last timestep, 3D, z-resolution = 5', ...
                             'CO2 saturation, last timestep, 3D, z-resolution = 15', ...
                             'CO2 saturation, last timestep, 3D, z-resolution = 40', ...
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
