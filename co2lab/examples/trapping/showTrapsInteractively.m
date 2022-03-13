%% Interactive spill point analysis of top surface grids
% This example demonstrates the use of the interactive viewer for viewing
% and finding structural traps for CO2 migration. Three different examples
% are presented: 
% 
% # a small conceptual model
% # a model of the Johansen formation from the CO2 Storage Atlas
% # a large model from the IGEMS project
%
% The following functionality is present in the interactive viewer:
%
% - Left clicking on a trap will colorize the trap and all other traps that
%   are downstream to that trap along with the path between them. Left
%   clicking shows an estimate of where and how any slowly injected CO2
%   will migrate from the click position. When left-clicking, two
%   additional plots are produced:
%
%    * The upper figure is a pie chart showing the approximate volumes
%      of both the clicked trap (direct volume) and the volume of traps
%      along the migration path (migration volume) as well as the volume of
%      the traps which are not on the migration path. By exploring the grid
%      interactively, one can try to find the best possible injection site
%      for a top surface grid with regards to long term migration.
%
%    * The second figure shows a logarithmic bar chart of the trap volumes
%      along with their position along the migration path. This is useful
%      for several reasons; one being that the quality of the dataset may
%      impact the correctness of far away traps, another being that far
%      away traps will likely meet a slower CO2 front, making the model of
%      infinitesimally slow migration better.
%
% - Right or ctrl-clicking on a trap will open a new plot showing the trap
%   in detail along with the region around it. The largest possible
%   structural trapping volume for the region will be indicated in red and
%   blue.
%
% - Middle mouse or shift-clicking on a trap will show all traps that are
%   upstream to the current trap much in the same manner as left mouse
%   click. This can be used to find possible injection sites which will
%   migrate over time to a large reservoir volume.
%
% - The toolbar in the figure window has several functions:
%
%    * Toggle display of unselected traps.
%
%    * Toggle display of spill regions. Spill regions determine which trap
%      gas injected in an area will migrate to. Regions which spill out of
%      the grid are marked in blue.
%
%    * Toggle lighting. Lighting can be computationally intensive in
%      MATLAB, but makes it easier to see local changes in topology.
% 
%    * Toggle contour lines. Only available when plotting atlas grids.
%
%    * Reset view. The selection algorithm works best from a top down
%      view.
%
%    * Setup VE simulation. The currently clicked cell will become an
%      injector in a vertically averaged CO2 migration simulation. The user
%      can setup injection rates and timesteps through a simple user
%      interface and stop simulations by closing the visualization window.

mrstModule add co2lab coarsegrid
% Select trapping algorithm
useCell = true;

switch useCell
    case true
        method = 'cell';
    otherwise
        method = 'node';
end

% Check if the script is run interactively, or as a standalone script
isScript = numel(dbstack) == 1;

%% Explore a synthetic grid
% We first create a simple synthetic grid to demonstrate the interactive
% trapping. This grid is created by producing a simple Cartesian geometry
% and perturbing the z coordinates by a periodic function based on x and y
% to get several local traps. Additionally, to create a hierachy of traps,
% we slant the grid.
%
% We plot the top surface grid in the interactive viewer. Initially the
% possible structural traps are shown, colorized by the volume of each
% structural trap on a white grid.
 
L = 1000; L_p = L/5; H = 10;
G = cartGrid([101 101 1],[L L H]);
G.nodes.coords(:,3) =   100 + G.nodes.coords(:,3) + ...
                              G.nodes.coords(:,1)*0.01...
                    -2*sin(pi*G.nodes.coords(:,1)/L_p).*...
                       sin(pi*G.nodes.coords(:,2)/L_p);
G = computeGeometry(G);

% Create top surface grid
Gt = topSurfaceGrid(G);

% Show interactive plot
h = interactiveTrapping(Gt, 'method', method, 'light', true);
view(180, 50)

disp('Showing synthetic dataset... Close to continue')
if isScript
    waitfor(h)
end

%% The Johansen formation
% By passing in the name of a CO2 Storage Atlas grid, the function
% automatically generates the grid for us. We select the Johansen formation
% and apply some coarsening to make the grid smaller. Note that full
% resolution should be used if possible to get the best accuracy.

h = interactiveTrapping('Johansenfm', 'coarsening', 2, ...
   'method', method, 'light', true, 'contour', true);
view(-70, 80)

disp('Showing CO2 Atlas dataset... Close to continue')
if isScript
    waitfor(h)
end

%% Grid from the IGEMS 
% The purpose of the IGEMS project was to study the impact that top-surface
% morphology has on storage capacity and the migration process. Alternative
% top-surface morphologies are created stochastically by combining
% different stratigraphic scenarios with different structural scenarios.
% For more information about the data set, see
% http://files.nr.no/igems/data.pdf.

G = readIGEMSIRAP('OSS', 1, 'coarse', [2 2]);
Gt = topSurfaceGrid(G);

% Activate trapping
h = interactiveTrapping(Gt, 'method', method);
view(180, 50)

disp('Showing IGEMS dataset... Close to continue')
if isScript
    waitfor(h)
end

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
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
