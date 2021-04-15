%% Visualizing in MRST
% MRST contains a suite of visualization routines which make it easy to
% create visualizations of grids and results. This tutorial show how to
% visualize grids, subsets of grids and details the different routines
% included in MRST. 

%% PLotting grids
% The |plotGrid| function is an essential part of MRST's visualization
% routines. It simply draws a grid to a figure with a reversed z axis.
G = cartGrid([10, 10, 3]);
G = computeGeometry(G);
clf;
plotGrid(G)
view(30,50)

%% MRST and patch
% All grid plotting routines are based on MATLAB's patch routine, which
% enables plotting of general polygons. Any keyword arguments will be
% passed on to patch, which makes it possible to alter many attributes.
% For instance, we can replot the grid with partially transparent edges and
% faces in another color:
clf;
plotGrid(G, 'EdgeAlpha', 0.1, 'FaceColor', 'blue')
view(30,50)

%% plotGrid and subsets
% plotGrid's second argument corresponds to a list of cells to be plotted.
% This can be either logical indices (a logical vector of length
% G.cells.num will plot Cell i if the ith element of the vector is true) or
% an explicit list of indices, i.e. the cell numbers to be plotted.
% To demonstrate this, we will plot all indices with equal values in a
% different color. Note that the plotting routines do not reset the figure
% between plots, making it easy to create compositions of different plots.
clf;
equal_index = mod(1:G.cells.num,2) == 0;
plotGrid(G,  equal_index, 'FaceColor', 'red')
plotGrid(G, ~equal_index, 'FaceColor', 'blue')
view(30,50)

%% Plotting cell properties
% To illustrate the visualization of cell properties, we will use
% routines from the |incomp| module to set up and solve a simple flow
% problem. This is only meant as an example for the visualization and will
% not be covered thourougly.
mrstModule add incomp
  
gravity off                                       % Disable gravity
rock = makeRock(G, 100*milli*darcy, 0.5);         % Uniform rock properties
fluid = initSimpleFluid( ...                      % A simple fluid
    'mu' , [   1,  10]*centi*poise     , ...      %  - viscosity
    'rho', [1014, 859]*kilogram/meter^3, ...      %  - density
    'n'  , [   2,   2]);                          %  - quadratic Corey relperms
W = verticalWell([], G, rock, 1, 1, [], ...       % Injector in cell (1,1)
            'Type', 'bhp' , 'Val', 1*barsa(), ... %   - bhp=1 bar
            'Radius', 0.1, 'InnerProduct', 'ip_tpf', ...
            'Comp_i', [0, 1]);
W = verticalWell(W, G, rock, 10, 10, [], ...      % Producer in cell (10,10)
            'Type', 'bhp' , 'Val', 0*barsa(), ... %   - bhp=0 bar
            'Radius', 0.1, 'InnerProduct', 'ip_tpf', ...
            'Comp_i', [1, 0]);
sol = initState(G, [], 0, [1, 0]);                % Initial: filled with phase 1

% Compute transmissibility and define function handles for solvers to avoid
% having to retype parameters that stay constant throught the simulation
T = computeTrans(G, rock);
psolve = @(state) incompTPFA(state, G, T, fluid, 'wells', W);
tsolve = @(state, dT) implicitTransport(state, G, dT, rock, fluid, 'wells', W);

%% Plot the pressure distribution
% Generally seeing only the grid is not that interesting. If we want to
% show actual values, we need to use plotCellData. Let us solve the
% reservoir pressure and plot it:
sol= psolve(sol);
clf;
plotCellData(G, sol.pressure)
colorbar, view(30,50)

%% Plot subset of data
% We can observe that the pressure has its largest values in the wells,
% which makes it seem that the pressure is almost homogeneous inside the
% domain. Let us then plot pressure values for the middle of the domain. We
% also plot an empty grid to see where we are actually plotting. In
% addition, we use the |plotWell| function to show wells.
[i,j,k] = ind2sub(G.cartDims, 1:G.cells.num);
clf;
plotGrid(G, 'FaceAlpha', 0, 'EdgeAlpha', .1)
plotCellData(G, sol.pressure, j == round(G.cartDims(2)/2))
plotWell(G, W); view(30,50)

%% plotFaces - plotting face values
% Some data is given at a per face basis instead of cells. For this purpose
% we have plotFaces and plotFaceData which are analogous to plotGrid and
% plotCellData respectively. To illustrate this function, we plot all faces
% having positive normal in the z-direction
clf;
plotFaces(G, find(G.faces.normals(:,3)>0));
view(30,50);

%% plotFaceData
% plotFaceData lets us plot face values, for instance when dealing with
% fluxes and faults. For an example, let us plot the faces colored by the z
% values of the face centroids.
clf;
plotFaceData(G, G.faces.centroids(:,3));
view(30,50);

%% Animated example
% We will now simulate a simple transport problem between the wells and
% create an animation that shows how the solution develops over time. To
% this end, we plot the empty grid and the wells, and in each time step
% only show the saturation in the cells where the saturation of the
% injected phase is above a threshold, here 0.05. Here, we solve and plot
% the solution simulaneously. We also store the time steps in a cell array
% so that they can be extracted after the simulation has finished

% Set up static parts of the plot
clf
plotGrid(G, 'FaceAlpha', 0, 'EdgeAlpha', .1)
plotWell(G, W);
view(30,50); pause(1)
hs = [];  % handle for saturation plot, empty initially

% Perform simulation
dT = 10*day;
solutions = cell(20,1);
for i = 1:20
    sol = tsolve(sol, dT);
    sol = psolve(sol);
    solutions{i} = sol;
    delete(hs)
    hs = plotCellData(G, sol.s(:,2), sol.s(:,2)>0.05);
    drawnow, pause(.5)
end

%% Alternatively, use plotGridVolumes for the same purpose
% The function |plotGridVolumes| is an alternative to |plotCellData|. While
% more computationally intensive, it plots interpolated surfaces
% corresponding to specific values of the saturation, giving a better
% indication of the saturation front and a less clear indication of the
% saturation mixture in sationary areas. These functions can be combined
% depending on the nature of the data set.

for i = 1:20
    clf;
    sol = solutions{i};
    plotGridVolumes(G, sol.s(:,2), 'basealpha', 2)
    plotGrid(G, 'FaceAlpha', 0, 'EdgeAlpha', .1);
    plotWell(G, W);
    view(30,50);
    pause(.5)
end

%%
displayEndOfDemoMessage(mfilename)

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
