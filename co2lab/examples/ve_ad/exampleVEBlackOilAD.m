%% VE simulation in a standard black-oil solver
% In this example we show how to set up a standard format black-oil model that
% can be used to simulate a VE model. For the actual simulation, we use the
% fully-implicit solver in MRST, based on automatic differentiation.

mrstModule add co2lab ad-props deckformat ad-core ad-blackoil

%% Parameters for the simulation
gravity on
[nx,ny,nz] = deal(40, 1, 1);     % Cells in Cartsian grid
[Lx,Ly,H]  = deal(2000,1000,15); % Physical dimensions of reservoir
total_time = 5*year;             % Total simulation time
nsteps     = 10;                 % Number of time steps in simulation
dt         = total_time/nsteps;  % Time step length
perm       = 100;                % Permeability in milli darcies
phi        = 0.1;                % Porosity
depth      = 1000;               % Initial depth
ipress     = 200;                % Initial pressure

%% Create input deck and construct grid
% Create an input deck that can be used together with the fully-implicit
% black oil solver. Since the grid is constructed as part of setting up the
% input deck, we obtain it directly.
deck = sinusDeckAdiVE([nx, ny, nz], [Lx, Ly, H], nsteps, dt, ...
                      -.1*pi/180, depth, phi, perm, ...
                      (H*phi*Lx*Ly)*0.2*day/year, ipress);

deck = convertDeckUnits(deck);

%% Initialize data structures
[x0, model, schedule] = initEclipseProblemAD(deck);

%% Show simulation grid
% This is stored as the |G| data member of the simulation |model|.
figure, plotGrid(model.G), view([0, -1, 0]), box on

%% Prepare simulation model, schedule and initial state
% Before we can run the simulation, we make sure that we have an initial
% hydrostatic pressure distribution.  We proceed by creating a simulation
% model object that the solver can work with, representing a two-phase
% oil/water system, where we let 'oil' represent the CO2 phase.  Finally,
% we convert the schedule from the input deck into MRST format.
z  = model.G.cells.centroids(:,3);
x0.pressure = ipress*barsa +(z(:)-z(end))*norm(gravity)*deck.PROPS.DENSITY(2);
x0.sGmax = x0.s(:,2);

%% Run the schedule setup in the file
% Before we can run the schedule, we make sure that we have an initial
% hydrostatic pressure distribution. Then we pick the schedule from the
% input deck and start the simulator.
[wellSols, states] = simulateScheduleAD(x0, model, schedule);

%% Plot results
%figure
Gt = topSurfaceGrid(model.G);
xc = Gt.cells.centroids(:,1);
zt = Gt.cells.z;
zb = zt + Gt.cells.H;
for nn=1:numel(states)
    clf
    state=states{nn};
    
    subplot(2,2,1),cla
    title('pressure')
    plotCellData(model.G,convertTo(state.pressure,barsa),'EdgeColor','none');
    colorbar('horiz'), caxis([100 200])
    
    subplot(2,2,2),cla
    title('saturation')
    plotCellData(model.G,state.s(:,1),'EdgeColor','none');
    colorbar('horiz'), caxis([0 1])
    
    % plot as VE
    subplot(2,2,3),cla
    plot(xc,convertTo(state.pressure,barsa)); set(gca,'YLim',[100 200]);

    subplot(2,2,4),cla,hold on
    patch(xc([1 1:end end]), [zt(end)-10; zt; zt(end)-10],.7*[1 1 1]);
    patch(xc([1 1:end end]), [zb(end)+10; zb; zb(end)+10],.7*[1 1 1]);
    patch(xc([1:end end:-1:1]), ...
      [zt + Gt.cells.H.*state.s(:,2); zt(end:-1:1)], getVEColors('plume'))
    patch(xc([1:end end:-1:1]), ...
      [zt + Gt.cells.H.*state.s(:,2); zb(end:-1:1)], getVEColors('brine'))
    set(gca,'YDir','reverse'), axis tight
    
    drawnow;
    pause(0.01)
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
