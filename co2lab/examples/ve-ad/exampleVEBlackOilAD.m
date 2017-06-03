%% VE simulation in a standard black-oil solver
% In this example we show how to set up a standard format black-oil model that
% can be used to simulate a VE model. For the actual simulation, we use the
% fully-implicit solver in MRST, based on automatic differentiation.

mrstModule add co2lab ad-props deckformat ad-core ad-props ad-blackoil

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
% Create an input deck that can be used together with the fully-implicit black
% oil solver. Since the grid is constructed as part of setting up the input
% deck, we obtain it directly.
[deck, G] = sinusDeckAdiVE([nx ny nz], [Lx Ly H], nsteps, dt, ...
                         -.1*pi/180, depth, phi, perm, ...
                         (H*phi*Lx*Ly)*0.2*day/year, ipress);

% Alternatively, we could read deck from file and construct the grid
% deck = readEclipseDeck( ...
%    fullfile(mrstPath('co2lab'),'data','decks','sinusDeckAdiVE.DATA');
% G = initEclipseGrid(deck);

figure, plotGrid(G),view([0 -1 0]), box on
                      
                 
%% Initialize data structures
% First, we convert the input deck to SI units, which is the unit system used by
% MRST. Second, we initialize the rock parameters from the deck; the resulting
% data structure may have to be post-processed to remove inactive cells. Then we
% set up the fluid object.
deck  = convertDeckUnits(deck);
rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);
fluid = initDeckADIFluid(deck);

%systemOW  = initADISystem({'Oil', 'Water'}, G, rock, fluid);

%% Prepare simulation model, schedule and initial state
% Before we can run the simulation, we make sure that we have an initial
% hydrostatic pressure distribution.  We proceed by creating a simulation model
% object that the solver can work with, representing a two-phase oil/water
% system, where we let 'oil' represent the CO2 phase.  Finally, we convert the
% schedule from the input deck into MRST format.
x0 = initEclipseState(G, deck, initEclipseFluid(deck));
z  = G.cells.centroids(:,3);
x0.pressure = ipress*barsa +(z(:)-z(end))*norm(gravity)*deck.PROPS.DENSITY(2);
x0.sGmax = x0.s(:,2);

model = TwoPhaseOilWaterModel(G, rock, fluid);
schedule = convertDeckScheduleToMRST(model, deck);

%% Run the schedule setup in the file
% Before we can run the schedule, we make sure that we have an initial
% hydrostatic pressure distribution. Then we pick the schedule from the
% input deck and start the simulator.
[wellSols, states] = simulateScheduleAD(x0, model, schedule);


%% Plot results
%figure
Gt = topSurfaceGrid(G);
xc = Gt.cells.centroids(:,1);
zt = Gt.cells.z;
zb = zt + Gt.cells.H;
for nn=1:numel(states)
    clf
    state=states{nn};
    
    subplot(2,2,1),cla
    title('pressure')
    plotCellData(G,state.pressure/barsa,'EdgeColor','none');
    colorbar('horiz'), caxis([100 200])
    
    subplot(2,2,2),cla
    title('saturation')
    plotCellData(G,state.s(:,1),'EdgeColor','none');
    colorbar('horiz'), caxis([0 1])
    
    % plot as VE
    subplot(2,2,3),cla
    plot(xc,state.pressure/barsa); set(gca,'YLim',[100 200]);

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
