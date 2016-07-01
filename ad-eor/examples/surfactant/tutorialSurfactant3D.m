%% Tutorial for a simple 3D case of a oil-water-surfactant system
% This example contains a $31\times31\times3$ fine grid containing two injectors
% in opposite corners and one producer in the middle of the domain. All wells
% are completed in the top layers of cells.
%
% The schedule being used contains first a period of injection with surfactant,
% followed by a water flooding phase without surfactant. Finally, the water rate
% is reduced for the final time steps.
%
% The data is read from deck (SURFACTANT3D.DATA)
%


%% We load the necessary modules
%

mrstModule add ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui

%% We load the input data and setup the grid, rock and fluid structures
%

current_dir = fileparts(mfilename('fullpath'));
fn = fullfile(current_dir, 'SURFACTANT3D.DATA');
gravity on

deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);

fluid = initDeckADIFluid(deck);
G = initEclipseGrid(deck);
G = computeGeometry(G);
rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);



%% Set up the initial state
% We want a layer of oil on top of the reservoir and water on the bottom.  To do
% this, we alter the initial state based on the logical height of each
% cell. Initially, the pressure is constant and the surfactant concentration is
% zero.

ijk = gridLogicalIndices(G);

state0 = initResSol(G, 300*barsa, [ .9, .1]);
state0.s(ijk{3} == 1, 2) = .9;
state0.s(ijk{3} == 2, 2) = .8;

% Enforce s_w + s_o = 1;
state0.s(:,1) = 1 - state0.s(:,2);

state0.c    = zeros(G.cells.num, 1);
state0.cmax = state0.c;
state0.ads = computeEffAds(state0.c, 0, fluid);
state0.adsmax = state0.ads;

clf
plotCellData(G, state0.s(:,2));
plotGrid(G, 'facec', 'none')
title('Oil concentration')
axis tight off
view(70, 30);
colorbar;


%% Set up the model
% The model object contains the grid, the fluid and rock properties and the
% modeling equations. See simulatorWorkFlowExample.
%
model = OilWaterSurfactantModel(G, rock, fluid, 'inputdata', deck);



%% Setup the schedule
% The wells and wells controls obtained from the input deck are parsed and a
% MRST schedule is set up.

schedule = convertDeckScheduleToMRST(model, deck);


%% Visualize some properties of the model we have setup
%

% We gathered visualizing command for this tutorial in the following script
example_name = '2D';
vizSurfactantModel;



%% Run the schedule
%
% We use the function simulateScheduleAD to run the simulation
% Options such as maximum non-linear iterations and tolerance can be set in
% the system struct.
%
% Here, we also send a resulthandler to save the output data, see ResultHandler.
%

resulthandler = ResultHandler('dataDirectory', pwd, 'dataFolder', 'cache', 'cleardir', true);
[wellSolsSurfactant, statesSurfactant] = simulateScheduleAD(state0, model, ...
                                                  schedule, 'OutputHandler', ...
                                                  resulthandler);
