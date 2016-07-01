%% 1D Tutorial For a Oil-Water-Surfactant system
% The input data is read from a deck using Eclipse format
% (SURFACTANT1D.DATA). The surfactant property (see file surfact.inc) are taken
% from SPE paper 145036.
%
% Surfactant is added to water in order to decrease the surface tension so
% that residual oil is mobilized. See more detail about the modeling equation
% in ad-eor/docs
%
% In this example, water and surfactant are injected at the left-hand side and
% oil is produced at the right-hand side at a given pressure.
%
% In a first period, only water is injected. Then, for a second period,
% surfactant is added to the water.

%% We load the necessary modules
%

mrstModule add ad-core ad-blackoil ad-eor ad-props deckformat mrst-gui


%% We load the input data and setup the grid, rock and fluid structures
% 

current_dir = fileparts(mfilename('fullpath'));
fn = fullfile(current_dir, 'SURFACTANT1D.DATA');
gravity off

deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);

fluid = initDeckADIFluid(deck);
G = initEclipseGrid(deck);
G = computeGeometry(G);
rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);


%% Set up the initial state
% Constant pressure, residual water saturation, no surfactant
%

nc = G.cells.num;
state0 = initResSol(G, 300*barsa, [ .2, .8]); % residual water saturation is 0.2
state0.c    = zeros(G.cells.num, 1);
state0.cmax = state0.c;
state0.ads = computeEffAds(state0.c, 0, fluid); % adsorption
state0.adsmax = state0.ads;

%% Set up the model
% 
% The model object contains the grid, the fluid and rock properties and the
% modeling equations. See simulatorWorkFlowExample.
%


model = OilWaterSurfactantModel(G, rock, fluid, ...
                                                  'inputdata', deck, ...
                                                  'extraStateOutput', true);

%% Convert the deck schedule into a MRST schedule by parsing the wells
%

schedule = convertDeckScheduleToMRST(model, deck);


%% Visualize some properties of the model we have setup
%

% We gathered visualizing command for this tutorial in the following script
example_name = '1D';
vizSurfactantModel;



%% Run the schedule
%
% We use the function simulateScheduleAD to run the simulation
% Options such as maximum non-linear iterations and tolerance can be set in
% the system struct.

[wellSolsSurfactant, statesSurfactant] = simulateScheduleAD(state0, model, ...
                                                  schedule);

figure()
plotToolbar(G, statesSurfactant, 'startplayback', true, 'plot1d', true)
