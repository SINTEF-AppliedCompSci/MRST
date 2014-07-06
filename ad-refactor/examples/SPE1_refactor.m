clear
mrstModule add ad-fi deckformat mrst-gui ad-refactor

% Read and process file.
current_dir = fileparts(mfilename('fullpath'));
fn    = fullfile(current_dir, 'odeh_adi.data');

deck = readEclipseDeck(fn);

% The deck is given in field units, MRST uses metric.
deck = convertDeckUnits(deck);

G = initEclipseGrid(deck);
G = computeGeometry(G);

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

% Create a special ADI fluid which can produce differentiated fluid
% properties.
fluid = initDeckADIFluid(deck);

% The case includes gravity
gravity on


% The initial state is a pressure field that is constant in each layer, a
% uniform mixture of water (Sw=0.12) and oil (So=0.88) with no initial free
% gas (Sg=0.0) and a constant dissolved gas/oil ratio ("Rs") throughout the
% model.
%
% The pressure and Rs values are derived through external means.
clear prod
[k, k] = ind2sub([prod(G.cartDims(1:2)), G.cartDims(3)], ...
                  G.cells.indexMap);  %#ok

p0    = [ 329.7832774859256 ; ...  % Top layer
          330.2313357125603 ; ...  % Middle layer
          330.9483500720813 ];     % Bottom layer

p0    = convertFrom(p0(k), barsa);
s0    = repmat([ 0.12, 0.88, 0.0 ], [G.cells.num, 1]);
rs0   = repmat( 226.1966570852417 , [G.cells.num, 1]);
rv0   = 0; % dry gas

state = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);   clear k p0 s0 rs0

schedule = convertDeckScheduleToMRST(G, rock, deck);

clear boModel
clear nonlinear

boModel = ThreePhaseBlackOilModel(G, rock, fluid, 'inputdata', deck);
%%


ellipSolver = BackslashSolverAD();
linsolve = CPRSolverAD('ellipticSolver', ellipSolver);


timer = tic();
[wellSols, states] = simulateScheduleAD(state, boModel, schedule, 'linearSolver', linsolve);
time_ms = toc(timer);

%%


scheduleOnestep = schedule;
scheduleOnestep.step.val = [1*day; sum(scheduleOnestep.step.val)/10];
scheduleOnestep.step.control = [1; 1];

% timestepper = SimpleTimeStepSelector('maxTimestep', 5*day, 'verbose', true);
timestepper = IterationCountTimeStepSelector('maxTimestep', 5*day, 'verbose', true);

nonlinear = NonLinearSolver('timeStepSelector', timestepper, 'verbose', true);

[ws, s, reports] = simulateScheduleAD(state, boModel, scheduleOnestep, 'nonlinearSolver', nonlinear);


