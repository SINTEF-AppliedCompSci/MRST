%% Read the problem from a deckfile
% The problem is defined in 'INPUT_NUMGRAD.DATA' which is a simple
% $10\times1\times10$ Cartesian grid with uniform permeability. We read the
% deck and create the grid, rock and fluid structures from the resulting
% output. This requires the deckformat module.
mrstModule add deckformat ad-fi ad-refactor

current_dir = fileparts(mfilename('fullpath'));
fn    = fullfile(current_dir, 'simple10x1x10.data');
deck = readEclipseDeck(fn);

% Convert to MRST units (SI)
deck = convertDeckUnits(deck);

% Create grid
G = initEclipseGrid(deck);

% Set up the rock structure
rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

% Create fluid
fluid = initDeckADIFluid(deck);

% Get schedule
schedule = deck.SCHEDULE;

% Enable this to get convergence reports when solving schedules
verbose = false;


G = computeGeometry(G);
T = computeTrans(G, rock);

gravity on

state = initResSol(G, deck.PROPS.PVCDO(1), [.15, .85]);

% scalFacs.pressure = 100*barsa;
% scalFacs.rate     = 100/day;

scalFacs.pressure = 1;
scalFacs.rate     = 1;
%% Run the whole schedule
% This is done to get values for the wells for all timesteps. Since the
% case is fairly small,
timer = tic;
system = initADISystem({'Oil', 'Water'}, G, rock, fluid);
[wellSols, states] = runScheduleADI(state, G, rock, system, schedule);
t_forward = toc(timer);

%% Create objective functions

objective_adjoint = @(tstep)NPVOW(G, wellSols, schedule, 'ComputePartials', true, 'tStep', tstep);
objective_numerical = @(wellSols)NPVOW(G, wellSols, schedule);

%% Compute derivatives using the adjoint formulation

timer = tic;
getEquations = @eqsfiOWExplicitWells;
adjointGradient = runAdjointADI(G, rock, fluid, schedule, objective_adjoint, system,  'Verbose', verbose, 'ForwardStates', states);
t_adjoint = toc(timer);

%% Find gradients numerically

timer = tic;
numericalGradient = computeNumGrad(state, G, rock, system, schedule, objective_numerical, 'scaling', scalFacs, 'Verbose', verbose);
t_gradient = toc(timer);

%%

model = TwoPhaseOilWaterModel(G, rock, fluid, 'inputdata', deck);
schedulenew = convertDeckScheduleToMRST(G, rock, deck);

% Ensure that solver produces multiple substeps
% ms = schedulenew.step.val(end)/2;
ms = inf;
stepsel = SimpleTimeStepSelector('maxTimestep', ms);       
nonlinear = NonLinearSolver('timeStepSelector', stepsel);

usemini = true;
[wellSols, states, report] = simulateScheduleAD(state, model, schedulenew,...
                'OutputMinisteps', usemini, 'NonLinearSolver', nonlinear);
            

%%
if usemini
    schedulemini = convertReportToSchedule(report, schedulenew);
else
    schedulemini = schedulenew;
end

obj = @(tstep)NPVOW(G, wellSols, schedulemini, 'ComputePartials', true, 'tStep', tstep);

grad = simulateAdjointAD(state, states, model, schedulemini, obj);
% vertcat(adjointGradient{:})
% grad(end-5:end)
%% Plot the gradients


wellNames = {wellSols{1}.name};

ga = cell2mat(adjointGradient);
ga_n = full(cell2mat(grad'));
gn = cell2mat(numericalGradient);
clf;
subplot(1,2,1)
plot(ga(1,:),'-ob'), hold on
plot(gn(1,:),'-xr')
plot(ga_n(1,:),'-*g')

title(['Well 1 (', wellNames{1}, ')'])
xlabel('Control #')

subplot(1,2,2)
plot(ga(2,:),'-ob'), hold on
plot(gn(2,:),'-xr'), hold on
plot(ga_n(2,:),'-*g')

title(['Well 2 (', wellNames{2}, ')'])
xlabel('Control #')
legend({'Adjoint-old', 'Adjoint-new', 'Numerical'})

