%% Read the problem from a deckfile
% The problem is defined in 'INPUT_NUMGRAD.DATA' which is a simple
% $10\times1\times10$ Cartesian grid with uniform permeability. We read the
% deck and create the grid, rock and fluid structures from the resulting
% output. This requires the deckformat module.
mrstModule add deckformat ad-fi ad-refactor

current_dir = fileparts(mfilename('fullpath'));
deck = readEclipseDeck(fullfile(current_dir, 'simple10x1x10.data'));

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

scalFacs.pressure = 100*barsa;
scalFacs.rate     = 100/day;

% scalFacs.pressure = 1;
% scalFacs.rate     = 1;
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

% Create result handler that writes to disk
clear output
output = ResultHandler('writeToDisk', true, 'storeInMemory', false);

[wellSols, states, report] = simulateScheduleAD(state, model, schedulenew,...
                'OutputMinisteps', usemini, 'NonLinearSolver', nonlinear, ...
                'OutputHandler', output);
            

%%
if usemini
    schedulemini = convertReportToSchedule(report, schedulenew);
else
    schedulemini = schedulenew;
end

obj = @(tstep)NPVOW(G, wellSols, schedulemini, 'ComputePartials', true, 'tStep', tstep);

% d = states;
d = output;

adjointGradientClass = computeGradientAdjointAD(state, d, model, schedulemini, obj, 'scaling', scalFacs);
% vertcat(adjointGradient{:})
% grad(end-5:end)
%%
numericalGradientClass = computeGradientPerturbationAD(state, model, schedulenew, objective_numerical, 'scaling', scalFacs);
%% Plot the gradients

wellNames = {wellSols{1}.name};

adjointNew = cell2mat(adjointGradientClass);
adjointOld = cell2mat(adjointGradient);

numericalNew = cell2mat(numericalGradientClass);
numericalOld = cell2mat(numericalGradient);


% numericalNew = full(cell2mat(grad));
% gn = cell2mat(numericalGradient);
clf;

for i = 1:2
    subplot(1,2,i)
    hold on
    plot(adjointNew(i,:),'-.b'), 
    plot(adjointOld(i,:),'-*k');
    plot(numericalNew(i,:),'-og');
    plot(numericalOld(i,:),'-+r');
    
    title(['Well ', num2str(i), ' (', wellNames{i}, ')'])
    xlabel('Control #')

    xlabel('Control #')
    legend({'Adjoint-new', 'Adjoint-old', 'Numerical-new', 'Numerical-old'})

end
