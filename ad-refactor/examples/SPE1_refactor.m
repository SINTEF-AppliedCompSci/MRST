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

boModel = ThreePhaseBlackOilModel(G, rock, fluid, 'inputdata', deck, ...
                                                  'useCNVConvergence', false, ...
                                                  'nonlinearTolerance', 1e-10);
% boModel.useCNVConvergence = false;
% boModel.nonlinearTolerance = 1e-12;


ellipSolver = BackslashSolverAD();
linsolve = CPRSolverAD('ellipticSolver', ellipSolver);

%%
timer = tic();
[wellSols_base, states_base, report_base] = simulateScheduleAD(state, boModel, schedule, 'linearSolver', linsolve);
time_base = toc(timer);
%%
% res.dt = schedule.step.val(2:end);
clear res
res.dt = cumsum(schedule.step.val(2:end));
res.itcount = report_base.Iterations(2:end);
res.target = nan;
res.time = time_base;
rbase = res;
%%
results = [];
for i = [4, 8, 10, 15, 20]

    scheduleOnestep = schedule;
    scheduleOnestep.step.val = [1*day; sum(scheduleOnestep.step.val)];
    scheduleOnestep.step.control = [1; 1];

    % timestepper = SimpleTimeStepSelector('maxTimestep', 5*day, 'verbose', true);
    % timestepper = IterationCountTimeStepSelector('maxTimestep', 5*day, 'verbose', true);
    timestepper = GustafssonLikeStepSelector('targetIterationCount', i,...
                                             'minRelativeAdjustment', sqrt(eps),...
                                             'maxRelativeAdjustment', inf);

    nonlinear = NonLinearSolver('timeStepSelector', timestepper, 'verbose', true);

    timer = tic();
    [ws, s, reports] = simulateScheduleAD(state, boModel, scheduleOnestep, 'nonlinearSolver', nonlinear);
    time_dynamic = toc(timer);
    deltat = cellfun(@(x) x.LocalTime, reports.ControlstepReports{2}.StepReports);
    itcount = cellfun(@(x) x.Iterations, reports.ControlstepReports{2}.StepReports);
    
    
    good = cellfun(@(x) x.Converged, reports.ControlstepReports{2}.StepReports);
    deltat = deltat(good);
    itcount = itcount(good);
    
    res.dt = deltat;
    res.itcount = itcount;
    res.target = i;
    res.time = time_dynamic;
    results = [results; res];
end
%%
close all

pl = @(r, c) stairs(r.dt, r.itcount, 'linewidth', 2, 'color', c);
pl(rbase, 'k')

hold on
nr = numel(results);
c = lines(nr);
for i = 1:nr
    pl(results(i), c(i, :))
end
legend(['base'; arrayfun(@(x) num2str(x.target), results, 'unif', false)])
%%

