%% Simulate a large example using parts of SPE10
% This example is a larger example demonstrating the solver on a medium
% size grid (66,000 cells) with a relatively large amount of time steps
% (100). This example will take some time, especially if MLDIVIDE is used
% as the elliptic solver. Be wary that increasing the number of layers may
% let the simulations take a very long time.
%
% Furthermore, note that the 'solvefiADI' solver and especially the CPR
% preconditioner, uses direct indexing into sparse matrices.  In MATLABs
% prior to release R2011a, this is relatively inefficient.  Starting from
% 2011a however, the bottleneck of direct indexing has been largely
% removed.
mrstModule add ad-fi deckformat spe10 ad-refactor

% Read and process file.
current_dir = fileparts(mfilename('fullpath'));
fn    = fullfile(current_dir, 'SPE10-S3.DATA.txt');

deck = readEclipseDeck(fn);

% The deck is given in field units, MRST uses metric.
deck = convertDeckUnits(deck);


% Create a special ADI fluid which can produce differentiated fluid
% properties.
fluid = initDeckADIFluid(deck);

% Oil rel-perm from 2p OW system.
% Needed by equation implementation function 'eqsfiOWExplictWells'.
fluid.krO = fluid.krOW;

% The case includes gravity
gravity on


% The initial state is provided as a binary file. The initial state
% contains a uniform mixture of water (.12) and oil (.88).
% load initialState;
%% Set up permeability, grid and wells
% We will simulate on the top 5 layers.
layers = 1:15;

[G, W, rock] = SPE10_setup(layers);


% SPE10 contains zero and extremely low porosities. For the purpose of this
% tutorial, we will mask away these values. An alternative would be to set
% these cells to inactive by using extractSubgrid and removing the
% corresponding cells.
low = 1e-4;
rock.poro(rock.poro < low) = low;

%% Plot the permeability
clf;
plotCellData(G, log10(rock.perm(:,1)));

%% Set up solution structures.

% The initial reservoir is at 6000 psi and is fully oil saturated. The well
% solution gets its initial pressure from the bottom hole pressure values
% provided.

for i = 1:numel(W)
    if strcmpi(W(i).name(1), 'p')
        W(i).sign = -1;
    else
        W(i).sign = 1;
    end
end

initSat = [0 1];
state0 = initResSol(G, 6000*psia, initSat);
state0.wellSol = initWellSolLocal(W, state0);

for i = 1:numel(W)
    state0.wellSol(i).bhp = W(i).val;
end

% Set up a Water / Oil system using CPR preconditioner. Alternatively we
% could have used a specialized elliptic solver using the option 'cprEllipticSolver'
% to exploit the nature of the variables involving pressure.

system = initADISystem({'Water', 'Oil'}, G, rock, fluid, 'cpr', true, 'cprRelTol', 2e-2);
% If an alternative solver for the pressure subproblem was installed, it
% could be added using
%
%    system.nonlinear.cprEllipticSolver = @(A,b) solver(A,b)
%
% This can greatly speed up the solution process, as the default option
% uses MATLABs direct solver @mldivide which is very expensive for a
% preconditioner.

mrstModule add agmg
system.nonlinear.cprEllipticSolver = @(A,b) agmg(A,b);
%% Simulate 1000 days of production and save iteration count and time
% We provide the solver with time steps for roughly 1000 days of
% production. A few smaller steps are done to get better accuracy during
% the initial injection period. After this we do 10 day intervals to step
% rapidly through the schedule. While this converges at every time step,
% implicit solvers will still get improved accuracy by doing smaller time
% steps. Numerical diffusion can, for instance, be problematic when doing
% large time steps.

dt = [        0.001*day;
      repmat( 0.1  *day, [  5, 1]);
      repmat( 1    *day, [ 10, 1]);
      repmat(10    *day, [100, 1])];
  
% dt = dt(1:5);
  
nstep = numel(dt);

states = cell(nstep,1);
its = zeros(nstep,1);
time = zeros(nstep,1);

state = state0;
%%
for t = 1 : nstep
    fprintf('Step %d/%d: %.2f -> %.2f [days]\n\n', t, nstep, ...
            convertTo(sum(dt(1:t-1)), day), ...
            convertTo(sum(dt(1:t  )), day));

    timer = tic();

    [state, it] = solvefiADI(state, dt(t), W, G, system);

    states{t} = state;
    its(t) = it;
    time(t) = toc(timer);
    
    time(t)

    fprintf('\n\n');
end

%%
% schedule.control.W = W;
schedule = struct();
schedule.control.W = W;
schedule.step.val = dt;
schedule.step.control = ones(size(dt));

clear boModel
clear nonlinear

boModel = twoPhaseOilWaterModel(G, rock, fluid, 'deck', deck);
%%

relTol = 5e-3;

mrstModule add coarsegrid mrst-experimental

cdims = round(G.cartDims./[10 10 5]);
% p = partitionUI(G, cdims);
p = partitionMETIS(G, computeTrans(G, rock), prod(cdims));
CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);
CG = storeInteractionRegion(CG);

multiscaleSolver = multiscaleVolumeSolverAD(CG);

owModel = twoPhaseOilWaterModel(G, rock, fluid, 'deck', deck);
linsolve = CPRSolverAD('ellipticSolver', multiscaleSolver, 'relativeTolerance', relTol);

%%
% linsolve = CPRSolverAD();
% nonlinear = nonlinearSolver();
% [state, status] = nonlinear.solveTimestep(state, 1*day, boModel)

timer = tic();
[wellSols, states] = runScheduleRefactor(state, boModel, schedule, 'verbose', true, 'linearSolver', linsolve);
time_ms = toc(timer);

%%
amgsolver = AGMGSolverAD();

basicCPR = CPRSolverAD('ellipticSolver', amgsolver, 'relativeTolerance', relTol);
timer = tic();
[wellSols, states] = runScheduleRefactor(state, boModel, schedule, 'verbose', true, 'linearSolver', basicCPR);
time_cpr = toc(timer);

%%
bar([time_cpr, time_ms])