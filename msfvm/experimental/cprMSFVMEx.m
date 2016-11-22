%% Simulate a large example using parts of SPE10
% This example is a larger example demonstrating the solver on a medium
% size grid (66000 cells) with a relatively large amount of time steps
% (100). This example will take some time, especially if mldivide is used
% as the elliptic solver. Be vary that increasing the number of layers may
% let the simulations take a very long time.
mrstModule add ad-fi deckformat spe10

% Read and process file.
current_dir = fileparts(mfilename('fullpath'));
fn    = fullfile(current_dir, 'SPE10-S3.DATA.txt');

deck = readEclipseDeck(fn);

% The deck is given in field units, MRST uses metric.
deck = convertDeckUnits(deck);


% Create a special ADI fluid which can produce differentiated fluid
% properties.
fluid = initDeckADIFluid(deck);

% The case includes gravity
gravity on


% The initial state is provided as a binary file. The initial state
% contains a uniform mixture of water (.12) and oil (.88).
% load initialState;
%% Set up permeability, grid and wells
% We will simulate on the top 5 layers.
layers = 1:6;

[G, W, rock] = getSPE10setup(layers);


% SPE10 contains zero and extremely low porosities. For the purpose of this
% tutorial, we will mask away these values. An alternative would be to set
% these cells to inactive by using extractSubgrid and removing the
% corresponding cells.
low = 0.01;
rock.poro(rock.poro < low) = low;

%% Plot the permeability
clf;
plotCellData(G, log10(rock.perm(:,1)));

%% Set up solution structures.

% The initial reservoir is at 6000 psi and is fully oil saturated. The well
% solution gets its initial pressure from the bottom hole pressure values
% provided.
initSat = [0 1 0];
state0 = initResSol(G, 6000*psia, initSat);
state0.wellSol = initWellSolLocal(W, 6000*psia);

for i = 1:numel(W)
    state0.wellSol(i).pressure = W(i).val;
    % Set well sign
    if strcmpi(W(i).name(1), 'p')
        W(i).sign = -1;
    else
        W(i).sign = 1;
    end
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


%% MSFV setup
mrstModule add coarsegrid msfvm
cdim = G.cartDims./5;
cdim(3) = 2;

p = partitionUI(G, cdim);
CG = generateCoarseGrid(G, p);
DG = partitionUIdual(CG, cdim);
[~, DG] = createPermutationMatrix(sparse(G.cells.num, G.cells.num), DG, CG, 'Speedup', true);



%% Simulate 1000 days of production and save iteration count and time
% We provide the solver with time steps for roughly 1000 days of
% production. A few smaller steps are done to get better accuracy during
% the initial injection period. After this we do 10 day intervals to step
% rapidly through the schedule. While this converges at every time step,
% implicit solvers will still get improved accuracy by doing smaller time
% steps. Numerical diffusion can, for instance, be problematic when doing
% large time steps.

dt = [.001*day; .1*day*ones(5,1); 1*day*ones(10,1); 10*day*ones(100,1)];
nstep = 67;

states = cell(nstep,1);
statesfv = cell(nstep,1);
its = zeros(nstep,1);
itsfv = zeros(nstep,1);
timefv = zeros(nstep,1);
time = zeros(nstep,1);

msfvPressureSolve(NaN, DG, [], [])
% solve fv
systemfv = system;
systemfv.nonlinear.cprEllipticSolver = @(A, b) msfvPressureSolve(CG, DG, A, b, 'Speedup', true, 'UseCorrection', false);
systemfv.stepFunction = @(varargin) msfvCPRWrapper(system.stepFunction, varargin{:});

statefv = state0;
for t = 1 : nstep
    timer = tic();

    [statefv, it] = solvefiADI(statefv, dt(t), W, G, systemfv);

    statesfv{t} = statefv;
    itsfv(t) = it;
    timefv(t) = toc(timer);
end
%%
% Solve regular
state = state0;
for t = 1 : nstep
    timer = tic();

    [state, it] = solvefiADI(state, dt(t), W, G, system);

    states{t} = state;
    its(t) = it;
    time(t) = toc(timer);
end
%%
sum(time(1:50))/sum(timefv(1:50))

%%
figure(1);clf; plotToolbar(G, [states{1:50}])
figure(3);clf; plotToolbar(G, [statesfv{1:50}])

%%
clf;
subplot(1,2,1)
bar([sum(time(1:50)), sum(timefv(1:50))])

title('SPE10 subset solution time')
set(gca,'Xtick',1:2,'XTickLabel',{'CPR SPE10', 'MSFV-CPR SPE10'})
axis tight

subplot(1,2,2)
bar([630, 462])

title('SPE1 solution time')
set(gca,'Xtick',1:2,'XTickLabel',{'CPR SPE1', 'MSFV-CPR SPE1'})
axis tight

%%
clf
subplot(1,2,1)
plotCellData(G, states{50}.s(:,1), states{50}.s(:,1) > 0)
axis tight off
title('Standard CPR')

subplot(1,2,2)
plotCellData(G, statesfv{50}.s(:,1), statesfv{50}.s(:,1) > 0)
axis tight off
title('MSFV-CPR')
