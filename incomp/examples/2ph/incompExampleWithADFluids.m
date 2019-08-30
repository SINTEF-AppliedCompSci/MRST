%% Interopability between incomp solvers and ad-props fluids
% We can use the fluid models that are used in the ad-blackoil and ad-core
% modules with the various incompressible solvers. This example
% demonstrates this on a simple synthetic test problem.
mrstModule add ad-blackoil ad-core ad-props mrst-gui
%% Set up model
% We set up a 1D Buckley-Leverett problem to have something to test.

% Build grid
G = cartGrid([50, 1, 1], [100, 10, 10]*meter);
G = computeGeometry(G);

% Homogenous rock properties
rock = makeRock(G, 5*darcy, .3);

% Default oil-water fluid with unit values
fluid = initSimpleADIFluid('phases', 'WO', 'n', [2 2]);
% Comment out the next line to get a very strong capillary pressure term
% (not suitable for explicitTransport).
% fluid.pcOW = @(s, varargin) -50*s*barsa;

% Set up model and initial state.
model = TwoPhaseOilWaterModel(G, rock, fluid);
state0 = initResSol(G, 50*barsa, [0, 1]);

% Set up drive mechanism: constant rate at x=0, constant pressure at x=L
pv = poreVolume(G, rock);
time = 10000*day;
injRate = sum(pv)/time;
W = [];
W = addWell(W, G, rock, 1, 'type', 'rate', 'val', injRate, 'comp_i', [1, 0]);
W = addWell(W, G, rock, G.cells.num, 'type', 'bhp', 'val', 0, 'comp_i', [1, 0]);

%% Solve with ad-blackoil fully-implicit solver
% We first solve the problem with the fully-implicit solvers.
n  = 100;
dT = time/n;
schedule = simpleSchedule(repmat(dT,1,n), 'W', W);
[~, states, report] = simulateScheduleAD(state0, model, schedule);
%% Solve using incomp module
% We can now solve the same problem with the incomp module. The different
% incompressible solvers will work directly with fluids from ad-props,
% under a few limiting assumptions:
%   - The two-phase system is assumed to be a water-oil system where the
%   gas phase will be ignored if present in fluid. Consequently, only pcOW,
%   krO, krOW and krW will be used and pcOG and krG is ignored.
%   - It is not recommended that the fluid model is compressible as the
%   incomp solvers do not account for these nonlinear effects. We still
%   evaluate the b-factors at the given pressure in state (or 1 atmosphere
%   if no pressure is present in state) to obtain the densities.
%   - The explicit transport solver does not handle spatial changes in
%   phase viscosity when estimating CFL.

mrstModule add incomp
% Choices for pressure
pressure_discretization = 'tpfa';
% pressure_discretization = 'mpfa';
% pressure_discretization = 'mimetic';

% Choices for transport
transport_discretization = 'implicit';
% transport_discretization = 'explicit';
state = state0;
states_incomp = cell(n, 1);

switch pressure_discretization
    case 'tpfa'
        T = computeTrans(G, rock);
        psolve = @(state) incompTPFA(state, G, T, fluid, 'W', W);
    case 'mpfa'
        mrstModule add mpfa
        T = computeMultiPointTrans(G, rock);
        psolve = @(state) incompMPFA(state, G, T, fluid, 'W', W);
    case 'mimetic'
        IP = computeMimeticIP(G, rock);
        psolve = @(state) incompMimetic(state, G, IP, fluid, 'W', W);
    otherwise
        error('Unknown pressure solver %s', p_solve);
end

switch transport_discretization
    case 'explicit'
        tsolve = @(state, dT) explicitTransport(state, G, dT, rock, fluid, 'W', W);
        substeps = 10;
    case 'implicit'
        tsolve = @(state, dT) implicitTransport(state, G, dT, rock, fluid, 'W', W);
        substeps = 1;
    otherwise
        error('Unknown transport solver %s', transport_discretization);
end
% Solve with a loop
for i = 1:n
    state = psolve(state);
    for j = 1:substeps
        state = tsolve(state, dT/substeps);
    end
    states_incomp{i} = state;
end
%% Compare the results
% We compare the results of the two solvers. The results will agree up to
% the choice of time integration (our AD solver is fully-implicit, while
% the incomp module uses either IMPES (implicit pressure, explicit
% transport) or IMPSAT (implicit pressure, implicit saturation
% sequentially).
solvename = ['Incomp (', pressure_discretization, '-', transport_discretization, ')'];
x = G.cells.centroids(:, 1);
for i = 1:n
    figure(1); clf; hold on
    plot(x, states{i}.s(:, 1));
    plot(x, states_incomp{i}.s(:, 1), 'o');
    ylim([0, 1]);
    legend('Fully-implicit AD', solvename)
end
