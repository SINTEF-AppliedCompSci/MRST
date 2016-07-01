%% Periodic Boundary Conditions for AD code
%
% We demonstrate the use of periodic boundary conditions. This is
% implemented as part of this module only using 'hacked' versions of the
% ad-blackoil models.
% 
% Note that there is a limitation in the implementation of the periodic
% boundary conditions, such each dimension must have at least two cells.
% 

mrstModule add ad-blackoil ad-core ad-props upscaling steady-state


%% Setup example
% 
% We set up a very simple incompressible model with periodic boundary
% conditions.

% Create grid
G = cartGrid([50, 2, 2], [100, 10, 10]*meter);
G = computeGeometry(G);

% Make the grid periodic
[Gp, bcp] = makePeriodicCartesianGrid(G);

% Homogenous rock properties
rock = struct('perm', 1000*ones(G.cells.num, 1)*milli*darcy, ...
              'poro', ones(G.cells.num, 1));

% Default fluid with unit values
fluid = initSimpleADIFluid();

% The fluid is incompressible. We have to tell the solver explicitly that
% this is the case to avoid a singular Jacobian matrix.
fluid.isIncomp = true;

% Set up model. We use a hacked version of the TwoPhaseOilWaterModel which
% allows for periodic boundary conditions.
model = TwoPhaseOilWaterModel_BCP(Gp, rock, fluid);

% Initial state with an oil slug in the middle of the domain
ijk = gridLogicalIndices(Gp);
inx = ijk{1} >= 25 & ijk{1} <= 40;
state0 = initResSol(Gp, 50*barsa, [1, 0]);
state0.s(inx,1) = 0;
state0.s(inx,2) = 1;
state0.wellSol = initWellSolAD([], model, state0);

% Set pressure drop over boundary
bcp.value(bcp.tags == 1) = -400*barsa;


%% Run simulation

solver = NonLinearSolver();
dT = 50*day;
n = 35;
states = cell(n+1, 1);
states{1} = state0;
for i = 1:n
    state = solver.solveTimestep(states{i}, dT, model, 'bcp', bcp);
    states{i+1} = state;
end


%% Plot animation of the solution
%
% Observe how the oil slug exits the boundary on the right and enters on
% the left as the grid is periodic.
% 

inx = ijk{2}==1 & ijk{3}==1;
x = G.cells.centroids(inx,1);
figure(1); clf;
for i=1:n
    cla;
    plot(x, states{i}.s(inx,2), 'r');
    axis([0 max(x) 0 1]);
    drawnow;
    pause(0.10); if i==1, pause(1); end
end



