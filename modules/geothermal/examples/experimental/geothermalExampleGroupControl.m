%% Add necessary modules
mrstModule add ruden-geothermal                              % Project
mrstModule add ad-core ad-props geothermal compositional     % Physics
mrstModule add dfm
mrstModule add test-suite                                    % Test suite
mrstModule add coarsegrid libgeometry upr                    % Gridding
mrstModule add linearsolvers
mrstModule add mrst-gui                                      % Visualization

mrstVerbose on % Set to `off` to supress command window output

%%
theta = @(n) linspace(0, 2*pi-2*pi/n, n)';
x = @(r,n) r.*[cos(theta(n)), sin(theta(n))];

bdr = x(100*meter, 50);
xwi = x(10*meter, 4);
xwp = x(50*meter, 8);

rng(20220501)
G = pebiGrid2D(100/20, [100,100], ...
    'cellConstraints', num2cell([xwi; xwp], 2)', ...
    'polyBdr', bdr, ...
    'CCrefinement', true, ...
    'CCFactor', 0.2);

G = computeGeometry(G);

%%
subGroups = false;

rock = makeRock(G, 100*milli*darcy, 0.1);
rock = addThermalRockProps(rock);
fluid = initSimpleADIFluid('phases', 'W', 'n', 1, 'rho', 1, 'mu', 1);
fluid = addThermalFluidProps(fluid, 'useEOS', true);

model = GeothermalModel(G, rock, fluid);

% mrstModule add ad-blackoil
% fluid = initSimpleADIFluid('phases', 'W', 'n', 1, 'rho', 1000, 'mu', 1*centi*poise,  'c', 1e-3/barsa);
% model = GenericBlackOilModel(G, rock, fluid, 'oil', false, 'gas', false);

W = [];
nw = nnz(G.cells.tag);
group = cell(nw,1);

ilim = struct('bhp', 500*barsa);
time = 1*year;
rate = 2*sum(poreVolume(G, rock))/time;
nwi = size(xwi,1);
for i = 1:nwi
    cells = findEnclosingCell(G, xwi(i,:));
    W = addWell(W, G, rock, cells, 'type', 'rate', 'val', rate/nwi, 'lims', ilim, 'Name', ['I-', num2str(i)]);
    if subGroups
        group{i} = ['Injectors-', num2str(2-rem(i,2))];
    else
        group{i} = 'Injectors';
    end
end

nwp = size(xwp,1);
plim = struct('bhp', 1*atm);
for i = 1:nwp
    cells = findEnclosingCell(G, xwp(i,:));
    W = addWell(W, G, rock, cells, 'type', 'rate', 'val', -rate/nwp, 'lims', plim, 'Name', ['P-', num2str(i)]);
    group{i+nwi} = 'Producers';
end
[W.group] = deal(group{:});

W = addThermalWellProps(W, G, rock, fluid, 'T', convertFromCelcius(80));

[W.WI] = deal(W(1).WI);

wnames = {W.name};
groups = [];
groups = addFacilityGroup(groups,  wnames(1:4), ...
    'name', 'Injectors', ...
    'type', 'rate', ...
    'val', rate ...
);

groups = addFacilityGroup(groups, wnames(5:end), ...
        'name', 'Producers', ...
        'type', 'rate', ...
        'val', nan ...
);
groups(2).val = 'Injectors';

[groups.T] = deal(convertFromCelcius(80));
       
dt = rampupTimesteps(time, 10*day);

bc = addBC([], boundaryFaces(G), 'pressure', 50*barsa);
bc = addThermalBCProps(bc, 'T', convertFromCelcius(20));

schedule = simpleSchedule(dt, 'W', W, 'bc', bc);
schedule.control(1).groups = groups;

state0 = initResSol(G, 50*barsa, 1);
state0.T = repmat(convertFromCelcius(20), G.cells.num, 1);

%%
lsol = selectLinearSolverAD(model);
nlsol = NonLinearSolver('LinearSolver', lsol);
problem = packSimulationProblem(state0, model, schedule, 'group-test', 'NonLinearSolver', nlsol);
simulatePackedProblem(problem, 'restartStep', 1);

%%
close all
[wellSols, states, reports] = getPackedSimulatorOutput(problem);
[groupSols, wellSols] = computeWellGroupSols(wellSols, schedule.control(1).W);
% wellSols = computeWell
plotToolbar(model.G, states); axis equal tight; colormap(hot);
plotWellSols(wellSols);