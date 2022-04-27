

theta = @(n) linspace(0, 2*pi-2*pi/n, n)';
x = @(r,n) r.*[cos(theta(n)), sin(theta(n))];

bdr = x(100*meter, 50);
xwi = x(10*meter, 4);
xwp = x(50*meter, 8);

G = pebiGrid2D(100/10, [100,100], ...
    'cellConstraints', num2cell([xwi; xwp], 2)', ...
    'polyBdr', bdr, ...
    'CCrefinement', true, ...
    'CCFactor', 0.2);

G = computeGeometry(G);

%%
subGroups = true;

rock = makeRock(G, 100*milli*darcy, 0.1);
rock = addThermalRockProps(rock);
fluid = initSimpleADIFluid('phases', 'W', 'n', 1, 'rho', 1, 'mu', 1);
fluid = addThermalFluidProps(fluid, 'useEOS', true);

model = GeothermalModel(G, rock, fluid);

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

W = addThermalWellProps(W, 'T', convertFromCelcius(80));

if subGroups
    injGroup = struct('name', 'Injectors', ...
                      'children', {{'Injectors-1', 'Injectors-2'}}, ...
                      'type', 'rate', ...
                      'val', rate, ...
                      'T', convertFromCelcius(80));

    injGroup1 = struct('name', 'Injectors-1', ...
                      'type', 'rate', ...
                      'val', 'Injectors', ...
                      'frac', 0.5, ...
                      'T', nan);

    injGroup2 = struct('name', 'Injectors-2', ...
                      'type', 'rate', ...
                      'val', 'Injectors', ...
                      'frac', 0.5, ...
                      'T', nan);
else
    injGroup = struct('name', 'Injectors', ...
                      'type', 'rate', ...
                      'val', rate, ...
                      'T', convertFromCelcius(80));
end

prodGroup = struct('name', 'Producers', ...
                  'type', 'rate', ...
                  'val', 'Injectors', ...
                  'T', nan);

if subGroups
    groups = {injGroup, injGroup1, injGroup2, prodGroup};
else
    groups = {injGroup, prodGroup};
end
              
dt = rampupTimesteps(time, 10*day);

bc = addBC([], boundaryFaces(G), 'pressure', 50*barsa);
bc = addThermalBCProps(bc, 'T', convertFromCelcius(20));

schedule = simpleSchedule(dt, 'W', W, 'bc', bc);
schedule.control(1).groups = groups;

state0 = initResSol(G, 50*barsa, 1);
state0.T = repmat(convertFromCelcius(20), G.cells.num, 1);

%%
problem = packSimulationProblem(state0, model, schedule, 'group-test');
simulatePackedProblem(problem, 'restartStep', 1);

%%
close all
[wellSols, states, reports] = getPackedSimulatorOutput(problem);
[groupSols, wellSols] = computeWellGroupSols(wellSols, schedule.control(1).W);
% wellSols = computeWell
plotToolbar(model.G, states); axis equal tight; colormap(hot);
plotWellSols(wellSols);