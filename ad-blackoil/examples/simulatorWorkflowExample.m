grdecl = makeModel3([50, 50, 5], [1000, 1000, 5]*meter);
G = processGRDECL(grdecl);
G = computeGeometry(G);
%%
[ii, jj, kk] = gridLogicalIndices(G);

top = kk < G.cartDims(3)/3;
lower = kk > 2*G.cartDims(3)/3;
middle = ~(lower | top);

px = ones(G.cells.num, 1);
px(lower) = 300*milli*darcy;
px(middle) = 100*milli*darcy;
px(top) = 500*milli*darcy;

perm = [px, px, 0.1*px];

rock = makeRock(G, perm, 0.3);


%%
simTime = 10*year;
nstep = 25;
refine = 5;

pv = poreVolume(G, rock);

injRate = 1*sum(pv)/simTime;



offset = 10;
W = [];
W = verticalWell(W, G, rock, offset, offset, [],...
                'Name', 'P1', 'comp_i', [1 0], 'Val', 250*barsa, 'Type', 'bhp');
W = verticalWell(W, G, rock, G.cartDims(1) - offset, G.cartDims(2) - offset, [],...
                'Name', 'P2', 'comp_i', [1 0], 'Val', 250*barsa, 'Type', 'bhp');
W = verticalWell(W, G, rock, offset, G.cartDims(2) - offset, [], ...
                'Name', 'P3', 'comp_i', [1 0], 'Val', 250*barsa, 'Type', 'bhp');

W = verticalWell(W, G, rock, floor(G.cartDims(1)/2) + 4, floor(G.cartDims(1)/2) - 3, 1,...
                'Name', 'Injector', 'comp_i', [1 0], 'Val', injRate, 'Type', 'rate');



startSteps = repmat((simTime/nstep)/refine, refine, 1);
restSteps =  repmat(simTime/nstep, nstep - 1, 1);

timesteps = [startSteps; restSteps];

schedule = simpleSchedule(timesteps, 'Wells', W);


clf;
plotCellData(G, log10(rock.perm(:, 1)))
plotWell(G, W)
axis tight
view(50, 50)

%%
mrstModule add ad-core ad-blackoil ad-props

fluid = initSimpleADIFluid('mu',    [1, 10, 0]*centi*poise, ...
                           'rho',   [1000, 700, 0]*kilogram/meter^3, ...
                           'n',     [2, 2, 0]);

c = 0.001/barsa;
p_ref = 300*barsa;
fluid.bO = @(p) exp((p - p_ref)*c);

p0 = (100:10:500)*barsa;
plot(p0/barsa, fluid.bO(p0))
%%
model = TwoPhaseOilWaterModel(G, rock, fluid);
%%
sW = ones(G.cells.num, 1);
sW(G.cells.centroids(:, 3) < 5) = 0;

sat = [sW, 1 - sW];

state0 = initResSol(G, p_ref, sat);
%%
mrstVerbose on
[wellSols, states] = simulateScheduleAD(state0, model, schedule);


%%
mrstModule add mrst-gui
clf
plotToolbar(G, states)
view(50, 50);


simtime = cumsum(schedule.step.val);
plotWellSols(wellSols, simtime, 'field', 'qOs');

%%
mrstModule add coarsegrid
cdims = [5, 5, 2];
p0 = partitionUI(G, cdims);

clf
plotCellData(G, mod(p0, 13))
view(50, 50)
%%
G_fault = makeInternalBoundary(G, find(G.faces.tag > 0));
p = processPartition(G_fault, p0);

clf
plotCellData(G, mod(p, 13), 'EdgeColor', 'none')
plotGrid(G, p ~= p0, 'EdgeColor', 'r', 'LineWidth', 2, 'FaceColor', 'none')
view(50, 50)

%%
model_c = upscaleModelTPFA(model, p);
G_c    = model_c.G;
rock_c = model_c.rock;

schedule_c = upscaleSchedule(model_c, schedule);
state0_c = upscaleState(model_c, model, state0);

[wellSols_c, states_c] = simulateScheduleAD(state0_c, model_c, schedule_c);

%%
clf
plotToolbar(G_c, states_c)
view(50, 50);

plotWellSols({wellSols, wellSols_c}, simtime, 'DatasetNames', {'Fine scale', 'Upscaled'}, 'Field', 'qOs');
%%
mrstModule add diagnostics
D = computeTOFandTracer(states{end}, G, rock, 'Wells', schedule.control.W);

D_c = computeTOFandTracer(states_c{end}, G_c, rock_c, 'Wells', schedule_c.control.W);

%%
figure(1); clf
plotCellData(G, log(sum(D.tof, 2)));
view(50, 50);

figure(2); clf
plotCellData(G_c, log(sum(D_c.tof, 2)));
view(50, 50);

%%
figure(1); clf
plotCellData(G, D.ppart);
view(50, 50);

figure(2); clf
plotCellData(G_c, D_c.ppart);
view(50, 50);
