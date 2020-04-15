% Quarter five-spot example on Cartesian grid
mrstModule add ad-core ad-props ad-blackoil coarsegrid
mrstVerbose on

%% Model
n     = 105;
G     = computeGeometry(cartGrid([n,n,1], [1000,1000,10]*meter));
rock  = makeRock(G, 100*milli*darcy, 0.4);
x = G.cells.centroids(:,1:2);
p = partitionCartGrid(G.cartDims, [2,2,1]);
[ii, jj] = gridLogicalIndices(G);
m   = floor(n/2)+1;
mid = (ii == m) | (jj == m);
partition = 1*(ii > m) + 2*(jj > m) + 1;
partition(mid) = 5;
k = zeros(G.cells.num,1);
for i = 1:4
    cells = partition == i;
    xl = x(cells,:);    
    xl = (xl - mean(xl))./500*6;
    if i == 2 || i == 4
        xl(:,1) = -xl(:,1);
    end
    if i > 2
        xl(:,2) = -xl(:,2);
    end
    k(cells) = peaks(xl(:,1), xl(:,2))/3;
end
lperm = log10(rock.perm);
lperm = lperm + k;
rock.perm = 10.^lperm;

fluid = initSimpleADIFluid('phases', 'WO'                         , ...
                           'n'     , [2,2]                        , ...
                           'mu'    , [1,2]*centi*poise            , ...
                           'rho'   , [1000,800]*kilogram/(meter^3), ...
                           'c'     , [1e-6, 1e-5]/barsa           );
modelFI = GenericBlackOilModel(G, rock, fluid, 'gas', false);

%% Schedule
% Wells
time = 2*year;
rate = sum(poreVolume(G, rock))/time;
bhp  = 50*barsa;
injector = @(W,i,j) verticalWell(W, G, rock, i, j, [], ...
                                     'type'  , 'rate', ...
                                     'val'   , rate  , ...
                                     'comp_i', [1,0] );
producer = @(W,i,j) verticalWell(W, G, rock, i, j, [], ...
                                      'type'  , 'bhp', ...
                                      'val'   , bhp  , ...
                                      'comp_i', [1,0]);
W = [];
W = injector(W, m, m);
W = producer(W, 1, 1);
W = producer(W, n, 1);
W = producer(W, n, n);
W = producer(W, 1, n);
% Schedule
schedule = simpleSchedule(rampupTimesteps(time, 30*day), 'W', W);

%% Partition domain
mrstModule add coarsegrid
% Partition
partition = partitionCartGrid(G.cartDims, [5,5,1]);
% Initial state
state0    = initResSol(G, bhp, [0,1]);

%% Construct submodel
mrstModule add ddc
cells = partition == 13;
submodel = SubdomainModel(modelFI, cells);

%% Plot submodel
close all
colors = lines(2);
mappings = submodel.mappings;
plotGrid(modelFI.G, 'faceAlpha', 0.1, 'edgeAlpha', 0.1);
plotGrid(modelFI.G, mappings.cells.keep    , 'faceColor', colors(2,:), 'edgeAlpha', 0.1);
plotGrid(modelFI.G, mappings.cells.internal, 'faceColor', colors(1,:));
axis equal tight

%% Simulate submodel
[substate0, mappings] = getSubState(state0, mappings);
subschedule = getSubSchedule(schedule, mappings);
[subsWellSols, substates, subreports] = simulateScheduleAD(substate0, submodel, subschedule);

%% Inspect submodel results
mrstModule add mrst-gui
% states = cellfun(@(substate) mapState(state0, substate, mappings, 'mapWellSol', false), substates, 'UniformOutput', false);
% plotToolbar(G, states); axis equal tight

%% Simulate entire problem with sequential splitting
mrstModule add blackoil-sequential
pmodel   = PressureModel(modelFI);
tmodel   = TransportModel(modelFI);
modelSeq = SequentialPressureTransportModel(pmodel, tmodel);
[wsSeq, statesSeq, reportsSeq] = simulateScheduleAD(state0, modelSeq, schedule);

%% Simulate with domain decomposition in transport
% tmodelDD = ParallelDomainDecompositionModel(tmodel, partition);
tmodelDD = DomainDecompositionModel(tmodel, partition, 'parallel', true);
modelseqDD = SequentialPressureTransportModel(pmodel, tmodelDD, 'parentModel', modelFI);
[wsSeqDD, statesSeqDD, reportsSeqDD] = simulateScheduleAD(state0, modelseqDD, schedule);

%% Plot results
close all
statesDiff = cellfun(@(s1, s2) compareFields(s1, s2, 'relative', false), statesSeq, statesSeqDD);
figure(), plotToolbar(G, statesDiff); axis equal tight
figure(), plotToolbar(G, statesSeqDD); axis equal tight

%%
its   = getIterations(reportsSeq.ControlstepReports, 'solver', 'TransportSolver');
itsDD = getIterations(reportsSeqDD.ControlstepReports, 'solver', 'TransportSolver');
t = cumsum(schedule.step.val)/day;
figure();
pargs = {'LineWidth', 2};
hold on
plot(t, cumsum(its.total), pargs{:})
plot(t, cumsum(itsDD.total), pargs{:})
hold off
xlim([0,t(end)])
xlabel('Time (days)')
title('Cumulative transport iterations')
legend({'Global', 'Domain decomposition'}, 'location', 'east')
box on;

%%
mrstModule add matlab_bgl
cmodel = upscaleModelTPFA(modelFI, partition);
st = upscaleState(cmodel, modelFI, statesSeq{end});
v = sum(st.flux,2);
v = v(cmodel.operators.internalConn);
order = getTopologicalPermutation(cmodel.G, v);
partition2 = order(partition);

%%
tmodelDD = DomainDecompositionModel(tmodel, partition2, 'strategy', 'multiplicative');
modelseqDD = SequentialPressureTransportModel(pmodel, tmodelDD, 'parentModel', modelFI);
[wsSeqDD, statesSeqDD, reportsSeqDD] = simulateScheduleAD(state0, modelseqDD, schedule);

%%

getPartition = @(varargin) getTopologicalFluxPartition(varargin{:}, 'blockSize', 50);
tmodelDD = DomainDecompositionModel(tmodel, partition, 'strategy', 'multiplicative', 'getPartition', getPartition);
tmodelDD.verboseSubmodel = true;
modelseqDD = SequentialPressureTransportModel(pmodel, tmodelDD, 'parentModel', modelFI);
[wsSeqDD, statesSeqDD, reportsSeqDD] = simulateScheduleAD(state0, modelseqDD, schedule);

%%

plotToolbar(G, statesSeqDD);