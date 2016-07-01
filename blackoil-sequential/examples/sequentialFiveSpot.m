%%
mrstModule add ad-core ad-blackoil blackoil-sequential spe10
mrstModule add ad-blackoil ad-core libgeometry spe10

cartDims = [50, 50, 1];

physDims = [1000, 1000, 1];
G = cartGrid(cartDims, physDims);
G = computeGeometry(G);
gravity reset on

rng(0);
% K = logNormLayers(cartDims);

phi = gaussianField(cartDims, [0.05, 0.3], 3, 8);
phi = phi(:);
A = 1e5;
K = phi.^3./((A.^2).*8*(1 - phi).^2);
rock = makeRock(G, K, phi);

pv = poreVolume(G, rock);
T = 30*year;
irate = sum(pv)/(T*4);

makeInj = @(W, name, I, J, compi) verticalWell(W, G, rock, I, J, [],...
    'Name', name, 'radius', 5*inch, 'sign', 1, 'Type', 'rate',...
    'Val', irate, 'comp_i', compi);
W = [];

pureGas = [0, 0, 1];
pureWat = [1, 0, 0];
initSat = [0, 1, 0];

% 
% W = makeInj(W, 'I1', 1,           1,           pureWat);
% W = makeInj(W, 'I3', cartDims(1), cartDims(2), pureGas);
% W = makeInj(W, 'I4', 1,           cartDims(2), pureWat);
% W = makeInj(W, 'I2', cartDims(1), 1,           pureGas);

W = makeInj(W, 'I1', 1,           1,           []);
W = makeInj(W, 'I3', cartDims(1), cartDims(2), []);
W = makeInj(W, 'I4', 1,           cartDims(2), []);
W = makeInj(W, 'I2', cartDims(1), 1,           []);

I = ceil(cartDims(1)/2);
J = ceil(cartDims(2)/2);


W = verticalWell(W, G, rock, I, J, [], 'Name', 'P1', 'radius', 5*inch, ...
    'Type', 'bhp', 'Val', 100*barsa, 'comp_i', [1, 1, 1]/3, 'Sign', -1);

[W_water, W_gas] = deal(W);
for i = 1:numel(W)
    if W(i).sign < 0
        % Skip producer
        continue
    end
    W_water(i).compi = [1, 0, 0];
    W_gas(i).compi   = [0, 0, 1];
end

dT_target = 90*day;
dt = rampupTimesteps(T, dT_target, 10);

schedule = struct();
schedule.control = [struct('W', W_water);... % Water control 1
                    struct('W', W_gas)]; % Gas control 2

schedule.step.val = dt;
schedule.step.control = (mod(cumsum(dt), 2*dT_target) > dT_target) + 1;

% schedule = simpleSchedule(dt, 'W', W);

% schedule.control(2) = schedule.control(1);
% schedule.control(2).W = W_gas;
% schedule.control(1).W = W_water;


fluid = initSimpleADIFluid('phases',    'WOG', ...
                           'rho',       [1000, 700, 250], ...
                           'n',         [2, 2, 2], ...
                           'c',         [0, 1e-4, 1e-3]/barsa, ...
                           'mu',        [1, 4, 0.25]*centi*poise ...
                           );

model = ThreePhaseBlackOilModel(G, rock, fluid);
state = initResSol(G, 0, initSat);


%%
figure; plotToolbar(G, rock)

%%
seqModel = getSequentialModelFromFI(model);

%%
timer = tic();
[wsSeq, statesSeq] = simulateScheduleAD(state, seqModel, schedule);
t_split = toc(timer);
%% Run the entire schedule
% solver.timeStepSelector.reset();

timer = tic();
[wsFIMP, statesFIMP] = simulateScheduleAD(state, model, schedule);
t_fi = toc(timer);

%%
plotWellSols({wsSeq, wsFIMP}, cumsum(schedule.step.val), 'datasetnames', {'Sequential', 'FIMP'}, 'linestyles', {'-'})




%% Compute solution with refined time steps
% We will now compute a solution with refined time steps. As the time-steps
% become smaller, the solution becomes more accurate. In order to achieve
% increased accuracy without manually changing the timesteps, we can use
% a automatic timestep selector based on saturation change targets. We let
% the solver aim for a maximum saturation change of 1% in each cell during
% the timesteps to get very small steps.
stepSel = StateChangeTimeStepSelector(...
          'targetProps', 's',...            % Saturation as change target
          'targetChangeAbs', 0.01,...       % Target change of 0.01
          'firstRampupStepRelative', 0.01); % Initial rampup step is dt0/100
solver = NonLinearSolver('timeStepSelector', stepSel);
% Simulate with small timesteps. As the resulting timesteps will be very
% small, it will take some time (about six minutes on the workstation where
% the example was written).
[wsFine, statesFine, repFine] = simulateScheduleAD(state, model, schedule, ...
                        'nonlinearsolver', solver, 'outputMinisteps', true);
