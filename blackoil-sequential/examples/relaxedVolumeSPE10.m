mrstModule add mrst-gui ad-props ad-core ad-blackoil blackoil-sequential msrsb spe10 coarsegrid
%% Set up simulation problem
totTime = 2000*day;

layers = 5;
nLayers = numel(layers);

% Offsets to move the wells if needed
ofs = [0, 0];
wloc     = [  1 + ofs(1),   60 - ofs(1),     1 + ofs(1),   60 - ofs(1),   30 ;
              1 + ofs(2),    1 + ofs(2),   220 - ofs(2),  220 - ofs(2),  110 ];

[G, W, rock] = getSPE10setup(layers, wloc);

mp = 0.01;
rock.poro(rock.poro < mp) = mp;

pv = poreVolume(G, rock);

% Either two-phase or three-phase
if 1
    compi = [.5, 0, .5];
    sat = [0, 1, 0];
    getmodel = @ThreePhaseBlackOilModel;
else
    compi = [1, 0];
    sat = [0, 1];
    getmodel = @TwoPhaseOilWaterModel;
end

for i = 1:numel(W)
    isinj = W(i).val > 400*barsa;
    if isinj
        W(i).val = sum(pv)/totTime;
        W(i).type = 'rate';
    else
        W(i).val = 4000*psia;
        W(i).type = 'bhp';
    end
    W(i).compi = compi;
    W(i).sign = 1 - 2*~isinj;
end

% Set up model with quadratic relperms and constant compressibility in oil
p0 = 300*barsa;
c = [1e-6, 1e-4, 1e-3]/barsa;
% c = [0,0,0];
fluid = initSimpleADIFluid('mu', [1, 5, 1]*centi*poise, ...
                           'c', c, ...
                           'rho', [1000, 700, 1], 'n', [2 2 2]);

% Fully implicit model
modelfi = getmodel(G, rock, fluid);
% Sequential pressure-transport model with same type


figure;
plotCellData(G, log10(rock.perm(:, 1)), 'edgec', 'none');
plotWell(G, W)
% Set up the schedule
nstep = 100;

dt = rampupTimesteps(totTime, totTime/nstep);
clear schedule;
schedule = simpleSchedule(dt, 'W', W);

state = initResSol(G, p0, sat);


%% Solve a sequential solve
% We first solve a pressure equation, fix total velocity and solve for n-1
% phases. To find the remaining saturation, we use the relation 
%
%    s_o + s_w + s_g = 1.
% 
% We naturally incur a bit of mass-balance error for the omitted phase
% pseudocomponent.
model = getSequentialModelFromFI(modelfi);
model.transportModel.useCNVConvergence = true;

[ws_seq, states_seq, rep_seq] = simulateScheduleAD(state, model, schedule);
%% Solve a second sequential solve
% We first solve a pressure equation, fix total velocity and solve for all
% phases/components. We relax the sum of saturations assumption to achieve
% this.
model.transportModel.conserveOil = true;
model.transportModel.conserveWater = true;
model.transportModel.useCNVConvergence = true;
[ws_seq2, states_seq2, rep_seq2] = simulateScheduleAD(state, model, schedule);
%% Fully-implicit
[wsfi, statesfi] = simulateScheduleAD(state, modelfi, schedule);
%% Plot well curves for fully implicit problem
plotWellSols({ws_seq, ws_seq2, wsfi}, cumsum(schedule.step.val), 'datasetnames', {'Sequential (relaxed mass)', 'Sequential (relaxed volume)', 'FI'})
