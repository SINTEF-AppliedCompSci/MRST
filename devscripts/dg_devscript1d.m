mrstModule add dg ad-props ad-core ad-blackoil blackoil-sequential ...
    mrst-gui matlab_bgl vemmech vista

%%

setup = getDGTestCase('simple1d');

%%

[ws, st, rep] = simulateScheduleAD(setup.state0, setup.modelDG{1}, setup.schedule);

%%


gravity reset off
% gravity reset on; gravity([0,-9.81]);

n = 40;
l = 1000*meter;
G = computeGeometry(cartGrid([n,1], [1, 0.01]*l));
G = computeVEMGeometry(G);
G = computeCellDimensions(G);

rock  = makeRock(G, 100*milli*darcy, 1);
fluid = initSimpleADIFluid('phases', 'WO'                   , ...
                           'rho'   , [1000, 1]*kilogram/meter^3, ...
                           'mu'    , [0.5, 0.5]*centi*poise     , ...
                           'n'     , [2, 2]                 );

modelfi = TwoPhaseOilWaterModel(G, rock, fluid);
modelFV = getSequentialModelFromFI(modelfi);
modelDG = modelFV;

%%

time = 5*year;
rate = 1.2*sum(poreVolume(G, rock))/time;
W = [];
W = addWell(W, G, rock, 1          , 'type', 'rate', 'val', rate    , 'comp_i', [1,0]);
W = addWell(W, G, rock, G.cells.num, 'type', 'bhp' , 'val', 50*barsa, 'comp_i', [1,0]);

dt    = 30*day;
dtvec = rampupTimesteps(time, dt, 0);
schedule = simpleSchedule(dtvec, 'W', W, 'src', src);

sW     = 0.0;
state0 = initResSol(G, 100*barsa, [sW,1-sW]);

%%

degree = [0, 1, 2, 3];
% degree = [1,2];

% degree = 1;

[jt, ot, mt] = deal(Inf);
% 
jt = 0.2;
% mt = 0.0;
ot = 0.01;

[wsDG, statesDG] = deal(cell(numel(degree),1));
for dNo = 1:numel(degree)
    disc = DGDiscretization(G, rock, ...
                                    'degree', degree(dNo), ...
                                    'basis'              , 'legendre', ...
                                    'useUnstructCubature', false, ...
                                    'jumpTolerance', jt, ...
                                    'outTolerance', ot, ...
                                    'meanTolerance', mt);
    modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid, ...
                                             'disc', disc, ...
                                             'dsMaxAbs', 0.2/(degree(dNo)+1));

    state0 = assignDofFromState(modelDG.transportModel.disc, state0);
    [wsDG{dNo}, statesDG{dNo}, rep] = simulateScheduleAD(state0, modelDG, schedule);
end

%%

close all

stateNo = round(linspace(1,numel(schedule.step.val), 5));

x = G.cells.centroids(:,1);
for sNo = stateNo
    figure
    hold on
    for dNo = 1:numel(degree)
        plot(x, statesDG{dNo}{sNo}.s(:,1), 'linew', 2)
    end
    hold off
end