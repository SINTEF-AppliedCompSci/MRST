N = 50;
l = 500;
G = computeVEMGeometry(cartGrid([N,1], [l, l/N]*meter));
G = computeCellDimensions(G);

rock = makeRock(G, 100*milli*darcy, 1);
n = 3;
fluid = initSimpleADIFluid('phases', 'WO'                  , ...
                           'rho'   , [1,1]*kilogram/meter^3, ...
                           'mu'    , [1,1]*centi*poise     , ...
                           'n'     , [n,n]);

modelFI = TwoPhaseOilWaterModel(G, rock, fluid);
modelDG = getSequentialModelFromFI(modelFI);

%%

time = 2*year;
rate = 1*sum(poreVolume(G, rock))/time;
radius = 0.05;
bhp = 400*barsa;
W = [];
W = addWell(W, G, rock, 1          , 'type', 'bhp', 'val', bhp, 'comp_i', [1,0], 'radius', radius);
W = addWell(W, G, rock, G.cells.num, 'type', 'bhp', 'val', 0*barsa  , 'comp_i', [1,0], 'radius', radius);

dt    = 30*day;
dtvec = rampupTimesteps(time, dt, 0);

schedule = simpleSchedule(dtvec, 'W', W);

sW     = 0.0;
state0 = initResSol(G, 100*barsa, [sW,1-sW]);

%%

degree = [2];
% statesDG = cell(numel(degree),1);
for dNo = 1:numel(degree)
    disc = DGDiscretization(modelDG.transportModel, 1, 'degree', degree(dNo), 'basis', 'legendre', 'limiter', 'cap');
    modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid, 'disc', disc);    

    state0 = assignDofFromState(modelDG.transportModel.disc, state0);
    [ws, statesDG{dNo}, rep] = simulateScheduleAD(state0, modelDG, schedule);
end

%%

NF   = 500;
GF   = computeVEMGeometry(cartGrid([NF,1], [l, l/NF]*meter));
rockF = makeRock(GF, 100*milli*darcy, 1);
modelFV = TwoPhaseOilWaterModel(GF, rockF, fluid);

WF = [];
WF = addWell(WF, GF, rockF, 1          , 'type', 'bhp', 'val', bhp, 'comp_i', [1,0], 'radius', radius);
WF = addWell(WF, GF, rockF, GF.cells.num, 'type', 'bhp', 'val', 0*barsa  , 'comp_i', [1,0], 'radius', radius);

dt    = 5*day;
dtvec = rampupTimesteps(time, dt, 0);

scheduleF = simpleSchedule(dtvec, 'W', WF);

state0F = initResSol(GF, 100*barsa, [sW,1-sW]);

[ws, statesF, rep] = simulateScheduleAD(state0F, modelFV, scheduleF);

%%

close all
figure('Position', [-2000,0,1500, 1000])

nsteps = 5;
steps = round(linspace(1,numel(schedule.step.val)-3,nsteps));
tvec = cumsum(schedule.step.val);

tvecF = cumsum(scheduleF.step.val);
clr = copper(numel(steps));

mrksz = [8, 8, 8];
mrks = {'--sq', '--o', '--^'};
lw = 1.5;

x = linspace(0,l,N)';
xF= linspace(0,l,NF)';

hold on
for sNo = 1:numel(steps)
    t = tvec(steps(sNo));
%     sE = exactBL(modelDG.transportModel, statesDG{1}{steps(sNo)});
%     plot(xf, sE(xf, t), 'color', clr(sNo,:), 'linew', lw);
    ix = abs(t - tvecF);
    ix = ix == min(ix);
    plot(xF, statesF{ix}.s(:,1), 'color', clr(sNo,:), 'linew', lw);
    for dNo = 1:numel(degree)
        plot(x, statesDG{dNo}{steps(sNo)}.s(:,1), mrks{dNo}, 'color', clr(sNo,:), 'linew', lw, 'markerfacec', 'w');
    end
end
hold off