mrstModule add dg vem vemmech ad-props ad-core ad-blackoil blackoil-sequential gasinjection mrst-gui reorder matlab_bgl

%%

gravity reset off

n = 10;
l = 1000*meter;
G = computeGeometry(cartGrid([n,1], [1,0.1]*l));
G = computeVEMGeometry(G);
G = computeCellDimensions(G);

rock = makeRock(G, 100*milli*darcy, 1);
fluid = initSimpleADIFluid('phases', 'WO'                   , ...
                           'rho'   , [1000, 1]*kilogram/meter^3, ...
                           'mu'    , [0.5, 0.5]*centi*poise     , ...
                           'n'     , [1, 1]                 );

modelfi = TwoPhaseOilWaterModel(G, rock, fluid);
modelFV = getSequentialModelFromFI(modelfi);
modelDG = modelFV;

%%

time = 2*year;
rate = 2*sum(poreVolume(G, rock))/time;
W = [];
W = addWell(W, G, rock, 1          , 'type', 'rate', 'val', rate    , 'comp_i', [1,0]);
W = addWell(W, G, rock, G.cells.num, 'type', 'bhp' , 'val', 50*barsa, 'comp_i', [1,0]);

dt    = 10*day;
dtvec = rampupTimesteps(time, dt, 0);
schedule = simpleSchedule(dtvec, 'W', W, 'src', src);

sW     = 0.0;
state0 = initResSol(G, 100*barsa, [sW,1-sW]);

%%

degree = [0, 1, 2, 3, 4, 5];
jt = 0.1;
mt = 0.0;
ot = 0.0;

nls = NonLinearSolver('maxIterations', 50, 'useLinesearch', true);
[wsDG, statesDG, disc] = deal(cell(numel(degree),1));
for dNo = 1:numel(degree)
    disc{dNo} = DGDiscretization(modelDG.transportModel                    , ...
                                    'degree'             , degree(dNo), ...
                                    'basis'              , 'legendre' , ...
                                    'useUnstructCubature', false      , ...
                                    'jumpTolerance'      , jt         , ...
                                    'outTolerance'       , ot         , ...
                                    'meanTolerance'      , mt         , ...
                                    'plotLimiterProgress', false       );
    modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid  , ...
                                       'disc'    , disc{dNo}          , ...
                                       'dsMaxAbs', Inf                , ...
                                       'nonlinearTolerance', 1e-3     );
    modelDG.pressureModel = PressureOilWaterModelSemiDG(G, rock, fluid, ...
                                       'disc', disc{dNo}              );
    modelDG.transportNonLinearSolver = nls;
                                   
    state0 = assignDofFromState(modelDG.transportModel.disc, state0);
    [wsDG{dNo}, statesDG{dNo}, rep] ...
        = simulateScheduleAD(state0, modelDG, schedule);
    
end

%%

[wsFV, statesFV, rep] = simulateScheduleAD(state0, modelFV, schedule);

%%

dsnDG = cellfun(@(d) ['dG(' num2str(d), ')'] , num2cell(degree), 'unif', false);

close all

for dNo = 1:numel(degree)
    figure('name', dsnDG{dNo});
    plotToolbar(G, statesDG{dNo}, 'plot1d', true);
    colormap(jet);
end

figure('name', 'FV');
plotToolbar(G, statesFV, 'plot1d', true);
colormap(jet)

dsn = horzcat('FV', dsnDG);

plotWellSols({wsFV, wsDG{:}}, schedule.step.val, 'datasetNames', dsn)
% plotWellSols({wsFV, wsDG{:}, wsDGReorder}, schedule.step.val)

%%

close all

ll = cellfun(@(d) ['dG(', num2str(d), ')'], num2cell(degree), 'unif', false);

cNo = round(G.cartDims(1)/2);
s = zeros(numel(schedule.step.val),numel(degree));
for dNo = 1:numel(degree)
    s(:,dNo) = cellfun(@(s) s.s(cNo,1), statesDG{dNo});
end

plot(s, 'linew', 2);
legend(ll, 'location', 'southeast');
    
%%
close all
figure('position', [-1000, 0, 800, 600])
% h = zeros(numel(degree),1);
saturation = cell(numel(degree), 1);
hold on
for dNo = 1:numel(degree)
    [h(dNo), saturation{dNo}, coords, keep, n] = ...
        plotSaturationDG(disc{dNo}, statesDG{dNo}{1}, 'plot1d', true, 'linew', 2);
end
ll = cellfun(@(d) ['dG(', num2str(d), ')'], num2cell(degree), 'unif', false);
legend(ll, 'location', 'southwest');
hold off
axis tight; box on
ylim([-0.1,1.1])

for sNo = 1:numel(statesDG{1})
    for dNo = 1:numel(degree)
        s = saturation{dNo}(statesDG{dNo}{sNo});
        set(h(dNo), 'YData', s);
    end
    pause(0.5);
end