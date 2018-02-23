mrstModule add dg vem vemmech ad-props ad-core ad-blackoil blackoil-sequential gasinjection

%%

n = 100;
l = 1000;
G = computeGeometry(cartGrid([n,1], [l,10]*meter));
% G = computeGeometry(cartGrid([n,n], [l,l]*meter));
G.nodes.coords = G.nodes.coords;
G = computeVEMGeometry(G);

rock = makeRock(G, 100*milli*darcy, 1);
fluid = initSimpleADIFluid('phases', 'WO'                        , ...
                           'rho'   , [1, 1]*kilogram/meter^3, ...
                           'mu'    , [1, 1]*centi*poise                 , ...
                           'n'     , [1, 1]                      );

modelfi = TwoPhaseOilWaterModel(G, rock, fluid);
modelFV = getSequentialModelFromFI(modelfi);
modelDG = modelFV;
                       
%%

time = 2*year;
rate = 1*sum(poreVolume(G, rock))/time;
W = [];
W = addWell(W, G, rock, 1          , 'type', 'rate', 'val', rate    , 'comp_i', [1,0]);
W = addWell(W, G, rock, G.cells.num, 'type', 'bhp' , 'val', 50*barsa, 'comp_i', [1,0]);

dt    = 30*day;
dtvec = rampupTimesteps(time, dt, 0);

schedule = simpleSchedule(dtvec, 'W', W);

sW     = 0.0;
state0 = initResSol(G, 100*barsa, [sW,1-sW]);

%%

degree = [0, 1, 2];
states = cell(numel(degree),1);
for dNo = 1:numel(degree)
    disc    = DGDiscretization(modelDG.transportModel, G.griddim, 'degree', degree(dNo), 'basis', 'legendre');
    modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid, 'disc', disc);    

    state0 = assignDofFromState(modelDG.transportModel.disc, state0);
    [ws, states{dNo}, rep] = simulateScheduleAD(state0, modelDG, schedule);
end

%%

[ws, statesFV, rep] = simulateScheduleAD(state0, modelFV, schedule);

%%

figure('Position', [0,0,1500,600])
x = linspace(0,l,n);

steps = round(linspace(1, numel(schedule.step.val)-5, 3));
steps = [2,12,20]
clr = lines(numel(states)+1);
clr = copper(numel(steps));
[h, hDG] = deal([]);
mrksz = [8, 8, 8];
mrks = {'-o', '^-', '-sq'};
lw = 1.5;

for sNo = 1:numel(steps)
%     subplot(1,numel(steps), sNo)
    hold on
    hFV = plot(x, statesFV{steps(sNo)}.s(:,1), '--', 'linew', 4, 'color', clr(sNo,:));
    for dNo = 1:numel(degree)
        hDG(dNo) = plot(x, states{dNo}{steps(sNo)}.s(:,1), mrks{dNo}, 'markers', mrksz(dNo), 'linew', lw, 'color', clr(sNo, :), 'markerfacecolor', clr(sNo,:));
    end
    if isempty(h)
%         h = [hFV, hDG];
    end
    ds = 0.1;
    ylim([-ds, 1+ds]);
%     xlabel('Distance from injector');
    ax = gca;
    ax.FontSize = 15;
    box on
    
end

dgNames = cellfun(@(c) ['dG(', num2str(c), ')'], num2cell(degree), 'unif', false);
lgnd = {'FV', dgNames{:}};
legend(h, lgnd, 'location', 'northeast')
% 
% yyaxis left
% lgnd = cellfun(@(ts) ['Timestep ', num2str(ts)], mat2cell(steps, 1, ones(1,numel(steps))), 'unif', false);
% legend(hT, lgnd);

ds = 0.1;
ylim([-ds, 1+ds]);
% xlabel('Distance from injector');
% ylabel('Water saturation');
ax = gca;
ax.FontSize = 15;
box on

%%

pth = mrstPath('dg');
print([pth, '/', 'dgExample1D'], '-dpng', '-r300');

%%

figure;
plotToolbar(G, state, 'plot1d', true);

figure
plotToolbar(G, state2, 'plot1d', true);