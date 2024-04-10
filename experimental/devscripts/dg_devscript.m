mrstModule add vem dg vemmech ad-props ad-core ad-blackoil blackoil-sequential gasinjection

%%

n = 50;
l = 1000;
G = computeGeometry(cartGrid([n,1], [l,10]*meter));
G = computeVEMGeometry(G);
G = computeCellDimensions(G);

rock = makeRock(G, 100*milli*darcy, 1);
fluid = initSimpleADIFluid('phases', 'WO'                   , ...
                           'rho'   , [1, 1]*kilogram/meter^3, ...
                           'mu'    , [1, 1]*centi*poise     , ...
                           'n'     , [1,1]                 );

modelfi = TwoPhaseOilWaterModel(G, rock, fluid);
modelFV = getSequentialModelFromFI(modelfi);
modelDG = modelFV;
                       
%%

time = 2*year;
rate = 1*sum(poreVolume(G, rock))/time;
W = [];
W = addWell(W, G, rock, 1          , 'type', 'rate', 'val', rate    , 'comp_i', [1,0]);
% W = addWell(W, G, rock, 1          , 'type', 'bhp', 'val', 2000*barsa, 'comp_i', [1,0]);
W = addWell(W, G, rock, G.cells.num, 'type', 'bhp' , 'val', 50*barsa, 'comp_i', [1,0]);

dt    = 15*day;
dtvec = rampupTimesteps(time, dt, 0);

schedule = simpleSchedule(dtvec, 'W', W);

sW     = 0.0;
state0 = initResSol(G, 100*barsa, [sW,1-sW]);

%%

degree = [0, 1, 2];
[wsDG, statesDG] = deal(cell(numel(degree),1));
for dNo = 1:numel(degree)
    disc    = DGDiscretization(modelDG.transportModel, 1, 'degree', degree(dNo), 'basis', 'legendre');%, 'limiter', 'none');
    
%     disc.dofPos = disc.updateDofPos();
    modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid, 'disc', disc);    

    state0 = disc.assignDofFromState(state0);
%     state0.degree = repmat(disc.degree, G.cells.num, 1);
    [wsDG{dNo}, statesDG{dNo}, rep] = simulateScheduleAD(state0, modelDG, schedule);
end

%%

% limiter = {'cap', 'none', 'tvb'};
% for lNo = 1:2
%     disc    = DGDiscretization(modelDG.transportModel, 1, 'degree', 1, 'basis', 'legendre', 'limiter', limiter{lNo});
%     modelDG.transportModel = TransportOilWaterModelDG(G, rock, fluid, 'disc', disc);    
% 
%     state0 = assignDofFromState(modelDG.transportModel.disc, state0);
%     [ws, statesDG{lNo}, rep] = simulateScheduleAD(state0, modelDG, schedule);
% end
    
%%

[wsFV, statesFV, rep] = simulateScheduleAD(state0, modelFV, schedule);

%%

figure('Position', [0,0,1500,600])
x = linspace(0,l,n);

steps = round(linspace(1, numel(schedule.step.val)-5,5));
% steps = [2,12,20]
clr = lines(numel(statesDG)+1);
clr = copper(numel(steps));
[h, hDG] = deal([]);
mrksz = [8, 8, 8];
mrks = {'-', 'o-', '-^'};
lw = 1.5;

for sNo = 1:numel(steps)
%     subplot(1,numel(steps), sNo)
    hold on
    hFV = plot(x, statesFV{steps(sNo)}.s(:,1), '--', 'linew', 4, 'color', clr(sNo,:));
    for dNo = 1:numel(degree)
        hDG(dNo) = plot(x, statesDG{dNo}{steps(sNo)}.s(:,1), mrks{dNo}, 'markers', mrksz(dNo), 'linew', lw, 'color', clr(sNo, :), 'markerfacecolor', clr(sNo,:));
    end
    if isempty(h)
%         h = [hFV, hDG];
%     xlabel('Distance from injector');
    end
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

% ds = 0.1;
% ylim([-ds, 1+ds]);
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
plotToolbar(G, statesDG{2}, 'plot1d', true);

figure
plotToolbar(G, statesFV, 'plot1d', true);

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
