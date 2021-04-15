mrstModule add ad-core ad-eor ad-blackoil ad-props sequential matlab_bgl

gravity reset on

n = 100;
G = computeGeometry(cartGrid([n,1,1], [1000,1,1]));
rock = makeRock(G, 100*milli*darcy, 1);

%%

fluid = initSimpleADIFluid('n'     , [2, 2, 2], ...
                           'rho'   , [1000, 800, 100]*kilogram/meter^3, ...
                           'phases', 'WOG', ...
                           'mu'    , [1, 10, 2]*centi*poise);

sOres_i= 0.2;
fluid = addSolventProperties(fluid, 'n', 2, ...
                                    'rho', 100*kilogram/meter^3, ...
                                    'mixPar', 1, ...
                                    'mu'    , 1*centi*poise, ...
                                    'sOres_i', sOres_i, ...
                                    'sOres_m', 0.0);
                                
model = FourPhaseSolventModel(G, rock, fluid);
model.extraStateOutput = true;

T = 2*year;
rate = 0.2*sum(poreVolume(G, rock))/T;
W = [];
W = addWell(W, G, rock, 1, 'type', 'rate', 'val', rate, 'comp_i', [1,0,0,0]);
W = addWell(W, G, rock, G.cells.num, 'type', 'bhp', 'val', 50*barsa, 'comp_i', [1,0,0,0]);

schedule = makeWAGschedule(W, 2, 'time', T);


sO = sOres_i + 0.0;
sG = 1e-1;
state0 = initResSol(G, 100*barsa, [1-sO-sG sO sG 0]);
state0.wellSol = initWellSolAD(W, model, state0);

nls = NonLinearSolver('useLineSearch', false);


%%

[ws, states, reports] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nls);

%%

mrstModule add mrst-gui

figure; clf
plotToolbar(G, states, 'plot1d', true)
ylim([0 1]);

%%

%%

% clc

pv = poreVolume(G, rock);
n = numel(ws);
nw = numel(W_G);
[resM, wellM] = deal(zeros(n,4));
rhoSurf = model.getSurfaceDensities;

wc = vertcat(W_G.cells);
for sNo = 1:n    
    
    w = ws{sNo};
    dt = schedule.step.val(sNo);
    for wNo = 1:numel(w)
        wellM(sNo,:) = wellM(sNo,:) ...
                      + dt.*[w(wNo).qWs, w(wNo).qOs, w(wNo).qGs, w(wNo).qSs].*rhoSurf;
    end
    for phNo = 1:4
        resM(sNo,phNo) = sum(states{sNo}.s(:,phNo).*states{sNo}.rho(:,phNo).*pv);
    end
    
end

wellM = cumsum(wellM,1);
totM = resM - wellM;

close all

hold on
clr = lines(4);
for phNo = 1:4
    plot(resM(:,phNo), '-', 'color', clr(phNo,:), 'linew', 1);
    plot(wellM(:,phNo), '--', 'color', clr(phNo,:), 'linew', 3);
end
hold off

figure
plot(totM);

%%

x = 1:n;
l = lines(3);
clr = [l(1,:); 0,0,0; l(2:end,:)];
for sNo = 1:numel(states)
    clf
    
    s = states{sNo}.s; 
    hold on
    for phNo = 1:4
        plot(x, s(:,phNo), 'color', clr(phNo,:))
    end
    hold off
    ylim([0,1]);
    pause(0.1);
end

%%
% 
% close all
% ns = numel(states);
% ns = 400;
% pos = [500, 500, 2000,1000];
% fig = figure('position', pos);
% M = struct('cdata',[],'colormap',[]);
% ttl = {'S_o', 'S_s'};
% nx = 1000;
% x = linspace(0,100,nx);
% 
% for i = 1:ns
%     
%     clf;
%     ll = lines(3);
%     clr = [ll(1,:); [0,0,0]; ll(2:end,:)];
%     
%     lw = 5;
%     hold on
%     ph = [2,4];
%     for j = 1:numel(ph)
%         plot(x, states{i}.s(:,ph(j)), 'color', clr(ph(j),:), 'linewidth', lw)
%     end
%     box on
%     ylim([0,1])
%     xlabel('x [m]');
%     ylabel('S [-]');
%     ax = gca;
%     ax.FontSize = 20;
%     ax.XLabel.FontSize = 45;
%     ax.YLabel.FontSize = 45;
%     legend({'Oil', 'Solvent'}, 'fontSize', 45);
%     
%     drawnow;
%     
%     ax.Units = 'pixels';
%     
%     
%     m = 0; n = 0;
%     rect = [m, n, 2000-2*m, 1000-2*n];
%     
%     M(i) = getframe(fig, rect);
%     
%     ax.Units = 'normalized';
    
% end

%%
% 
% % n = numel(states);
% n = 400;
% % n = 4;
% pth = [mrstPath('mrst-solvent'), '/presentation/figures/displacement1D/'];
% vo = VideoWriter([pth, 'displacement1D_1.avi']);
% vo.FrameRate = 1000/30;
% open(vo);
% 
% writeVideo(vo, M);
% 
% close(vo)



% M0 = totM(1,:);

% err = totM - M0;

% 
% 
% wellMtot = zeros(n,nw);
% for phNo = 1:4
%     wellMtot = wellMtot + wellM{phNo};
% end
% resMtot = sum(resM,2);
% 
% 
% wellMcum = cellfun(@(m) cumsum(m,1), wellM, 'uniformOutput', false);
% wellMtotCum = cumsum(wellMtot,1);
% 
% errTotRel = (resMtot(1) + sum(wellMtotCum,2) - resMtot)./resMtot;
% 
% errTot = abs(resMtot(1) + sum(wellMtotCum(end,:)) - resMtot(end));
% 
% err = cell(4,1);
% for phNo = 1:4
%     err{phNo} = (resM(1, phNo) + sum(wellM{phNo},2) - resM(:,phNo))./resMtot;
% end
% 
% 
% fprintf(['Absolute error: \t %.2d \n'], ...
%          errTot);
         
%%

W_G(1).name = 'I';
W_G(2).name = 'P';

figure; clf
plotGrid(G, 'edgecolor', 'none'); axis equal off
plotWell(G, W_G, 'height', 0.5)
w = [-17,14];
view([-20,10])
light('Position', [2,-3,0])

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
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
