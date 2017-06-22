mrstModule add ad-core ad-eor ad-blackoil ad-props blackoil-sequential matlab_bgl

gravity reset on

n = 1000;
G = computeGeometry(cartGrid([n,1,1], [100,1,1]));
rock = makeRock(G, 100*milli*darcy, 1);

%%

fluid = initSimpleADIFluid('n'     , [2, 2, 2], ...
                           'rho'   , [1000, 800, 100]*kilogram/meter^3, ...
                           'phases', 'WOG', ...
                           'mu'    , [1, 10, 2]*centi*poise);

sOres_i= 0.2;
fluid = addSolventProperties(fluid, 'n', 2, ...
                                    'rho', 100*kilogram/meter^3, ...
                                    'mixPar', 0, ...
                                    'mu'    , 1*centi*poise, ...
                                    'sOres_i', sOres_i, ...
                                    'sOres_m', 0.0);
                                
model = FourPhaseSolventModel(G, rock, fluid);
model.extraStateOutput = true;

T =4*year;
rate = 1*sum(poreVolume(G, rock))/year;
[schedule, W_G, W_W] = makeWAGschedule(model, {1}, {G.cells.num}, 10, 'T', T, 'nStep', 1000, 'wRate', rate, 'gRate', rate);

s = sOres_i + 0.0;
state0 = initResSol(G, 100*barsa, [1-s s 0 0]);
state0.wellSol = initWellSolAD(W_G, model, state0);

nls = NonLinearSolver('useLineSearch', true);
nls = NonLinearSolver('useLineSearch', false);

%%

[ws, states, reports] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nls);

%%

mrstModule add mrst-gui

figure(1); clf
plotToolbar(G, states, 'plot1d', true)
ylim([0 1]);

%%

x = 1:n;
for i = 1:numel(states)
    s = states{i}.s;
    plot(x, s(:,1), x, s(:,2), x, s(:,3), x, s(:,4))
    ylim([0,1]);
    pause(0.01);
end

%%

close all
ns = numel(states);
ns = 400;
pos = [500, 500, 2000,1000];
fig = figure('position', pos);
M = struct('cdata',[],'colormap',[]);
ttl = {'S_o', 'S_s'};
nx = 1000;
x = linspace(0,100,nx);

for i = 1:ns
    
    clf;
    ll = lines(3);
    clr = [ll(1,:); [0,0,0]; ll(2:end,:)];
    
    lw = 5;
    hold on
    ph = [2,4];
    for j = 1:numel(ph)
        plot(x, states{i}.s(:,ph(j)), 'color', clr(ph(j),:), 'linewidth', lw)
    end
    box on
    ylim([0,1])
    xlabel('x [m]');
    ylabel('S [-]');
    ax = gca;
    ax.FontSize = 20;
    ax.XLabel.FontSize = 45;
    ax.YLabel.FontSize = 45;
    legend({'Oil', 'Solvent'}, 'fontSize', 45);
    
    drawnow;
    
    ax.Units = 'pixels';
    
    
    m = 0; n = 0;
    rect = [m, n, 2000-2*m, 1000-2*n];
    
    M(i) = getframe(fig, rect);
    
    ax.Units = 'normalized';
    
end

%%

% n = numel(states);
n = 400;
% n = 4;
pth = [mrstPath('mrst-solvent'), '/presentation/figures/displacement1D/'];
vo = VideoWriter([pth, 'displacement1D_1.avi']);
vo.FrameRate = 1000/30;
open(vo);

writeVideo(vo, M);

close(vo)

%%

clc

pv = poreVolume(G, rock);
n = numel(ws);
nw = numel(W_G);
wellW = zeros(n,nw);
wellO = zeros(n,nw);
wellS = zeros(n,nw);
wellM = cell(3,1);
resM = zeros(n,4);

wc = vertcat(W_G.cells);
for i = 1:n    
    
%     rhoW = [fluid.rhoWS, (states{i}.rho(wc(2:end),1))'];
%     rhoW = (states{i}.rho(wc,1))';
%     wellW(i,:) = [ws{i}.qWs].*schedule.step.val(i).*fluid.rhoWS;
    wellM{1}(i,:) = [ws{i}.qWs].*schedule.step.val(i).*fluid.rhoWS;
    
%     rhoO = [fluid.rhoOS, (states{i}.rho(wc(2:end),2))'];
%     rhoO = (states{i}.rho(wc,2))';
%     wellO(i,:) = [ws{i}.qOs].*schedule.step.val(i).*fluid.rhoOS;
    wellM{2}(i,:) = [ws{i}.qOs].*schedule.step.val(i).*fluid.rhoOS;
    
%     rhoO = [fluid.rhoOS, (states{i}.rho(wc(2:end),2))'];
%     rhoO = (states{i}.rho(wc,2))';
%     wellO(i,:) = [ws{i}.qOs].*schedule.step.val(i).*fluid.rhoOS;
    wellM{3}(i,:) = [ws{i}.qGs].*schedule.step.val(i).*fluid.rhoGS;
    
    rhoS = [fluid.rhoSS, (states{i}.rho(wc(2:end),4))'];
%     rhoS = (states{i}.rho(wc,4))';
%     wellS(i,:) = [ws{i}.qSs].*schedule.step.val(i).*fluid.rhoSS;
    wellM{4}(i,:) = [ws{i}.qSs].*schedule.step.val(i).*fluid.rhoSS;
    
    for phNo = 1:4
        resM(i,phNo) = sum(states{i}.s(:,phNo).*states{i}.rho(:,phNo).*pv);
    end
    
end

wellMtot = zeros(n,nw);
for phNo = 1:4
    wellMtot = wellMtot + wellM{phNo};
end
resMtot = sum(resM,2);


wellMcum = cellfun(@(m) cumsum(m,1), wellM, 'uniformOutput', false);
wellMtotCum = cumsum(wellMtot,1);

errTotRel = (resMtot(1) + sum(wellMtotCum,2) - resMtot)./resMtot;

errTot = abs(resMtot(1) + sum(wellMtotCum(end,:)) - resMtot(end));

err = cell(4,1);
for phNo = 1:4
    err{phNo} = (resM(1, phNo) + sum(wellM{phNo},2) - resM(:,phNo))./resMtot;
end


fprintf(['Absolute error: \t %.2d \n'], ...
         errTot);
         
%%

W_G(1).name = 'I';
W_G(2).name = 'P';

figure(1); clf
plotGrid(G, 'edgecolor', 'none'); axis equal off
plotWell(G, W_G, 'height', 0.5)
w = [-17,14];
view([-20,10])
light('Position', [2,-3,0])

