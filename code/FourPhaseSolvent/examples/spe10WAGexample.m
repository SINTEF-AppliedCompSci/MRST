%% Add necessary modules

mrstModule add ad-blackoil ad-props % AD-moduels
mrstModule add spe10                % SPE10 dataset
mrstModule add solvent              % Solvent model

gravity reset on

df = get(0, 'defaultfigureposition');
close all

%% Set up grid and rock properties

% We pick up layer 10 of SPE10, and extract the grid and rock properties.
% We will define our own fluid.
layers = 10;
[~, model, ~] = setupSPE10_AD('layers', layers);
G    = model.G;
rock = model.rock;

%% Define fluid properties

% We start by defining a three-phase fluid with water, oil and gas
fluid = initSimpleADIFluid('n'     , [2, 2, 2], ...
                           'rho'   , [1000, 800, 100]*kilogram/meter^3, ...
                           'phases', 'WOG', ...
                           'mu'    , [1, 3, 0.1]*centi*poise, ...
                           'c'     , [1e-6, 1e-5, 1e-4]/barsa);
                       
% The solvent model we will use, treats solvent gas as a fourth phase,
% which is either miscible or immiscible with the oil and gas, depending on
% the fraction of solvent concentration to total gas concentration,
% $S_s/(S_g + S_s)$, and the pressure.

sOres_i = 0.34; % Immiscible residual oil saturation
sOres_m = 0.12; % Miscible residual oil saturation
fluid   = addSolventProperties(fluid, 'n'      , 2                   , ...
                                    'rho'    , 100*kilogram/meter^3, ...
                                    'mixPar' , 2/3                 , ...
                                    'mu'     , 0.2*centi*poise     , ...
                                    'sOres_i', sOres_i             , ...
                                    'sOres_m', sOres_m             , ...
                                    'c'      , 1e-4/barsa          );

%% Inspect fluid

figure('Position', [df(1:2), 1000, 400]);

n = 100;
[sO,sS] = meshgrid(linspace(0,1,n)', linspace(0,1,n)');
sW      = 0;
sG      = 1-(sW+sO+sS);

ss  = deal(linspace(0,1-sOres_m,100)');
b   = sOres_i + 2*ss - 1;
sg  = (-b + sqrt(b.^2 - 4*(ss.*(sOres_m + ss - 1))))/2;
sOr = 1 - sg - ss;

sW = sW(:); sO = sO(:); sG = sG(:); sS = sS(:);
[sWres, sOres , sSGres ]  = computeResidualSaturations(fluid, 0 , sG , sS );
[krW, krO, krG, krS]      = computeRelPermSolvent(fluid, 0, sW, sO, sG, sS, sWres, sOres, sSGres, mobMult);
sO = reshape(sO, n,n); sG = reshape(sG, n,n); sS = reshape(sS, n,n);

phName = {'O', 'G', 'S'};

for phNo = 1:3

    subplot(1,3,phNo)
    
    relpermName = ['kr', phName{phNo}];
    kr = reshape(eval(relpermName),[n,n]);
    kr(sG<0) = nan;

    contourf(mapx(sG, sS, sO), mapy(sG,sS,sO), kr, 20, 'linecolor', 0.5.*[1,1,1])
    [mapx, mapy] = ternaryAxis('names', {'S_g', 'S_s', 'S_o'});
    plot(mapx(sg, ss, sOr), mapy(sg,ss,sOr), 'color', 0.99*[1,1,1], 'linewidth', 2)
    
    axis([0,1,0,sqrt(1-0.5^2)]); axis equal
    title(relpermName, 'position', [0.5,-0.2]); 
    
end
                                
%%
                                
model = FourPhaseSolventModel(G, rock, fluid);
model.extraStateOutput = true;

producers = {1, G.cartDims(1), G.cartDims(1)*G.cartDims(2) - G.cartDims(1) + 1, G.cartDims(1)*G.cartDims(2)};
injectors = {round(G.cartDims(1)*G.cartDims(2)/2 + G.cartDims(1)/2)};

time = 2*year;
pv = poreVolume(G, rock);
rate = sum(pv)/time;
bhp = 50*barsa;

W = [];
for pNo = 1:numel(producers)
    W = addWell(W, G, rock, producers{pNo}, ...
                'type', 'bhp', ...
                'val', bhp, ...
                'comp_i', [1,0,0,0], ...
                'sign', -1, ...
                'name', ['P', num2str(pNo)]);
end
for iNo = 1:numel(injectors)
    W = addWell(W, G, rock, injectors{iNo}, ...
                'type', 'rate', ...
                'val', rate/numel(injectors), ...
                'comp_i', [1,0,0,0], ...
                'sign', 1, ...
                'name', ['I', num2str(iNo)]);
end

%%

dt          = 30*day;
useRampUp   = true;
scheduleWAG = makeWAGschedule(W, 4, 'time'     , time     , ...
                                    'dt'       , dt       , ...
                                    'useRampup', useRampUp);
schedule = scheduleWAG;

tvec                  = rampupTimesteps(time, dt);
control               = 2*ones(numel(tvec),1);
schedule.step.val     = [tvec; scheduleWAG.step.val];
schedule.step.control = [control; scheduleWAG.step.control];

%%
sO = 0.6; sG = 0.2;
state0 = initResSol(G, 100*barsa, [1-sO-sG, sO sG, 0]);
state0.wellSol = initWellSolAD(W, model, state0);

%%

[ws, states, reports] = simulateScheduleAD(state0, model, schedule);

%%

plotToolbar(G, states)

%%

W = addWell([], G, rock, injectors{1}, 'type', 'rate', 'val', rate, 'comp_i', [1,0,0,0]);
for i = 1:numel(producers)
    W = addWell(W, G, rock, producers{i}, 'type', 'bhp', 'val', 0, 'comp_i', [1,0,0,0]);
end
dT = repmat(time/nStep, nStep,1);
step.val = dT;
control(1).W = W;
step.control = ones(nStep,1);
schedule.control = control;
schedule.step = step;

state0 = initResSol(G, 100*barsa, [0 1 0 0]);
state0.wellSol = initWellSolAD(W_G, model, state0);

%%

[wsW, statesW, reportsW] = simulateScheduleAD(state0, model, schedule);

%%

figure(1)
plotToolbar(G, states)
axis equal tight
colorbar
% 
% figure(2)
% plotToolbar(G, statesW)
% axis equal tight
% colorbar

%%

close all
ns = numel(states);
pos = [500, 500, 2000,1000];
fig = figure('position', pos);
M = struct('cdata',[],'colormap',[]);
jv = [4,1,2];
% ns = 4;
for i = 1:ns
    
    clf;
    
    subplot(1,3,1)
    plotCellData(G, states{i}.s(:,1), 'edgecolor', 'none');
    colorbar('fontsize', 20);
    colormap(jet)
    caxis([0,1])
    axis equal tight off
    title('S_w');
    ax = gca;
    ax.FontSize = 20;
    
    subplot(1,3,2)
    plotCellData(G, states{i}.s(:,4), 'edgecolor', 'none');
    colorbar('fontSize', 20);
    caxis([0,1])
    colormap(jet)
    axis equal tight off
    title('S_s');
    ax = gca;
    ax.FontSize = 20;
    
    subplot(1,3,3)
    plotCellData(G, states{i}.s(:,2), 'edgecolor', 'none');
    colorbar('fontSize', 20);
    colormap(jet)
    caxis([0,1])
    axis equal tight off
    title('S_o');
    ax = gca;
    ax.FontSize = 20;
    
    drawnow;
    ax = gca;
    ax.Units = 'pixels';
%     pos = ax.Position;
    
    m = 150; n = 150;
    rect = [m, n, 2000-m-20, 1000-n-40];
    
    M(i) = getframe(fig, rect);
    
    ax.Units = 'normalized';
    
end

%%

pth = [mrstPath('mrst-solvent'), '/presentation/figures/spe10/spe10_1'];
vo = VideoWriter(pth);
n = numel(states);
vo.FrameRate = n/30;
open(vo);

writeVideo(vo, M);

close(vo)

%%

close all
pos = [500, 500, 2000,1000];
fig = figure('position', pos);
axis off

movie(M)

%%

wellS = 0;
for i = 1:n
    wellS = wellS + [ws{i}.qSs]*step.val(i);
end

% %%
% 
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

err = (resMtot(1) + sum(wellMtotCum,2) - resMtot)./resMtot;

errTot = abs(resMtot(1) + sum(wellMtotCum(end,:)) - resMtot(end));

fprintf(['Absolute error: \t %.2d \n'], ...
         errTot);
     
%%

path = [mrstPath('mrst-solvent'), '/presentation/figures/'];
savepng = @(name) print([path, 'spe10/', name], '-dpng', '-r300');
% savepng = @(name) [];

close all

figure(1); clf;
perm = rock.perm(:,1);
plotCellData(G, log10(perm), 'edgecolor', 'none');
axis equal tight off
% view([90,90])
ax = gca;
pos = ax.Position;
logColorbar%('southoutside', 'position', [pos(1), pos(2)+0.07, pos(3), 0.05*pos(4)]);

savepng('spe10perm');

figure(2); clf;
plotCellData(G, rock.poro, 'edgecolor', 'none');
axis equal tight off
% view([90,90])
ax = gca;
pos = ax.Position;
colorbar%('southoutside', 'position', [pos(1), pos(2)+0.07, pos(3), 0.05*pos(4)]);

savepng('spe10poro');
