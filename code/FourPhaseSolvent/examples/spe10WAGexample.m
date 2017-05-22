mrstModule add mrst-solvent spe10

%% Use setupSPE10_AD to Fetch the SPE10 model
% We pick up only one layer 
%
layers = 10;
[~, model, ~] = setupSPE10_AD('layers', layers);
G = model.G;
rock = model.rock;

%%

gravity reset on

fluid = initSimpleADIFluid('n'     , [2, 2, 2], ...
                           'rho'   , [1000, 800, 100]*kilogram/meter^3, ...
                           'phases', 'WOG', ...
                           'mu'    , [1, 10, 2]*centi*poise);
                       
fluid = addSolventProperties(fluid, 'n', 2, ...
                                    'rho', 100*kilogram/meter^3, ...
                                    'mixPar', 2/3, ...
                                    'mu'    , 1*centi*poise, ...
                                    'sOres_i', 0.3, ...
                                    'sOres_m', 0.05);
                                
model = FourPhaseSolventModel(G, rock, fluid);
model.extraStateOutput = true;

% injectors = {round(G.cartDims(1)*G.cartDims(2)/2 + G.cartDims(1)/2)};
% producers = {1, G.cartDims(1), G.cartDims(1)*G.cartDims(2) - G.cartDims(1) + 1, G.cartDims(1)*G.cartDims(2)};

injectors = {1};
producers = {G.cells.num};

T = 1*year;
pv = sum(poreVolume(G, rock));
rate = 2*pv/T;
nStep = 100;

[schedule, W_G, W_W] = makeWAGschedule(model, injectors, producers, 4, 'T', T, 'gRate', rate, 'wRate', rate, 'nStep', nStep);
state0 = initResSol(G, 100*barsa, [0 1 0 0]);
state0.wellSol = initWellSolAD(W_G, model, state0);

%%

[ws, states, reports] = simulateScheduleAD(state0, model, schedule);

%%

W = addWell([], G, rock, injectors{1}, 'type', 'rate', 'val', rate, 'comp_i', [1,0,0,0]);
for i = 1:numel(producers)
    W = addWell(W, G, rock, producers{i}, 'type', 'bhp', 'val', 0, 'comp_i', [1,0,0,0]);
end
dT = repmat(T/nStep, nStep,1);
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

figure(2)
plotToolbar(G, statesW)
axis equal tight
colorbar

%%

% close all
% n = numel(states);
% pos = [500, 500, 2000,1000];
% fig = figure('position', pos);
% M = struct('cdata',[],'colormap',[]);
% jv = [4,1,2];
% ttl = {'S_s', 'S_w', 'S_o'};
% for i = 1:n
%     
%     clf;
%     
%     subplot(1,3,1)
%     plotCellData(G, states{i}.s(:,4), 'edgecolor', 'none');
%     colorbar;
%     colormap(jet)
%     caxis([0,1])
%     axis equal tight off
%     title('S_s');
%     
%     subplot(1,3,2)
%     plotCellData(G, states{i}.mob(:,2), 'edgecolor', 'none');
%     colorbar;
%     caxis([0,130])
%     colormap(jet)
%     axis equal tight off
%     title('\lambda_o');
%     
%     subplot(1,3,3)
%     plotCellData(G, states{i}.s(:,2), 'edgecolor', 'none');
%     colorbar;
%     colormap(jet)
%     caxis([0,1])
%     axis equal tight off
%     title('S_o');
%      
%     drawnow;
%     ax = gca;
%     ax.Units = 'pixels';
% %     pos = ax.Position;
%     
%     m = 150; n = 100;
%     rect = [m, n, 2000-2*m, 1000-2*n];
%     
%     M(i) = getframe(fig, rect);
%     
%     ax.Units = 'normalized';
%     
% end
% 
% %%
% 
% vo = VideoWriter('solvent.avi');
% vo.FrameRate = n/15;
% open(vo);
% 
% writeVideo(vo, M);
% 
% close(vo)
% 
% %%
% 
% close all
% pos = [500, 500, 2000,1000];
% fig = figure('position', pos);
% axis off
% 
% movie(M)
% 
% %%

wellS = 0;
for i = 1:n
    wellS = wellS + [ws{i}.qSs]*step.val(i);
end

% %%
% 
% clc
% 
% pv = poreVolume(G, rock);
% n = numel(ws);
% wellW = zeros(n,5);
% wellO = zeros(n,5);
% wellS = zeros(n,5);
% 
% wc = vertcat(W_G.cells);
% for i = 1:n    
%     
%     
%     rhoW = [fluid.rhoWS, (states{i}.rho(wc(2:end),1))'];
% %     rhoW = (states{i}.rho(wc,1))';
%     wellW(i,:) = [ws{i}.qWs].*schedule.step.val(i).*fluid.rhoWS;
%     
%     rhoO = [fluid.rhoOS, (states{i}.rho(wc(2:end),2))'];
% %     rhoO = (states{i}.rho(wc,2))';
%     wellO(i,:) = [ws{i}.qOs].*schedule.step.val(i).*fluid.rhoOS;
%     
%     rhoS = [fluid.rhoSS, (states{i}.rho(wc(2:end),4))'];
% %     rhoS = (states{i}.rho(wc,4))';
%     wellS(i,:) = [ws{i}.qSs].*schedule.step.val(i).*fluid.rhoSS;
%     
% end
% wellWT = cumsum(wellW,1);
% wellOT = cumsum(wellO,1);
% wellST = cumsum(wellS,1);
% 
% resWs = 0;
% resOs = sum(pv.*fluid.rhoOS);
% resSs = 0;
% 
% % resWs = sum(states{1}.s(:,1).*states{1}.rho(:,1).*pv);
% % resOs = sum(states{1}.s(:,2).*states{1}.rho(:,2).*pv);
% % resSs = sum(states{1}.s(:,4).*states{1}.rho(:,4).*pv);
% 
% resWe = sum(states{end}.s(:,1).*states{end}.rho(:,1).*pv);
% resOe = sum(states{end}.s(:,2).*states{end}.rho(:,2).*pv);
% resSe = sum(states{end}.s(:,4).*states{end}.rho(:,4).*pv);
% 
% resWs + sum(wellWT(end,:)) - resWe
% resOs + sum(wellOT(end,:)) - resOe
% resSs + sum(wellST(end,:)) - resSe
% 
