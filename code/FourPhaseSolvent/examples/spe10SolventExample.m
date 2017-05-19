mrstModule add mrst-solvent spe10

%% Use setupSPE10_AD to Fetch the SPE10 model
% We pick up only one layer 
%
layers = 35;
[~, model, ~] = setupSPE10_AD('layers', layers);
% We recover the grid and rock properties from the model
G = model.G;
rock = model.rock;

%%

gravity reset on

fluid = initSimpleADIFluid('n'     , [2, 2, 2], ...
                           'rho'   , [1000, 800, 100]*kilogram/meter^3, ...
                           'phases', 'WOG', ...
                           'mu'    , [1, 10, 2]*centi*poise);
%'c'     , [1e-7, 1e-6, 1e-6]/barsa, ...
                       
fluid = addSolventProperties(fluid, 'n', 2, ...
                                    'rho', 100*kilogram/meter^3, ...
                                    'mixPar', 2/3, ...
                                    'mu'    , 1*centi*poise);
%'c', 1e-6/barsa, ...                                

model = FourPhaseSolventModel(G, rock, fluid);
model.extraStateOutput = true;

T = 1*year;
pv = poreVolume(G, rock);
injRate = 1*sum(pv)/T;

nStep = 100;
dT = T/nStep;
%%

WA = verticalWell([], G, rock, 1, 1, [], ...
                 'comp_i', [0, 0, 0, 1], ...
                 'type'  , 'rate', ...
                 'val'   , injRate);
     
WA = verticalWell(WA, G, rock, 60, 220, [], ...
                 'comp_i', [0, 0, 0, 1], ...
                 'type', 'bhp', ...
                 'val', 0      );

WB = verticalWell([], G, rock, 1, 1, [], ...
                 'comp_i', [1, 0, 0, 0], ...
                 'type'  , 'rate', ...
                 'val'   , injRate);

WB = verticalWell(WB, G, rock, 60, 220, [], ...
                 'comp_i', [0, 0, 0, 1], ...
                 'type'  , 'bhp', ...
                 'val'   , 0);
             
control(1).W = WA;
control(2).W = WB;

injStart = 0;
injStop = 0.25;

dT_A = rampupTimesteps(injStop*T, dT);
dT_B = rampupTimesteps((1 - injStop)*T, dT);

step.val = [dT_A; dT_B];
step.control = [1*ones(numel(dT_A),1); 2*ones(numel(dT_B),1)];
schedule.control = control;
schedule.step = step;

state0 = initResSol(G, 100*barsa, [0 1 0 0]);
state0.wellSol = initWellSolAD(WA, model, state0);

%%

[ws, states, reports] = simulateScheduleAD(state0, model, schedule);


%%

plotToolbar(G, states)

%%

close all
n = numel(states);
pos = [500, 500, 2000,1000];
fig = figure('position', pos);
M = struct('cdata',[],'colormap',[]);
jv = [4,1,2];
ttl = {'S_s', 'S_w', 'S_o'};
for i = 1:n
    
    clf;
    
    subplot(1,3,1)
    plotCellData(G, states{i}.s(:,4), 'edgecolor', 'none');
    colorbar;
    colormap(jet)
    caxis([0,1])
    axis equal tight off
    title('S_s');
    
    subplot(1,3,2)
    plotCellData(G, states{i}.mob(:,2), 'edgecolor', 'none');
    colorbar;
    caxis([0,130])
    colormap(jet)
    axis equal tight off
    title('\lambda_o');
    
    subplot(1,3,3)
    plotCellData(G, states{i}.s(:,2), 'edgecolor', 'none');
    colorbar;
    colormap(jet)
    caxis([0,1])
    axis equal tight off
    title('S_o');
     
    drawnow;
    ax = gca;
    ax.Units = 'pixels';
%     pos = ax.Position;
    
    m = 150; n = 100;
    rect = [m, n, 2000-2*m, 1000-2*n];
    
    M(i) = getframe(fig, rect);
    
    ax.Units = 'normalized';
    
end

%%

vo = VideoWriter('solvent.avi');
vo.FrameRate = n/15;
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

qs = 0;
for i = 1:n
    qs = qs + [ws{i}.qSs]*step.val(i);
end

%%
n = numel(ws);
qs = zeros(n,2);
for i = 1:n
    qs(i,:) = [ws{i}.qSs].*step.val(i).*fluid.rhoSS;
end
qsc = cumsum(qs,1);
qr = sum(states{end}.s(:,4).*states{end}.rho(:,4).*pv);
