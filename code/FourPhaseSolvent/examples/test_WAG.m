mrstModule add ad-core ad-eor ad-blackoil ad-props blackoil-sequential matlab_bgl

gravity reset on

n = 100;
G = computeGeometry(cartGrid([n,1,1]));
rock = makeRock(G, 100*milli*darcy, 1);

T = 1*year;
pv = poreVolume(G, rock);
injRate = 2*sum(pv)/T;
nStep = 100;
dT = T/nStep;

%%

fluid = initSimpleADIFluid('n'     , [2, 2, 2], ...
                           'rho'   , [1000, 800, 100]*kilogram/meter^3, ...
                           'phases', 'WOG', ...
                           'mu'    , [1, 10, 2]*centi*poise);

fluid = addSolventProperties(fluid, 'n', 2, ...
                                    'rho', 100*kilogram/meter^3, ...
                                    'mixPar', 0, ...
                                    'mu'    , 1*centi*poise, ...
                                    'sOres_i', 0.3, ...
                                    'sOres_m', 0.1, ...
                                    'Msat', @(sG, sS) sS.*0);
                                
model = FourPhaseSolventModel(G, rock, fluid);
model.extraStateOutput = true;

%%

WA = addWell([], G, rock, 1, ...
                 'comp_i', [0, 0, 0, 1], ...
                 'type'  , 'rate', ...
                 'val'   , injRate);
     
WA = addWell(WA, G, rock, G.cells.num, ...
                 'comp_i', [0, 0, 0, 1], ...
                 'type', 'bhp', ...
                 'val', 0      );

WB = addWell([], G, rock, 1, ...
                 'comp_i', [1, 0, 0, 0], ...
                 'type'  , 'rate', ...
                 'val'   , injRate);
     
WB = addWell(WB, G, rock, G.cells.num, ...
                 'comp_i', [0, 0, 0, 1], ...
                 'type', 'bhp', ...
                 'val', 0      );
             
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

[wsS4, statesS4, reportsS4] = simulateScheduleAD(state0, model, schedule);


%%

mrstModule add mrst-gui

figure(1); clf
plotToolbar(G, statesS4, 'plot1d', true)
ylim([0 1]);
% 
% figure(2); clf
% plotToolbar(G, statesW, 'plot1d', true)
% ylim([0 1]);


%%

n = numel(wsS4);
qs = zeros(n,2);
for i = 1:n
    qs(i,:) = [wsS4{i}.qWs].*dT(i).*fluid.rhoWS;
end
qsc = cumsum(qs,1);
qr = sum(statesS4{end}.s(:,1).*statesS4{end}.rho(:,1).*pv);

%%

n = numel(wsS4);
qs = zeros(n,2);
for i = 1:n
    qs(i,:) = [wsS4{i}.qSs].*step.val(i).*fluid.rhoSS;
end
qsc = cumsum(qs,1);
qr = sum(statesS4{end}.s(:,4).*statesS4{end}.rho(:,4).*pv);

%%

n = numel(wsS4);

[ql, qr] = deal(0);

for i = 1:n

    ql = ql + statesS4{i}.flux(2,1).*statesS4{i}.rho(1,1).*step.val(i);
    qr = qr + statesS4{i}.flux(200,1).*statesS4{i}.rho(G.cells.num-1,1).*step.val(i);
    
end

res = sum(pv(2:G.cells.num-1).*statesS4{end}.s(2:G.cells.num-1,1).*statesS4{end}.rho(2:G.cells.num-1,1));

