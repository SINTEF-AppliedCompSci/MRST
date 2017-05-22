mrstModule add ad-core ad-eor ad-blackoil ad-props blackoil-sequential matlab_bgl

gravity reset on

n = 100;
G = computeGeometry(cartGrid([n,1,1]));
rock = makeRock(G, 100*milli*darcy, 1);

%%

fluid = initSimpleADIFluid('n'     , [2, 2, 2], ...
                           'rho'   , [1000, 800, 100]*kilogram/meter^3, ...
                           'phases', 'WOG', ...
                           'mu'    , [1, 10, 2]*centi*poise);

fluid = addSolventProperties(fluid, 'n', 2, ...
                                    'rho', 100*kilogram/meter^3, ...
                                    'mixPar', 2/3, ...
                                    'mu'    , 1*centi*poise, ...
                                    'sOres_i', 0.3, ...
                                    'sOres_m', 0.1);
                                
model = FourPhaseSolventModel(G, rock, fluid);
model.extraStateOutput = true;

[schedule, W_G, W_W] = makeWAGschedule(model, {1}, {G.cells.num}, 4, 'T', 4*year, 'nStep', 1000);
state0 = initResSol(G, 100*barsa, [0 1 0 0]);
state0.wellSol = initWellSolAD(W_G, model, state0);

%%

[ws, states, reports] = simulateScheduleAD(state0, model, schedule);

%%

mrstModule add mrst-gui

figure(1); clf
plotToolbar(G, states, 'plot1d', true)
ylim([0 1]);

%%

n = numel(wsS4);
qs = zeros(n,2);
for i = 1:n
    qs(i,:) = [wsS4{i}.qWs].*dT(i).*fluid.rhoWS;
end
qsc = cumsum(qs,1);
qr = sum(statesS4{end}.s(:,1).*statesS4{end}.rho(:,1).*pv);

%%

pv = poreVolume(G, rock);

n = numel(ws);
qs = zeros(n,2);
for i = 1:n
    qs(i,:) = [ws{i}.qSs].*schedule.step.val(i).*fluid.rhoSS;
end
qsc = cumsum(qs,1);
qr = sum(states{end}.s(:,4).*states{end}.rho(:,4).*pv);

%%
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
         

