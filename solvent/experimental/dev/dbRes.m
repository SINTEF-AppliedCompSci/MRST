mrstModule add ad-core ad-eor ad-blackoil ad-props sequential matlab_bgl

gravity reset on

n = 100;
% G = computeGeometry(cartGrid([n,1,1], [100,1,1]));
G = computeGeometry(cartGrid([1,1,n], [1,1,100]));
rock = makeRock(G, 100*milli*darcy, 1);

%%

fluid = initSimpleADIFluid('n'     , [2, 2, 2], ...
                           'rho'   , [1000, 800, 100]*kilogram/meter^3, ...
                           'phases', 'WOG', ...
                           'mu'    , [1, 10, 2]*centi*poise, ...
                           'c'     , [1e-6, 1e-6, 1e-5]/barsa);

sOres_i= 0.3;
fluid = addSolventProperties(fluid, 'n', 2, ...
                                    'rho', 100*kilogram/meter^3, ...
                                    'mixPar', 2/3, ...
                                    'mu'    , 2*centi*poise, ...
                                    'sOres_i', sOres_i, ...
                                    'sOres_m', 0.0, ...
                                    'c'      , 1e-5/barsa);
                                
model = FourPhaseSolventModel(G, rock, fluid);
% model.extraStateOutput = true;

T = 1*year;
nStep = 100;
rate = sum(poreVolume(G, rock))/T;

W = addWell([], G, rock, 1, 'type', 'rate', 'val', rate, 'comp_i', [0,0,0,1]);
W = addWell(W, G, rock, G.cells.num, 'type', 'bhp', 'val', 0, 'comp_i', [1,0,0,0]);
% 
dT = T/nStep;
dT = rampupTimesteps(T, dT, 0);
step.val = dT;
step.control = ones(numel(dT),1);
schedule.step = step;
schedule.control(1).W = W;

sO = sOres_i;
sW = 1 - sOres_i;

state0 = initResSol(G, 100*barsa, [sW, sO, 0, 0]);
state0.wellSol = initWellSolAD(W, model, state0);

nls = NonLinearSolver('useLineSearch', true);
nls = NonLinearSolver('useLineSearch', false);

%%

[ws, states, reports] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nls);

%%

plotToolbar(G, states, 'plot1d', true);


%%

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
    err{phNo} = (resM(1, phNo) + sum(wellMcum{phNo},2) - resM(:,phNo))./resMtot;
end


fprintf(['Absolute error: \t %.2d \n'], ...
         errTot);

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
