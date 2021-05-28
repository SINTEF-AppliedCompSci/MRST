mrstModule add deckformat ad-props ad-blackoil ad-core

%%

gravity reset on

% filename = 'SIMPLE1D';
filename = 'SPE1';

datafldr = fullfile(mrstPath('solvent'), 'code', 'FourPhaseSolvent', 'examples', 'dev', 'data');
deckpath = fullfile(datafldr, filename);
deck     = readEclipseDeck([deckpath, '.DATA']);
deck     = convertDeckUnits(deck);

G = initEclipseGrid(deck);
G = computeGeometry(G);

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

fluid = initDeckADIFluid(deck);
fluid = addSolventProperties(fluid, 'overwrite', false, 'smin', 1e-16);

switch filename
    case 'SIMPLE1D'
        state0 = initResSol(G, 50*barsa, [0.2, 0.6, 0.2, 0]);
%         state0.rs = 200*ones(G.cells.num,1);
        disgas = false;
%         disgas = true;
        unts = 'metric';
        
    case 'SPE1'
        [~, ~, ~, ~, state0] = setupSPE1();
        state0.s(:,4) = zeros(G.cells.num,1);
        disgas = true;
        unts = 'field';
end

dynEPS = true;
hysSat = false;
model  = FourPhaseSolventModel(G, rock, fluid, 'disgas', disgas, 'dynamicEndPointScaling', dynEPS, 'extraStateOutput', true, 'hystereticResSat', hysSat);

schedule = convertDeckScheduleToMRST(model, deck);
dt       = rampupTimesteps(sum(schedule.step.val),    30*day, 10);
schedule.step.val = dt;
schedule.step.control = ones(numel(dt),1);

%%

model.dpMaxRel = 0.05;
model.dsMaxAbs = 0.05;

[wellSols, states, reports] = simulateScheduleAD(state0, model, schedule);

%%

% close all

time = cumsum(schedule.step.val);
[wellSolsOPM, timeOPM] = convertSummaryToWellSolsSolvent(deckpath, unts);


plotWellSols({wellSols, wellSolsOPM}, {time, timeOPM});
% plotWellSols({wellSols}, {time});

%%
% 
% compd = 1:(size(smry.data, 2));
% 

% 

% 
% convertSummaryToWellSols(
% 
% names = {'INJ', 'PROD'};
% 
% % propsECL  = {'WBHP', 'WOPR', 'WGPR', 'WNPR', 'WNIR'       };
% % propsMRST = {'bhp' , 'qOs' , 'qGs' , 'qSs' , 'qSs'     };
% propsECL  = {'WBHP', 'WOPR'};
% propsMRST = {'bhp' , 'qOS'};
% % unts      = [psia  , 1000*ft^3/day];
% % unts      = [barsa, stb/day, 1000*ft^3/day, meter^3/day, meter^3/day];
% 
% unts      = [psia];%, stb/day, 1000*ft^3/day, meter^3/day, meter^3/day];
% post      = {@(p) p};%, @(p) abs(p), @(p) abs(p), @(p) abs(p), @(p) abs(p)};
% 
% for pNo = 1:numel(propsECL)
%    
%     figure('name', propsECL{pNo})
%     
%     for wNo = 1:numel(names)
%    
%         subplot(1,2,wNo);
%         hold on
%         prpMRST = cellfun(@(ws) ws(wNo).(propsMRST{pNo}), wellSols);
%         prpMRST = post{pNo}(prpMRST);
%         prpECL = convertFrom(smry.get(names(wNo), propsECL{pNo}, compd), unts(pNo))';
%         plot(time, prpMRST)
%         plot(timeComp, prpECL, 'o')
% %         plot(prpMRST)
% %         plot(prpECL, 'o')
%         hold off
%         
%     end
%     
% end

%%

[statesOPM, restart] = convertRestartToStates(deckpath  , G    , ...
                                'use_opm'            , true, ...
                                'includeWellSols'    , true , ...
                                'wellSolsFromRestart', true , ...
                                'includefluxes'      , true , ...
                                'consistentWellSols' , false , ...
                                'includeMobilities'  , true );

% 
% %%
% 
% sW = 0; sO = 0.2; sG = 0.4; sS = 1-(sW + sO + sG);
% p = 9000*psia; mobMult = 1;
% [krW, krO, krG, krS] = computeRelPermSolvent(model, p, sW, sO, sG, sS, fluid.sWr, fluid.sOr_i, fluid.sGc_i, mobMult);
% 
%%

%%

close all

smry = readEclipseSummaryUnFmt(deckpath);

compd = 1:(size(smry.data, 2));
filename = 'PROD';

timeComp = smry.get(':+:+:+:+', 'TIME', compd)*day';

time = cumsum(schedule.step.val);

names = {'INJ', 'PROD'};

% propsECL  = {'WBHP', 'WOPR', 'WGPR', 'WNPR', 'WNIR'       };
% propsMRST = {'bhp' , 'qOs' , 'qGs' , 'qSs' , 'qSs'     };
propsECL  = {'WBHP'};
propsMRST = {'bhp' };
% unts      = [psia  , 1000*ft^3/day];
% unts      = [barsa, stb/day, 1000*ft^3/day, meter^3/day, meter^3/day];

unts      = [psia];%, stb/day, 1000*ft^3/day, meter^3/day, meter^3/day];
post      = {@(p) p};%, @(p) abs(p), @(p) abs(p), @(p) abs(p), @(p) abs(p)};

for pNo = 1:numel(propsECL)
   
    figure('name', propsECL{pNo})
    
    for wNo = 1:numel(names)
   
        subplot(1,2,wNo);
        hold on
        prpMRST = cellfun(@(ws) ws(wNo).(propsMRST{pNo}), wellSols);
        prpMRST = post{pNo}(prpMRST);
        prpECL = convertFrom(smry.get(names(wNo), propsECL{pNo}, compd), unts(pNo))';
        plot(time, prpMRST)
        plot(timeComp, prpECL, 'o')
%         plot(prpMRST)
%         plot(prpECL, 'o')
        hold off
        
    end
    
end

%%

save('bo.mat', 'eqs', 'rho', 'mob', 'sat', 'rs', 'rv', 'bW', 'bO', 'bG')

%%

save('s.mat', 'eqs', 'rho', 'mob', 'sat', 'rs', 'rv')

%%
% 
% dt = rampupTimesteps(sum(schedule.step.val),    30*day, 10);
% schedule.step.val = dt;
% schedule.step.control= ones(numel(dt),1);
% 
% model.dpMaxRel = 0.05;
% model.dsMaxAbs = 0.1;
% model.maximumPressure = 10e3*barsa;
% model.minimumPressure = 0*barsa;
% nls = NonLinearSolver('maxTimestepCuts', 15);
% [wellSols, states, reports] = simulateScheduleAD(state0, model, schedule, 'nonlinearsolver', nls);
% 
% %%
% 
% 
% smry = readEclipseSummaryUnFmt(deckpath);
% 
% compd = 1:(size(smry.data, 2));
% name = 'PROD';
% bhpComp = convertFrom(smry.get(name, 'WBHP', compd), psia)';
% timeComp = smry.get(':+:+:+:+', 'TIME', compd)';
% 
% %%
% 
% [G, rock] = setupSPE1();
% rock.perm(:,3) = rock.perm(:,1);
% model = FourPhaseSolventModel(G, rock, fluid);
% 
% W = schedule.control(1).W;
% for wNo = 1:numel(W)
%     if W(wNo).sign > 0
%         W(wNo).compi = [0,0,0,1];
%     end
% end
% schedule.control(1).W = W;
% 
% %%
% 
% 
% 
% [ws, states, reports] = simulateScheduleAD(state0, model, schedule);
% 
% % W = [];
% % W = verticalWell(W, G, rock, 1, 1, [], 'type', 'rate');
% % W = verticalWell(W, G, rock, 10, 10, []);
% % 


schedule.control.W(1).compi = [1,0,0];
schedule.control.W(2).compi = [1,0,0];
schedule.control.W(1).val = 1e-4;

% state0.s = state0.s(:,1:3);

modelBO = ThreePhaseBlackOilModel(G, rock, fluid, 'disgas', true);
modelBO.dpMaxRel = 0.05;
modelBO.dsMaxAbs = 0.1;

[wellSolsBO, statesBO, reports] = simulateScheduleAD(state0, modelBO, schedule);

%%

% schedule.control.W(1).compi = [1,0,0,0];
% schedule.control.W(2).compi = [1,0,0,0];
% schedule.control.W(1).val = 1e-4;

state0 = initResSol(G, 50*barsa, [0.2, 0.6, 0.2, 0]);

% state0.s(:,4) = 1 - sum(state0.s(:,1:3),2);


model.dpMaxRel = 0.05;
model.dsMaxAbs = 0.05;

nls = NonLinearSolver('useLinesearch', false);
[wellSols, states, reports] = simulateScheduleAD(state0, model, schedule, 'nonlinearSolver', nls);


close all
p = (50:100:450)*barsa;
for pNo = 1:numel(p)
    plotSolventFluidProps(model, {'kr'}, {'O', 'G', 'S'}, 'pressure', p(pNo));
end

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
