mrstModule add ad-core ad-eor ad-blackoil ad-props sequential matlab_bgl

gravity reset on

%%

pth = mrstPath('mrst-solvent');

load([pth, '/code/FourPhaseSolvent/examples/data/norne.mat']);

G = computeGeometry(G);

mrstModule add deckformat

%%


inj = [9, 15; ...
       26, 15; ...
       36, 80; ...
       10, 85; ...
       24, 30; ...
       14, 52; ...
       18, 80;
       23, 66];

nInj = size(inj,1);
injectors = cell(nInj,1);
for i = 1:nInj
    W = verticalWell([], G, rock, inj(i, 1), inj(i, 2), []);
    injectors{i} = W.cells;
end
   
prod = [10, 66; ...
        12, 32; ...
        22, 49; ...
        13, 91; ...
        37, 95;
        35, 64];
nProd = size(prod,1);
producers = cell(nProd, 1);
for i = 1:size(prod,1)
    W = verticalWell([], G, rock, prod(i, 1), prod(i, 2), []);
    producers{i} = W.cells;
end
    
%%

gravity reset on

fluid = initSimpleADIFluid('n'     , [2, 2, 2], ...
                           'rho'   , [1000, 800, 100]*kilogram/meter^3, ...
                           'phases', 'WOG', ...
                           'mu'    , [1, 10, 2]*centi*poise, ...
                           'c'     , [1e-6, 1e-6, 1e-5]/barsa);

sOres_i = 0.38;
sOres_m = 0.08;
fluid = addSolventProperties(fluid, 'n', 2, ...
                                    'rho', 100*kilogram/meter^3, ...
                                    'mixPar', 2/3, ...
                                    'mu'    , 1*centi*poise, ...
                                    'sOres_i', 0.38, ...
                                    'sOres_m', 0.0 , ...
                                    'c'     , 1e-6/barsa);
                                
model = FourPhaseSolventModel(G, rock, fluid);
model.extraStateOutput = true;

Tperiod = 1*year;
Ttot = 4*Tperiod;
pv = 2*sum(poreVolume(G, rock));
rate = (pv/Ttot)/nInj;
nStep = 500;
nCycles = 5;
[scheduleWAG, W_G, W_W] = makeWAGschedule(model, injectors, producers, nCycles, ...
    'T', 2*Tperiod, 'gRate', rate, 'wRate', rate, 'nStep', 2*nStep);

dT = rampupTimesteps(Tperiod, Tperiod/nStep,0);
step.val = [dT; scheduleWAG.step.val; dT];
step.control = [2*ones(numel(dT),1); scheduleWAG.step.control; 2*ones(numel(dT),1)];

schedule.control = scheduleWAG.control;
schedule.step = step;

state0 = initResSol(G, 100*barsa, [0 1 0 0]);
state0.wellSol = initWellSolAD(W_G, model, state0);

%%

getHandler = @(name) ResultHandler('dataFolder', ['norneWAG', ...
    '_nStep', num2str(nStep), '_nCycles', num2str(nCycles), ...
    '_sOres_i', num2str(sOres_i), '_sOres_m', num2str(sOres_m)], ...
    'dataPrefix', [name, '_step'], 'cleardir', false);
getHandler = @(name) [];

stateHandler = getHandler('state');
repHandler = getHandler('rep');

%%

[ws, states, reports] = simulateScheduleAD(state0, model, schedule, ...
    'OutputHandler', stateHandler, 'ReportHandler', repHandler);

%%

plotToolbar(G, states);

%%
sHandler = ResultHandler('dataFolder', 'norneWAG_nStep400_', 'dataPrefix', [name, '_full_', num2str(step), '_step'], 'cleardir', false);

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
