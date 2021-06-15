mrstModule add ad-core ad-eor ad-blackoil ad-props sequential matlab_bgl

gravity reset on

n = 3;
G = computeGeometry(cartGrid([n,1,1]));
rock = makeRock(G, 100*milli*darcy, 1);

T = 1*year;
pv = poreVolume(G, rock);
injRate = 2*sum(pv)/T;
nStep = 10;
dT = rampupTimesteps(T, T/nStep);

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
                                    'sOres_m', 0.1);
                                
model = FourPhaseSolventModel(G, rock, fluid);
model.extraStateOutput = true;

W = [];
W = addWell(W, G, rock, 1, 'comp_i', [1, 0, 0, 0], 'type', 'rate', 'val', injRate);
W = addWell(W, G, rock, G.cells.num, 'comp_i', [1, 0, 0, 0], 'type', 'bhp', 'val', 0);
schedule = simpleSchedule(dT, 'W', W);

%%

state0 = initResSol(G, 100*barsa, [0 1 0 0]);
state0.wellSol = initWellSolAD(W, model, state0);

[wsFP, statesFP, reportsFP] = simulateScheduleAD(state0, model, schedule);

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
