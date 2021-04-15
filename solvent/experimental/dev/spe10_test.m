mrstModule add solvent spe10 ad-props ad-blackoil mrst-gui

%% Use setupSPE10_AD to Fetch the SPE10 model
% We pick up only one layer 

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
                                    'sOres_i', 0.4 , ...
                                    'sOres_m', 0.23);
                                
model = FourPhaseSolventModel(G, rock, fluid);
model.extraStateOutput = true;


producers = {1, G.cartDims(1), G.cartDims(1)*G.cartDims(2) - G.cartDims(1) + 1, G.cartDims(1)*G.cartDims(2)};
injectors = {round(G.cartDims(1)*G.cartDims(2)/2 + G.cartDims(1)/2)};

W = [];
for pNo = 1:numel(producers)
    W = addWell(W, G, rock, producers{pNo}, ...
                'type', 'bhp', ...
                'val', 300*barsa, ...
                'comp_i', [0,0,0,1], ...
                'sign', 1);
end
for iNo = 1:numel(injectors)
    W = addWell(W, G, rock, injectors{iNo}, ...
                'type', 'bhp', ...
                'val', 50*barsa, ...
                'comp_i', [1,0,0,0], ...
                'sign', -1);
end

T  = 5*year;
dT = 30*day;
dt = rampupTimesteps(T, dT);

schedule = simpleSchedule(dt, 'W', W);

state0         = initResSol(G, 100*barsa, [0 1 0 0]);
state0.wellSol = initWellSolAD(W, model, state0);

%%

[ws, states, rep] = simulateScheduleAD(state0, model, schedule);

%%
plotToolbar(G, states);

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
