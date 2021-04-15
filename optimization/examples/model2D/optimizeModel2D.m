% optimizeModel2D - optimize NPV for example-model of this folder

mrstModule add ad-core ad-blackoil ad-props optimization spe10
setupModel2D

% Create model-object of class TwoPhaseOilWaterModel
model  = TwoPhaseOilWaterModel(G, rock, fluid);
% Set initial state and run simulation:
state0 = initResSol(G, 200*barsa, [0, 1]);

%% Set up box limits for scaling and define function evaluation handle
li = [0 400]/day;        % Injector limits  
lp = [100 250]*barsa;    % Producer limits 
scaling.boxLims = [li;li;lp;lp];  % control scaling  
scaling.obj     = 3.2e7;      % objective scaling    
% Get initial scaled controls 
u_base = schedule2control(schedule, scaling);
% Define objective function
d   = 0.05;    % yearly discount factor
ro  = 60;      % oil revenue/price ($/stb)
rwp =  6;      % water production handling costs ($/stb)
rwi =  6;      % water injection cost ($/stb) 
npvopts = {'OilPrice',             ro , ...
           'WaterProductionCost', rwp , ...
           'WaterInjectionCost',  rwi , ...
           'DiscountFactor',        d};

obj = @(model, states, schedule, varargin)NPVOW(model, states, schedule,varargin{:}, npvopts{:});
% Get function handle for objective evaluation
f = @(u)evalObjective(u, obj, state0, model, schedule, scaling);

%% Run optimization with default options
[v, u_opt, history] = unitBoxBFGS(u_base, f);
schedule_opt = control2schedule(u_opt, schedule, scaling);
pth = fullfile(mrstPath('optimization'), 'examples', 'model2D', 'schedule_opt.mat');
save(pth, 'schedule_opt')


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
