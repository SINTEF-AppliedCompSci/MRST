%% Simple script for validation of parameter sensitivities by comparing 
mrstModule add ad-core ad-blackoil ad-props optimization spe10 deckformat

%% Setup simple model
nxyz = [ 10,  10,  2];
Dxyz = [400, 400, 10];
rng(0)
G    = computeGeometry(cartGrid(nxyz, Dxyz));
rock = getSPE10rock(1:nxyz(1), 101+(1:nxyz(2)), 1:nxyz(3));

% fluid
pRef  = 200*barsa;
fluid = initSimpleADIFluid('phases', 'WO',... 
                           'mu' , [.3, 3]*centi*poise,...
                           'rho', [1014, 859]*kilogram/meter^3, ...
                           'n', [2 2]);
fluid .krPts  = struct('w', [0 0 1 1], 'ow', [0 0 1 1]);
scaling = {'SWL', .1, 'SWCR', .2, 'SWU', .9, 'SOWCR', .1, 'KRW', .9, 'KRO', .8};

c = 5e-5/barsa;
p_ref = 200*barsa;
fluid.bO = @(p) exp((p - p_ref)*c);
model = GenericBlackOilModel(G, rock, fluid, 'gas', false);
model = imposeRelpermScaling(model, scaling{:});

%% wells/schedule
W = [];
% Injectors (lower-left and upper-right)
[wx, wy] = deal([1, nxyz(1)], [1, nxyz(2)]);
for k  = 1:2
    W = verticalWell(W, G, rock, wx(k), wy(k), 1:nxyz(3), 'Type' , 'rate', ...
                     'Val', 300*meter^3/day, 'Name', sprintf('I%d', k), ...
                     'comp_i', [1 0], 'Sign' , 1);
end
% Producers (upper-left and -right)
[wx, wy] = deal([1, nxyz(1)], [nxyz(2), 1]);
for k  = 1:2
    W = verticalWell(W, G, rock, wx(k), wy(k), 1:nxyz(3), 'Type' , 'bhp', ...
                     'Val', 100*barsa, 'Name', sprintf('P%d', k), ...
                     'comp_i', [1 0], 'Sign' , 1);
end
% Set up 4 control-steps each 150 days
schedule = simpleSchedule(rampupTimesteps(2*year, 30*day, 5), 'W', W);

%% run simulation
state0 = initState(G, W, 200*barsa, [0, 1]); 
% The accuracy in the gradient depend on the acuracy on the CNV tolerance
model.toleranceCNV = 1e-6;
[ws, states, r] = simulateScheduleAD(state0, model, schedule);

%% make a perturbed state for reference case
states_ref = states;
pflds = {'qWs', 'qOs'};
fac   = .1;
rng(0);
for tk = 1:numel(schedule.step.val)
    for wk = 3:numel(W)
        for fk = 1:numel(pflds)
            states_ref{tk}.wellSol(wk).(pflds{fk}) = ...
                (1+2*fac*(rand))*states{tk}.wellSol(wk).(pflds{fk});
        end
    end
end
%% parameter options
setup = struct('model', model, 'schedule', schedule, 'state0', state0);
nc =  model.G.cells.num;
nf =  numel(model.operators.T);
% select configuration for sensitivity computations
config = {...%name      include     scaling    boxlims     lumping     subset
          'porevolume',       1,   'linear',       [],          [],    1:100
          'conntrans',        1,   'log',          [],          [],       [] 
          'transmissibility', 1,   'log'   ,       [],          [],       [] 
          'swl',              1,   'linear',       [],  ones(nc,1),       []
          'swcr',             1,   'linear',       [],  ones(nc,1),       []
          'swu',              1,   'linear',       [],  ones(nc,1),       []
          'sowcr',            1,   'linear',       [],  ones(nc,1),       []
          'krw',              1,   'linear',       [],  ones(nc,1),       []
          'kro',              1,   'linear',       [],  ones(nc,1),       []
          'sw',               1,   'linear',   [0 .3],          [],       []  
          'pressure'          1,   'linear',       [],          [],       []};
parameters = [];
for k = 1:size(config,1)
    if config{k, 2} ~= 0
        parameters = addParameter(parameters, setup, 'name',    config{k,1}, ...
                                                     'scaling', config{k,3}, ...
                                                     'boxLims', config{k,4}, ...
                                                     'lumping', config{k,5}, ...
                                                     'subset',  config{k,6});
    end
end

%% Setup function handle to evaluateMatch
u = getScaledParameterVector(setup, parameters);
% Define weights for objective
weighting =  {'WaterRateWeight',  (300/day)^-1, ...
              'OilRateWeight',    (300/day)^-1, ...
              'BHPWeight',        (500*barsa)^-1};
          
obj = @(model, states, schedule, states_ref1, tt, tstep, state) matchObservedOW(model, states, schedule, states_ref,...
       'computePartials', tt, 'tstep', tstep, weighting{:}, 'state', state, 'from_states', false);   

f = @(u)evaluateMatch(u, obj, setup ,parameters,  states_ref, 'enforceBounds', false);       
       
%% Check gradient in random direction and compare to numerical
gp    = [];
[v,g] = f(u);
% optimal perturbation factor depends on combination of parameters, run
% with a few different perturbations
epsilons = 10.^(-(6:10));
for fac = epsilons
    rng(0);
    du  = fac*(rand(size(u)));
    vp = f(u+du);
    fprintf('Directional gradient obtained by perturbation: %e\n', (vp-v)/norm(du));
    fprintf('Directional gradient obtained by adjoint:      %e\n', g'*du/norm(du));
    gp = [gp, (vp-v)/norm(du)]; %#ok
end
ga =  g'*du/norm(du);
fprintf('Directional gradients\n')
fprintf('==============================================\n')
fprintf('Perturbation:    Relative gradient differences\n')
fprintf('       %5.0e            %7.2e         \n', [epsilons; abs(ga-gp)/abs(ga)])

%% Copyright Notice
%
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
