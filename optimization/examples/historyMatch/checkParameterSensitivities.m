%% Simple script for validation of parameter sensitivities

mrstModule add ad-core ad-blackoil ad-props optimization spe10 deckformat


%% Setup simple model
nxyz = [ 10,  10,  2];
Dxyz = [400, 400, 10];
rng(0)
G    = computeGeometry(cartGrid(nxyz, Dxyz));
rock = getSPE10rock(1:nxyz(1), 101+(1:nxyz(2)), 1:nxyz(3));

% fluid
pRef  = 200*barsa;
fluid = initSimpleADIFluid('phases', 'WO',... %Fluid phase
                           'mu' , [.3, 3]*centi*poise     , ...%Viscosity
                           'rho', [1014, 859]*kilogram/meter^3, ...%Surface density [kg/m^3]
                           'n', [2 2]);
fluid = fluid;
fluid .krPts  = struct('w', [0 0 1 1], 'ow', [0 0 1 1]);
scaling = {'SWL', .1, 'SWCR', .2, 'SWU', .9, 'SOWCR', .1, 'KRW', .9, 'KRO', .8};

c = 1e-5/barsa;
p_ref = 200*barsa;
fluid.bO = @(p) exp((p - p_ref)*c);
model = GenericBlackOilModel(G, rock, fluid, 'gas', false, 'OutputStateFunctions', false);
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
prob = struct('model', model, 'schedule', schedule, 'state0', state0);

n_cells =  model.G.cells.num;
n_faces =  length(model.operators.T);
% Fluid Parameters
 parameters{1} = ModelParameter(prob, 'name', 'swl','lumping',ones(n_cells,1));
 parameters{2} = ModelParameter(prob, 'name', 'swcr','lumping',ones(n_cells,1));
 parameters{3} = ModelParameter(prob, 'name', 'swu','lumping',ones(n_cells,1));
 parameters{4} = ModelParameter(prob, 'name', 'kro','lumping',ones(n_cells,1));
 parameters{5} = ModelParameter(prob, 'name', 'krw','lumping',ones(n_cells,1));
% Well, porevolume and transmisibility
 parameters{6} = ModelParameter(prob, 'name', 'conntrans');
 parameters{7} = ModelParameter(prob, 'name', 'porevolume','lumping',[ones(round(n_cells/3),1);2*ones(n_cells-round(n_cells/3),1)]);
 parameters{8} = ModelParameter(prob, 'name', 'transmissibility','lumping',[1:n_faces]');

%% 
values = applyFunction(@(p)p.getParameterValue(prob), parameters);
% scale values
u = cell(size(values));
for k = 1:numel(u)
    u{k} = parameters{k}.scale(values{k});
end
u = vertcat(u{:});


% Defining the weights to evaluate the match
weighting =  {'WaterRateWeight',  (300/day)^-1, ...
              'OilRateWeight',    (300/day)^-1, ...
              'BHPWeight',        (500*barsa)^-1};
          
obj = @(model, states, schedule, states_ref1, tt, tstep, state) matchObservedOW(model, states, schedule, states_ref,...
       'computePartials', tt, 'tstep', tstep, weighting{:}, 'state',state, 'from_states', false);   

f = @(u)evaluateMatch_simple(u, obj, state0, model, schedule, 1 ,parameters,  states_ref);       
       
%% Check gradient in random direction and compare to numerical
[v,g] = f(u);
% optimal perturbation factor depends on combination of parameters 
fac = 1e-7;
du  = fac*(rand(size(u))-.5);
du  = max(du, -u);
du  = min(du, 1-u);
vp = f(u+du);       
fprintf('\nDirectional gradient obtained by perturbation: %e\n', (vp-v)/norm(du));
fprintf('Directional gradient obtained by adjoint:      %e\n', g'*du/norm(du));

%%
