%% Simple script for validation of parameter sensitivities by comparing 
mrstModule add ad-core ad-blackoil ad-props optimization spe10 deckformat

%% Setup simple model
nxyz = [ 10,  10,  1];
Dxyz = [400, 400, 10];
rng(0)
G    = computeGeometry(cartGrid(nxyz, Dxyz));
rock = getSPE10rock(1:nxyz(1), 101+(1:nxyz(2)), 1:nxyz(3));
rock.perm = ones(size(rock.perm))*1*darcy;

% fluid
pRef  = 200*barsa;
fluid = initSimpleADIFluid('phases', 'WO',... 
                           'mu' , [.3, 3]*centi*poise,...
                           'rho', [1014, 859]*kilogram/meter^3, ...
                           'n', [2 2]);
fluid .krPts  = struct('w', [0 0 1 1], 'ow', [0 0 1 1]);

c = 5e-5/barsa;
p_ref = 200*barsa;
fluid.bO = @(p) exp((p - p_ref)*c);
model = GenericBlackOilModel(G, rock, fluid, 'gas', false);

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
%schedule = simpleSchedule([1]*day, 'W', W);

%% run simulation
state0 = initState(G, W, 200*barsa, [0, 1]); 
% The accuracy in the gradient depend on the acuracy on the CNV tolerance
model.toleranceCNV = 1e-8;
[ws, states, r] = simulateScheduleAD(state0, model, schedule);

%% make a perturbed output as reference case
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
% transmissibility
parameters{1} = ModelParameter(setup, 'name', 'transmissibility', ...
                                      'type', 'multiplier');
                                                

%% Setup function handle to evaluateMatch
u = getScaledParameterVector(setup, parameters);
% Define weights for objective
weighting =  {'WaterRateWeight',  (300/day)^-1, ...
              'OilRateWeight',    (300/day)^-1, ...
              'BHPWeight',        (500*barsa)^-1};
          
% 1. gradient case - objective is sum of mismatches squared          
obj1 = @(model, states, schedule, states_ref1, tt, tstep, state) matchObservedOW(model, states, schedule, states_ref,...
       'computePartials', tt, 'tstep', tstep, weighting{:}, 'state', state, 'from_states', false, 'mismatchSum', true);   
f1 = @(u)evaluateMatch(u, obj1, setup ,parameters,  states_ref, 'enforceBounds', false);     

% 2. sensitivity matrix case - objective computes vector of all mismatches
obj2 = @(model, states, schedule, states_ref1, tt, tstep, state) matchObservedOW(model, states, schedule, states_ref,...
       'computePartials', tt, 'tstep', tstep, weighting{:}, 'state', state, 'from_states', false, 'mismatchSum', false);
f2 = @(u)evaluateSensitivityMatrix(u, obj2, setup ,parameters,  states_ref, 'enforceBounds', false);       

%% Check gradient in random direction and compare to numerical
% compute (negative) sum of squared mismatches and (negative) gradient
[v1, g] = f1(u);

% compute vector of squared mismatches and Jacobian
[v2, J] = f2(u);

% check that both produce the same gradient, i.e J'*v2 = -g
fprintf('Relative gradient mismatch: %e\n', norm(J'*v2 + g)/norm(g))
