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

obj = @(wellSols, schedule, varargin)NPVOW(G, wellSols, schedule, varargin{:}, npvopts{:});
% Get function handle for objective evaluation
f = @(u)evalObjective(u, obj, state0, model, schedule, scaling);

%% Run optimization with default options
[v, u_opt, history] = unitBoxBFGS(u_base, f);
schedule_opt = control2schedule(u_opt, schedule, scaling);
pth = fullfile(mrstPath('optimization'), 'examples', 'model2D', 'schedule_opt.mat');
save(pth, 'schedule_opt')
