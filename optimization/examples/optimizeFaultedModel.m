grdecl = makeModel3([20, 20, 4], [500, 500, 8]*meter);
G = processGRDECL(grdecl);
G = computeGeometry(G);

% Set up permeability based on K-indices
[I, J, K] = gridLogicalIndices(G);

layer2 = K == 2;

px       = 100*milli*darcy*ones(G.cells.num,1);
px(K==2) = 500*milli*darcy;
px(K==3) = 50*milli*darcy;

% Introduce anisotropy by setting K_x = 10*K_z.
perm = [px, px, 0.1*px];
rock = makeRock(G, perm, 0.2);


%% Define wells and simulation schedule

simTime = 640*day;
pv = poreVolume(G, rock);
injRate = 1*sum(pv)/simTime;

% Place wells
[nx, ny] = deal(G.cartDims(1), G.cartDims(2));
ppos = [nx-2, 4; 4   , ny-4];
ipos = [4   , 5; nx-3, ny-3];
offset = 5;
W = [];
% Injectors            
W = verticalWell(W, G, rock, ipos(1,1), ipos(1,2), [], 'sign', 1,...
                'Name', 'I1', 'comp_i', [1 0], 'Val', injRate/2, 'Type', 'rate');

W = verticalWell(W, G, rock, ipos(2,1), ipos(2,2), [], 'sign', 1, ...
                'Name', 'I2', 'comp_i', [1 0], 'Val', injRate/2, 'Type', 'rate');
% Producers
W = verticalWell(W, G, rock, ppos(1,1), ppos(1,2), [], 'sign', -1, ...
                'Name', 'P1', 'comp_i', [0 1], 'Val', -injRate/2, 'Type', 'lrat');
W = verticalWell(W, G, rock, ppos(2,1), ppos(2,2), [], 'sign', -1, ...
                'Name', 'P2', 'comp_i', [0 1], 'Val', -injRate/2, 'Type', 'lrat');
            
%% Set schedule
% Set up 4 control-steps each 140 days
ts = { [1 2 5 7 10 15 20 20 20 20 20 20]'*day, ...
                      repmat(160/7, 7, 1)*day, ...
                      repmat(160/7, 7, 1)*day, ...
                      repmat(160/7, 7, 1)*day};
       
numCnt = numel(ts);
[schedule.control(1:numCnt).W] = deal(W);
schedule.step.control = rldecode((1:4)', cellfun(@numel, ts));
schedule.step.val     = vertcat(ts{:});

%% Set fluid
fluid = initSimpleADIFluid('mu',    [.3, 5, 0]*centi*poise, ...
                           'rho',   [1000, 700, 0]*kilogram/meter^3, ...
                           'n',     [2, 2, 0]);

c = 1e-5/barsa;
p_ref = 300*barsa;
fluid.bO = @(p) exp((p - p_ref)*c);

%% Initial state:
state0 = initResSol(G, p_ref, [0 1]);

%% Initialize model and run simulation
model = TwoPhaseOilWaterModel(G, rock, fluid);
[wellSols, states] = simulateScheduleAD(state0, model, schedule);

%% Set prices (/stb), discount rate (/year) and plot evolutiuon of NPV 
npvopts     = {'OilPrice',             60.0 , ...
               'WaterProductionCost',   7.0 , ...
               'WaterInjectionCost',    7.0 , ...
               'DiscountFactor',        0.1 };
%                    'ComputePartials',      false, ...
%                    'tStep' ,               []);
vals  = NPVOW(G, wellSols, schedule, npvopts{:});
vals  = cell2mat(vals);
times = cumsum(schedule.step.val);
plot(convertTo(times,day), cumsum(vals));

%% Set up box limits for scaling and define function evaluation 
li = [  10, 300]/day;  % Injector limits  
lp = [-300, -10]/day;  % Producer limits 
scaling.boxLims = [li;li;lp;lp];  % control scaling  
scaling.obj     = sum(vals);   % objective scaling    
% Get initial scaled controls
uInit = schedule2control(schedule, scaling);
% Define objective function with above options
obj = @(wellSols, schedule, varargin)NPVOW(G, wellSols, schedule, varargin{:}, npvopts{:});
% Get function handle for objective evaluation
f = @(u)evalObjective(u, obj, state0, model, schedule, scaling);

%% Define scaled linear equality and inequality constraints, and run optimization
% Sum of all rates should be zero to preserve reservoir pressure
linEq = struct('A', [1 1 1 1], 'b', 0);
% Impose total water injection constraint <= 500/day
linIneq = struct('A', [1 1 0 0], 'b', 500/day);  
% Constraints must be scaled!
linEqS   = setupConstraints(linEq,   schedule, scaling);
linIneqS = setupConstraints(linIneq, schedule, scaling);
% Run optimization
[v, u, history] = unitBoxBFGS(uInit, f, 'linEq', linEqS, 'linIneq', linIneqS);


            



