%% Optimize NPV for a simple faulted model 
% This example covers the steps for setting up and running a simulation-based 
% optimization problem.

%% Construct faulted grid, and populate with layered permeability field  
mrstModule add ad-core ad-blackoil ad-props optimization

grdecl = makeModel3([20, 20, 4], [500, 500, 8]*meter);
G = processGRDECL(grdecl);
G = computeGeometry(G);

% Set up permeability based on K-indices
[I, J, K] = gridLogicalIndices(G);

px       = 100*milli*darcy*ones(G.cells.num,1);
px(K==2) = 500*milli*darcy;
px(K==3) = 50*milli*darcy;

% Introduce anisotropy by setting K_x = 5*K_z.
perm = [px, px, 0.1*px];
rock = makeRock(G, perm, 0.2);

%% Define wells 
% We define two vertical producers and two vertical injectors. Rates
% correspond to one reservoir pore volume injected during 640 days
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
figure,            
plotCellData(G, log(rock.perm(:,1)));
plotWell(G, W), view([1 1 1])
%% Define base schedule
% We set up 4 control-steps each 160 days. We explicitly set shorter time
% steps during start.
ts = { [1 2 5 7 10 15 20 20 20 20 20 20]'*day, ...
                      repmat(160/7, 7, 1)*day, ...
                      repmat(160/7, 7, 1)*day, ...
                      repmat(160/7, 7, 1)*day};
       
numCnt = numel(ts);
[schedule.control(1:numCnt).W] = deal(W);
schedule.step.control = rldecode((1:4)', cellfun(@numel, ts));
schedule.step.val     = vertcat(ts{:});

%% Set fluid properties (slightly compressible oil)
fluid = initSimpleADIFluid('mu',    [.3, 5, 0]*centi*poise, ...
                           'rho',   [1000, 700, 0]*kilogram/meter^3, ...
                           'n',     [2, 2, 0]);
c = 1e-5/barsa;
p_ref = 300*barsa;
fluid.bO = @(p) exp((p - p_ref)*c);

%% Initialize model and run base-case simulation
model  = TwoPhaseOilWaterModel(G, rock, fluid);
state0 = initResSol(G, p_ref, [0 1]);
schedule_base = schedule;
[wellSols_base, states_base] = simulateScheduleAD(state0, model, schedule_base);

%% Set prices ($/stb), discount rate (/year) and plot evolutiuon of NPV 
npvopts     = {'OilPrice',             60.0 , ...
               'WaterProductionCost',   7.0 , ...
               'WaterInjectionCost',    7.0 , ...
               'DiscountFactor',        0.1 };
v_base  = NPVOW(model, states_base, schedule_base, npvopts{:});
v_base  = cell2mat(v_base);
t_base  = cumsum(schedule_base.step.val);
figure, plot(convertTo(t_base,day), cumsum(v_base), '-o', 'LineWidth', 2);
title('Base run evolution NPV'), xlabel('days')

%% Set up box limits for scaling and define function evaluation handle
li = [  10, 300]/day;  % Injector limits  
lp = [-300, -10]/day;  % Producer limits 
scaling.boxLims = [li;li;lp;lp];  % control scaling  
scaling.obj     = sum(v_base);    % objective scaling    
% Get initial scaled controls 
u_base = schedule2control(schedule_base, scaling);
% Define objective function with above options
obj = @(model, states, schedule, varargin)NPVOW(model, states, schedule, varargin{:}, npvopts{:});
% Get function handle for objective evaluation
f = @(u)evalObjective(u, obj, state0, model, schedule, scaling);

%% Define linear equality and inequality constraints, and run optimization
% Constraints are applied to all control steps such for each step i, we
% enforce 
%       linEq.A*u_i   = linEq.b
%       linIneq.A*u_i = linIneq.b
%
% All rates should add to zero to preserve reservoir pressure (I1, I2, P1, P2)
linEq = struct('A', [1 1 1 1], 'b', 0);
% We also impose a total water injection constraint <= 500/day
linIneq = struct('A', [1 1 0 0], 'b', 500/day);  
% Constraints must be scaled!
linEqS   = setupConstraints(linEq,   schedule, scaling);
linIneqS = setupConstraints(linIneq, schedule, scaling);
% Run optimization with default options
[v, u_opt, history] = unitBoxBFGS(u_base, f, 'linEq', linEqS, 'linIneq', linIneqS);

%% Evaluate evolution of NPV for optimal schedule and compare with base
schedule_opt = control2schedule(u_opt, schedule, scaling);
[wellSols_opt, states_opt] = simulateScheduleAD(state0, model, schedule_opt);
v_opt  = NPVOW(model, states_opt, schedule_opt, npvopts{:});
v_opt  = cell2mat(v_opt);
t_opt  = cumsum(schedule_opt.step.val);
figure, plot(convertTo([t_base, t_opt],day), cumsum([v_base, v_opt]), '-o', 'LineWidth', 2);
title('Evolution NPV'), legend('Base', 'Optimal', 'Location', 'se'), xlabel('days')

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
