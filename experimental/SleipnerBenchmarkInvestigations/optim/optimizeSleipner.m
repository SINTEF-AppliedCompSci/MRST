% run setup
optimizeSleipnerSetup;
%[wellSols, states, sim_report] = simulateScheduleAD(initState, smodel, schedule);
schedule_base = schedule;
states_base   = states; 
v_base = matchToDataSens(smodel, wellSols, states_base, schedule_base, newplumes);
v_base = cell2mat(v_base);

%% Set up box limits for scaling and define function evaluation handle
dzLims  = [-1 1];
rhoLims = [.2 2.5];
permLims = [.2 2.5];
poroLims = [.2 2.5];
scaling.boxLims = [ones(model.G.cells.num, 1)*dzLims; rhoLims; permLims; poroLims];
scaling.obj     = sum(v_base);  

% Get initial scaled controls 
%u_base = schedule2control(schedule_base, scaling);

% Convert schedule-params to control vector (assume all control-steps have same)
c  = schedule.control(1);
u  = [c.dz; c.rhofac; c.permfac; c.porofac];
[umin, umax] = deal(scaling.boxLims(:,1), scaling.boxLims(:,2));
u  = (u-umin)./(umax-umin); 

% Define objective functi
obj = @(wellSols, states, schedule, varargin)...etr
            matchToDataSens(smodel, wellSols, states, schedule, newplumes, varargin{:});
% Get function handle for objective evaluation
f = @(u)evalObjectiveAndSens(u, obj, initState, smodel, schedule, scaling);

[v, u_opt, history] = unitBoxBFGS(u, f);

% recompute opt
[v_opt, der_opt, wellSols_opt, states_opt] = evalObjectiveAndSens(u_opt, obj, initState, smodel, schedule, scaling);
% scale back:
v_opt = -v_opt*scaling.obj;

[umin, umax] = deal(scaling.boxLims(:,1), scaling.boxLims(:,2));
u_opt = u_opt.*(umax-umin)+umin;

fprintf('rho mutiplier: %4.2f, perm mutiplier: %4.2f, poro mutiplier: %4.2f\n', u_opt(end-2), u_opt(end-1),u_opt(end))

figure
subplot(1,3,1)
title('Optimal dz')
plotCellData(smodel.G, u_opt(1:end-3), 'EdgeColor', 'none')
axis tight, colorbar
subplot(1,3,2)
title('Mismatch base')
mm0 = model.G.cells.H.*states_base{10}.s(:,2)/(1-fluid.res_water)-newplumes{10}.h;
plotCellData(smodel.G, mm0, 'EdgeColor', 'none')
caxis([-3 5])
axis tight, colorbar
subplot(1,3,3)
title('Mismatch optimal')
mm0 = model.G.cells.H.*states_opt{10}.s(:,2)/(1-fluid.res_water)-newplumes{10}.h;
plotCellData(smodel.G, mm0, 'EdgeColor', 'none')
caxis([-3 5])
axis tight, colorbar
