% run setup
optimizeSleipnerSetup;
%[wellSols, states, sim_report] = simulateScheduleAD(initState, smodel, schedule);
schedule_base = schedule;
states_base   = states; 

%% Define objective function
matchToPlume = true;
if matchToPlume
    obj = @(wellSols, states, schedule, varargin)...
            matchToDataSens(smodel, wellSols, states, schedule, newplumes, varargin{:});
else
    % run perturbed model
    s_pert = schedule;
    zp = 0*(rand(model.G.cells.num,1)-.5);
    [qp, rp, prmp, porp] = deal(1, 1, 1, 1.2);
    for i=1:numel(schedule.control)
        s_pert.control(i).dz       = zp;
        s_pert.control(i).W.val    = qp*s_pert.control(i).W.val;
        s_pert.control(i).rhofac   = rp;
        s_pert.control(i).permfac  = prmp;
        s_pert.control(i).porofac  = porp;
    end
    [wellSols_pert, states_pert] = simulateScheduleAD(initState, smodel, s_pert);
    obj = @(wellSols, states, schedule, varargin)...
        matchToSim(smodel, wellSols, states, schedule, states_pert, varargin{:});
end
v_base = obj(wellSols, states, schedule);
v_base = cell2mat(v_base);

%% Set up box limits for scaling and define function evaluation handle
dzLims  = [-1 1];
qLims   = [300 550]/day;
rhoLims = [.5 2];
permLims = [.5 2];
poroLims = [.5 2];
scaling.boxLims = [ones(model.G.cells.num, 1)*dzLims; qLims; rhoLims; permLims; poroLims];
scaling.obj     = max(sum(v_base), 1);  

% Get initial scaled controls 
%u_base = schedule2control(schedule_base, scaling);

% Convert schedule-params to control vector (assume all control-steps have same)
q_mean = mean(arrayfun(@(x)x.W.val, schedule_base.control));
c  = schedule.control(1);
u  = [c.dz; q_mean; c.rhofac; c.permfac; c.porofac];
[umin, umax] = deal(scaling.boxLims(:,1), scaling.boxLims(:,2));
u  = (u-umin)./(umax-umin); 

%% Get function handle for objective evaluation
f = @(u)evalObjectiveAndSens(u, obj, initState, smodel, schedule, scaling);

[v, u_opt, history] = unitBoxBFGS(u, f);

% recompute opt
[v_opt, der_opt, wellSols_opt, states_opt] = evalObjectiveAndSens(u_opt, obj, initState, smodel, schedule, scaling);
% scale back:
v_opt = -v_opt*scaling.obj;

[umin, umax] = deal(scaling.boxLims(:,1), scaling.boxLims(:,2));
u_opt = u_opt.*(umax-umin)+umin;

fprintf('rate multiplier: %4.2f\n', u_opt(end-3)/q_mean);
fprintf('rho mutiplier: %4.2f, perm mutiplier: %4.2f, poro mutiplier: %4.2f\n', u_opt(end-2), u_opt(end-1),u_opt(end))

figure
subplot(1,3,1)
title('Optimal dz')
plotCellData(smodel.G, u_opt(1:end-4), 'EdgeColor', 'none')
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
