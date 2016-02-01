% script producing results for sim-matching example paper
optimizeSleipnerSetup;

% perturbed simulation
q_mean = mean(arrayfun(@(x)x.W.val, schedule.control));
[qp, rp, prmp, porp] = deal(1, .8, 1.2, 0.75);
zp = 4*(rand(model.G.cells.num,1)-.5);
dzLims  = [-3 3];
qLims   = [.995 1.005]*q_mean;
rhoLims = [.5 2];
permLims = [.5 2];
poroLims = [.5 2];

% cases match-instances
matchCases = {...
    [1 1 1 1 1 1 1 1 1 1]', ...
    [0 0 0 0 1 0 0 0 0 1]', ...
    [0 0 0 0 0 0 0 0 0 1]'};

for caseNo = 1:numel(matchCases)
    matchSteps = matchCases{caseNo};
    
    %[wellSols, states, sim_report] = simulateScheduleAD(initState, smodel, schedule);
    schedule_base = schedule;
    states_base   = states;
    
    model.nonlinearTolerance = 1e-9;
    
    
    %% Define objective function
    matchToPlume = false;
    if matchToPlume
        obj = @(wellSols, states, schedule, varargin)...
            matchToDataSens(smodel, wellSols, states, schedule, newplumes, varargin{:});
    else
        % run perturbed model
        s_pert = schedule;
        %zp = 2*(rand(model.G.cells.num,1)-.5);
        %[qp, rp, prmp, porp] = deal(1.2, .8, 1.14, 0.76);
        for i=1:numel(schedule.control)
            s_pert.control(i).dz       = zp;
            s_pert.control(i).W.val    = qp*s_pert.control(i).W.val;
            s_pert.control(i).rhofac   = rp;
            s_pert.control(i).permfac  = prmp;
            s_pert.control(i).porofac  = porp;
        end
        [wellSols_pert, states_pert] = simulateScheduleAD(initState, smodel, s_pert);
        obj = @(wellSols, states, schedule, varargin)...
            matchToSim(smodel, wellSols, states, schedule, states_pert, 'matchSteps', matchSteps, varargin{:});
    end
    v_base = obj(wellSols, states, schedule);
    v_base = cell2mat(v_base);
    
    %% Set up box limits for scaling and define function evaluation handle
    scaling.boxLims = [ones(model.G.cells.num, 1)*dzLims; qLims; rhoLims; permLims; poroLims];
    scaling.obj     = max(sum(abs(v_base)), 1);
    
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
    
    % store results
    res = struct('perm', u_opt(end-1), ...
                 'poro', u_opt(end), ...
                 'rho', u_opt(end-2), ...
                 'rate', u_opt(end-3)/q_mean, ...
                 'dz', u_opt(1:end-4), ...
                 'v_base', sum(v_base), ...
                 'v_opt', v_opt, ...
                 'scaling', scaling, ...
                 'states_opt', {states_opt});
    save(['SimMatch', num2str(caseNo)], 'res') 
end

