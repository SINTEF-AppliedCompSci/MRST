function grad = computeGradientPerturbationAD(state0, model, schedule, getObjective, varargin)
    % Compute numerical gradient w.r.t wells
    opt = struct('Verbose', mrstVerbose(),...
                 'scaling', []);
    opt = merge_options(opt, varargin{:});

    if ~isempty(opt.scaling)
        scalFacs = opt.scaling;
    else
        scalFacs.rate = 1; scalFacs.pressure = 1;
    end
    
    solve = @(schedule) simulateScheduleAD(state0, model, schedule, 'Verbose', opt.Verbose);
    schedule0 = schedule;

    % Solve entire schedule to get baseline of the objective function
    wellSols = solve(schedule0);

    % Set up objective function storage
    computeObj = @(ws) sum(cell2mat(getObjective(ws)));
    val0 = computeObj(wellSols);


    grad = cell(1, numel(schedule0.control));
    % Run a schedule per well, for each control step, perturbing the
    % control variable ever so slightly to get a local finite difference
    % approximation of the gradient.
    for cn = 1:numel(schedule.control)
        
        ctrl = schedule0.control(cn);
        nWell = numel(ctrl.W);
        grad{cn} = zeros(nWell, 1);

        dispif(opt.Verbose, 'Solving for control %d of %d', cn, numel(schedule.control));
        for k = 1:nWell
            dispif(opt.Verbose, 'Processing well %d of %d\n', k, numel(ctrl.W));

            w = ctrl.W(k);

            if strcmpi(w.type, 'bhp')
                e = scalFacs.pressure*1e-7;
            else
                e = scalFacs.rate*1e-7;
            end

            schedule = schedule0;
            w.val = w.val + e;
            schedule.control(cn).W(k) = w;
            
            wellSols = solve(schedule);
            valk = sum( cell2mat(getObjective(wellSols)));
            grad{cn}(k) = (valk-val0)/e;
        end
    end
end


