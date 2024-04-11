function [misfit, varargout] = evaluatePlumeMatch(pvec, obj, setup, parameters, plumes, scaling)

    % strictly enforce [0, 1] bounds on pvec, and split it up into cells
    % corresponding to each parameter group
    nparam = cellfun(@(x) x.nParam, parameters);
    pvec = max(0, min(1, pvec));
    pvec = mat2cell(pvec, nparam, 1);
    
    % Create new setup, and set parameter values
    pval = cell(size(parameters));
    setupNew = setup;
    setupNew.model.FlowDiscretization = [];
    setupNew.model.FlowPropertyFunctions = []; % @@ Do we need to reset more of these?
    for k = 1:numel(parameters)
        pval{k} = parameters{k}.unscale(pvec{k});
        setupNew = parameters{k}.setParameter(setupNew, pval{k});
    end
    
    % Run simulation on the new setup, and computing the corresponding misfit
    [wellSols, states] = simulateScheduleAD(setupNew.state0, setupNew.model, setupNew.schedule);
    
    misfitVals = obj(setupNew.model, states, plumes);
    misfit = -1 * sum(vertcat(misfitVals{:})) / scaling;
    
    if nargout == 1
        return
    end

    %% if we got to this point, we also need to compute partials
    objh = @(tstep, model, state) obj(setupNew.model, state, plumes(tstep));
    names = applyFunction(@(x)x.name, parameters);
    scaledGradient = cell(numel(names), 1);
    
    setupNew.model = setupNew.model.validateModel();
    gradient = computeSensitivitiesAdjointAD(setupNew, states, parameters, objh);
    for k = 1:numel(names)
        scaledGradient{k} = -1 * parameters{k}.scaleGradient(gradient.(names{k}), pval{k});
    end
    
    %% Include additional varargouts, if requested
    varargout{1} = vertcat(scaledGradient{:}) / scaling;
    if nargout > 2
        [varargout{2:3}] = deal(wellSols, states);
    end
    if nargout > 4
        varargout{4} = setupNew;
    end
end
