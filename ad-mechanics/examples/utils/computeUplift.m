function obj = computeUplift(model, states, schedule, topnode, varargin)
% Compute the average of the vertical displacement at the top of the domain 
% This function is used in runAdjointExample

    opt = struct('tStep', [], 'ComputePartials', false);
    opt = merge_options(opt, varargin{:});
    
    % get index in the mechanical displacement field xd where the vertical
    % displacement of the top node is stored.
    G = model.G;
    nx = G.cartDims(1);
    ny = G.cartDims(2);
    isdirdofs = model.mechModel.operators.isdirdofs;
    u = (1 : (G.griddim*G.nodes.num))';
    indlift = G.griddim*(topnode - 1) + 2;
    u = u(~isdirdofs);
    indlift = find(u == indlift);
    
    
    numSteps = numel(schedule.step.val);
    lastStep = numSteps;
    tSteps = opt.tStep;
    if isempty(tSteps) 
        % do all
        tSteps = (1 : numSteps)';
    else
        numSteps = 1;
    end
    obj = repmat({[]}, numSteps, 1);

    for step = 1 : numSteps
        state = states{tSteps(step)};
        wellSol = state.wellSol;
        
        % This is specific to oil water fluid model
        % TODO: extend to other fluid models
        [p, sW, wellSol, xd] = model.getProps(state, 'pressure', ...
                                                     'water', ...
                                                     'wellSol', ...
                                                     'xd');
        [wellVars, wellVarNames, wellMap] = ...
            model.fluidModel.FacilityModel.getAllPrimaryVariables(wellSol);
        
        if opt.ComputePartials
            [p, sW, wellVars{:}, xd] = initVariablesADI(p, sW, wellVars{:}, xd);
        end
        obj{step} = xd(indlift);
        if tSteps(step) ~= lastStep
            obj{step} = double2ADI(0, obj{step});
        end
    end
    
end