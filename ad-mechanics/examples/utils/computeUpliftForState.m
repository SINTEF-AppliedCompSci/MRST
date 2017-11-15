function uplift = computeUpliftForState(model, state, topnode, ComputePartials)
% Compute the vertical displacement at the node given by topnode. This function
% is used in runAdjointExample

    
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
    

    wellSol = state.wellSol;
        
    % This is specific to oil water fluid model
    % TODO: extend to other fluid models
    [p, sW, wellSol, xd] = model.getProps(state, 'pressure', ...
                                                 'water', ...
                                                 'wellSol', ...
                                                 'xd');
    [wellVars, wellVarNames, wellMap] = ...
        model.fluidModel.FacilityModel.getAllPrimaryVariables(wellSol);
    
    if ComputePartials
        [p, sW, wellVars{:}, xd] = initVariablesADI(p, sW, wellVars{:}, xd);
    end
    uplift = xd(indlift);

    
end