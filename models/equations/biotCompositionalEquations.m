function [eqs, names, types, state] = biotCompositionalEquations(model, state0, state, dt, drivingForces, varargin)
            
    G     = model.G;
    fluid = model.fluid;
    op    = model.operators;
    
    pv         = op.pv;
    momentop   = op.momentop;
    mechdirop  = op.mechDirichletop;
    
    
    u  = model.getProp(state, 'u');
    bp = model.getProp(state, 'bp');
    lm = model.getProp(state, 'lambdamech');
    
    u0  = model.getProp(state0, 'u');
    bp0 = model.getProp(state0, 'bp');
    lm0 = model.getProp(state0, 'lambdamech');
    
    p = model.getProp(state, 'pressure');

    fac = 1 / (1e6 * mean(G.cells.volumes));
    meqs{1} = fac*momentop(u, p, lm);
    meqs{2} = fac*mechdirop(u, p, lm);
    meqs{3} = p - bp;

    mnames = {'momentum', 'mechbcs', 'coupling'};
    mtypes  = {'cells', 'bc', 'cells'};
    
    [feqs, flux, fnames, ftypes] = model.FlowDiscretization.componentConservationEquations(model, state, state0, dt);
    src = model.FacilityModel.getComponentSources(state);
    
    % Treat source or bc terms
    if ~isempty(drivingForces.bc) || ~isempty(drivingForces.src)
        [pressures, sat, mob, rho, rs, rv] = model.getProps(state, 'PhasePressures', 's', 'Mobility', 'Density', 'Rs', 'Rv');
        dissolved = model.getDissolutionMatrix(rs, rv);
        feqs = model.addBoundaryConditionsAndSources(feqs, fnames, ftypes, state, ...
                                                    pressures, sat, mob, rho, ...
                                                    dissolved, {}, ...
                                                    drivingForces);
    end

    % Add aquifer contributions if any.
    if ~isempty(model.AquiferModel)
        feqs = addAquifersContribution(model.AquiferModel, feqs, fnames, state, dt);
    end

    % Add sources
    feqs = model.insertSources(feqs, src);
    % Assemble equations
    for i = 1:numel(feqs)
        feqs{i} = model.operators.AccDiv(feqs{i}, flux{i});
    end
    
    % Get facility equations
    [weqs, wnames, wtypes, state] = model.FacilityModel.getModelEquations(state0, state, dt, drivingForces);
    
    eqs = [meqs, feqs, weqs];
    names = [mnames, fnames, wnames];
    types = [mtypes, ftypes, wtypes];

    
end

