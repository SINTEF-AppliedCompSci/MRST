function [eqs, names, types, state] = biotCompostionalEquations(model, state0, state, dt, drivingForces, varargin)
            
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
    
    [feqs, flux, fnames, ftypes] = model.FluxDiscretization.componentConservationEquations(model, state, state0, dt);
    src = model.FacilityModel.getComponentSources(state);
    % Assemble equations and add in sources
    [pressures, sat, mob, rho, X] = model.getProps(state, 'PhasePressures', 's', 'Mobility', 'Density', 'ComponentPhaseMassFractions');
    comps = cellfun(@(x, y) {x, y}, X(:, 1+model.water), X(:, 2+model.water), 'UniformOutput', false);
    
    
    feqs = model.addBoundaryConditionsAndSources(feqs, fnames, ftypes, state, pressures, sat, mob, rho, {}, comps, ...
                                                 drivingForces);
    
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

