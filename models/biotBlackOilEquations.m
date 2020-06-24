function [eqs, names, types, state] = biotBlackOilEquations(model, state0, state, dt, drivingForces, varargin)
            
    G     = model.G;
    fluid = model.fluid;
    op    = model.operators;
    
    pv         = op.pv;
    divKgradop = op.divKgradop;
    divuop     = op.divuop;
    momentop   = op.momentop;
    mechdirop  = op.mechDirichletop;
    fluiddirop = op.fluidDirichletop;
    
    extforce = drivingForces.extforce;
    
    u  = model.getProp(state, 'u');
    p  = model.getProp(state, 'p');
    lm = model.getProp(state, 'lambdamech');
    lf = model.getProp(state, 'lambdafluid');
    
    u0  = model.getProp(state0, 'u');
    p0  = model.getProp(state0, 'p');
    lm0 = model.getProp(state0, 'lambdamech');
    
    c = fluid.c;
    fac =  1 / (1e6 * mean(G.cells.volumes));
    
    % dohorriblehack = true;
    % if dohorriblehack
    %     if all(u0 == 0) & all(p0 == 0)
    %         divu0 = 0;
    %     else
    %         divu0 = divuop(u0, p0, lm0);
    %     end
    % end 
    
    divu  = model.getProp(state, 'Dilation');
    divu0 = model.getProp(state0, 'Dilation');
    
    eqs{1} = fac*momentop(u, p, lm, extforce);
    eqs{2} = 1/dt.*(pv.*c.*(p - p0) + (divu - divu0)) + divKgradop(p, lf);
    eqs{3} = fac*mechdirop(u, p, lm);
    eqs{4} = fluiddirop(p, lf);
    
    
    [eqs, flux, names, types] = model.FluxDiscretization.componentConservationEquations(model, state, state0, dt);
    src = model.FacilityModel.getComponentSources(state);
    % Treat source or bc terms
    if ~isempty(drivingForces.bc) || ~isempty(drivingForces.src)
        [pressures, sat, mob, rho, rs, rv] = model.getProps(state, 'PhasePressures', 's', 'Mobility', 'Density', 'Rs', 'Rv');
        dissolved = model.getDissolutionMatrix(rs, rv);
        eqs = model.addBoundaryConditionsAndSources(eqs, names, types, state, ...
                                                    pressures, sat, mob, rho, ...
                                                    dissolved, {}, ...
                                                    drivingForces);
    end

    % Add aquifer contributions if any.
    if ~isempty(model.AquiferModel)
        eqs = addAquifersContribution(model.AquiferModel, eqs, names, state, dt);
    end

    % Add sources
    eqs = model.insertSources(eqs, src);
    % Assemble equations
    for i = 1:numel(eqs)
        eqs{i} = model.operators.AccDiv(eqs{i}, flux{i});
    end
    
    % Get facility equations
    [weqs, wnames, wtypes, state] = model.FacilityModel.getModelEquations(state0, state, dt, drivingForces);
    eqs = [eqs, weqs];
    names = [names, wnames];
    types = [types, wtypes];
    names = {'momentum', 'mass', 'bcmech', 'bcfluid'};
    types = {'cellcol', 'cell', 'bcmech', 'bcfluid'};
    
end

