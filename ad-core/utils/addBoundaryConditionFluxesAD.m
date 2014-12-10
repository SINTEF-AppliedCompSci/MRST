function eqs = addBoundaryConditionFluxesAD(model, eqs, gdz, pressure, rho, mob, b, s, bc)
    if isempty(bc)
        return
    end
    % Setup the fluxes from the boundary condition
    [qRes, BCTocellMap] = getBoundaryConditionFluxesAD(model, gdz, pressure, rho, mob, b, s, bc);
    
    for i = 1:numel(qRes)
        % Subtract fluxes
        eqs{i}  = eqs{i} - BCTocellMap*qRes{i};
    end
end