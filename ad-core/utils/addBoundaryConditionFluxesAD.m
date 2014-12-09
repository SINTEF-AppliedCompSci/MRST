function eqs = addBoundaryConditionFluxesAD(model, eqs, gdz, pressure, rho, mob, b, s, bc)
    if isempty(bc)
        return
    end
    [qRes, bcface2cellMap] = getBoundaryConditionFluxesAD(model, gdz, pressure, rho, mob, b, s, bc);
    
    for i = 1:numel(qRes)
        eqs{i}  = eqs{i} + bcface2cellMap*qRes{i};
    end
end