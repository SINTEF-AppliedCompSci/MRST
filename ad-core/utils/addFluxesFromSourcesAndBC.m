function eqs = addFluxesFromSourcesAndBC(model, eqs, pressure, rho, mob, b, s, forces)
    if isempty(forces.bc) && isempty(forces.src)
        return
    end
    
    if ~isempty(forces.bc)
        % Setup the fluxes from the boundary condition
        [qBC, BCTocellMap] = getBoundaryConditionFluxesAD(model, pressure, rho, mob, b, s, forces.bc);

        for i = 1:numel(qBC)
            % Subtract fluxes
            eqs{i}  = eqs{i} - BCTocellMap*qBC{i};
        end
    end
    
    if ~isempty(forces.src)
        % Fluxes from source terms
        [qSRC, cells] = getSourceFluxesAD(model, mob, b, s, forces.src);
        for i = 1:numel(qSRC)
            % Subtract fluxes
            eqs{i}(cells)  = eqs{i}(cells) - qSRC{i};
        end
    end
end