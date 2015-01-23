function [eqs, qBC, BCTocellMap, qSRC, srcCells] = addFluxesFromSourcesAndBC(model, eqs, pressure, rho, mob, b, s, forces)
    
    [qBC, qSRC] = deal(cell(numel(mob), 1));
    [BCTocellMap, srcCells] = deal([]);
    
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
        [qSRC, srcCells] = getSourceFluxesAD(model, mob, b, s, forces.src);
        for i = 1:numel(qSRC)
            % Subtract fluxes
            eqs{i}(srcCells)  = eqs{i}(srcCells) - qSRC{i};
        end
    end
end