function [bc, changed] = convertPressureBoundaryConditionsToFlux(G, state, bc)
    if isempty(bc)
        changed = false;
    else
        isDir = strcmpi(bc.type, 'pressure');
        changed = any(isDir);
        if any(isDir)
            % Convert Dirichlet boundary conditions to flux
            % boundary conditions for the transport
            dirFace = bc.face(isDir);
            q = sum(state.flux(dirFace, :), 2);
            sgn = 1 - 2*(G.faces.neighbors(dirFace, 2) == 0);
            bc.value(isDir) = sgn.*q;
            [bc.type{isDir}] = deal('resflux');
        end
    end
end