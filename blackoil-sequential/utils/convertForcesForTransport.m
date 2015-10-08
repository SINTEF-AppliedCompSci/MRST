function forces = convertForcesForTransport(state, forces)
    if ~isfield(forces, 'bc') || isempty(forces.bc)
        return
    end
    isDir = strcmpi(forces.bc.type, 'pressure');
    [forces.bc.type{isDir}] = deal('flux');
    
    qBC = state.flux(forces.bc.face(isDir), :);
    qT = sum(qBC, 2);

    forces.bc.value(isDir) = qT;
    forces.bc.sat(isDir, :) = qBC./qT;
    forces.bc.sat(isnan(forces.bc.sat)) = 0;
end