function forces = convertForcesForTransport(state, forces)
    if ~isfield(forces, 'bc') || isempty(forces.bc)
        return
    end
    isDir = strcmpi(forces.bc.type, 'pressure');
    [forces.bc.type{isDir}] = deal('flux');
    
    qBC = state.flux(forces.bc.face(isDir), :);
    qT = sum(qBC, 2);
    netq = sum(abs(qBC), 2);
    
    sat0 = forces.bc.sat;
    forces.bc.value(isDir) = qT;
    forces.bc.sat(isDir, :) = bsxfun(@rdivide, abs(qBC), netq);

    bad = any(isnan(forces.bc.sat), 2);
    forces.bc.sat(bad, :) = sat0(bad, :);
end