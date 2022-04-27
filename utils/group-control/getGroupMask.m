function mask = getGroupMask(model, state, names)
    % Get logical mask into wellSol for a given group

    map  = model.getProp(state, 'FacilityWellMapping');
    
    if ~iscell(names), names = {names}; end
    mask = false;
    for i = 1:numel(names)
        mask = mask | arrayfun(@(W) strcmpi(W.group, names{i}), map.W);
    end

end