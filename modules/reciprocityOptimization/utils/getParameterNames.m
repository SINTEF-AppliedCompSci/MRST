function names = getParameterNames(params_cell)
    % Extracts all parameter names from a cell array of ModelParameter objects
    names = cellfun(@(x) x.name, params_cell, 'UniformOutput', false);
end