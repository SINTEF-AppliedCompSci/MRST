function [p_scaled, p_unscaled] = getParameterVectors(setup, params)
    % Get scaled and unscaled parameter vectors for a problem
    p_scaled = getScaledParameterVector(setup, params);
    p_unscaled = cellfun(@(x) x.value, params, 'UniformOutput', false);
end