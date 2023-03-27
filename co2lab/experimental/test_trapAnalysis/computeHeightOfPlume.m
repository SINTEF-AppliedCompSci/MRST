function [h, h_max] = computeHeightOfPlume(Gt, state, sw, sr)

    % computing height corresponding to the current saturation state
    if isfield(state, 'sGmax')
        % We operate with dissolution, implying a gradually receding 'max' saturation
        smax = state.sGmax;
    else
        % There is no dissolution, so the historical maximum is the current
        % maximum
        smax = state.smax(:,2);
    end

    [h, h_max] = upscaledSat2height(state.s(:,2), smax, Gt, 'resSat', [sw sr]);
end
