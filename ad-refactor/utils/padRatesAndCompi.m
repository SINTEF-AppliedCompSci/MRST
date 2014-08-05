function [q_new, W, isActive] = padRatesAndCompi(q_s, W, model)
    % Ensure that fluxes exist for all phases, but let the value be zero
    % for phases not present to make rest of implementation cleaner
    isActive = model.getActivePhases();
    dummy = zeros(size(double(q_s{1})));
    
    q_new = cell(1, 3);
    
    [q_new{:}] = deal(dummy);
    q_new(isActive) = q_s;
    
    for i = 1:numel(W)
        W(i).compi = padByIndex(W(i).compi, isActive);
    end
end

function tmp = padByIndex(val, index)
    tmp = zeros(1, 3);
    tmp(index) = val;
end