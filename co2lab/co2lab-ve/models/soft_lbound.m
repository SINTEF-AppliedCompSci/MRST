function res = soft_lbound(x, lower, delta)

    res = x;

    if isscalar(lower)
        lower = (x+lower)*0 + lower;
    end
    
    % indices affected by lower bound
    aix = x < (lower + delta);
    if sum(aix) == 0
        return;
    end

    % apply an adjustment with the property that derivative is continuous
    res(aix) = lower(aix) + delta * exp((x(aix) - (lower(aix) + delta)) / delta);
    
end
