function res = soft_ubound(x, upper, delta)

    res = x;

    if isscalar(upper)
        upper = (x+upper)*0 + upper;
    end
    
    % indices affected by upper bound
    aix = x > (upper - delta);
    if sum(aix) == 0
        return;
    end

    % apply an adjustment with the property that derivative is continuous
    res(aix) = upper(aix) - delta * exp(-(x(aix) - (upper(aix) - delta)) / delta);
    
end
