function x = incrementSubset(x, subs, v)
% Update a subset directly
    if isa(x, 'ADI')
        x(subs) = x(subs) + v;
    else
        x(subs, :) = x(subs, :) + v;
    end
end