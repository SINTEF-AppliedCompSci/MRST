function x = assignValue(x, v, inx)
    if isa(x, 'ADI')
        x.val(inx) = v;
    else
        x(inx)=v;
    end
end