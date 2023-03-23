function y = cellfunUniOut(fun, x)
    y = cellfun(fun, x, 'UniformOutput', false);
end