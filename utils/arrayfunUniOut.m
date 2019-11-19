function y = arrayfunUniOut(fun, x)
    y = arrayfun(fun, x, 'UniformOutput', false);
end