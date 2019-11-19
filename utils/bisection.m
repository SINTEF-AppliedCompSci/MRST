function [x, fx] = bisection(f, bot, top, tol)
% Find root by bisection method
% f: function handle
% bot, top: inintial X guess
    if f(bot) * f(top)>0
        error('Invalid boundary')
    end
    while abs(top - bot) > tol
        x  = (bot + top)/2;
        fx = f(x);
        if fx == 0
            bot = x;
            top = x;
        elseif f(bot) * fx < 0
            top = x ;
        else
            bot = x;
        end
    end
end