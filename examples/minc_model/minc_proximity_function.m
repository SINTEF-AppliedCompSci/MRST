function [ p, dp ] = minc_proximity_function( x, L )
%PROXIMITY Summary of this function goes here
%   Detailed explanation goes here

    min_x = min(L);
    
    dimension = length(L);
    switch dimension
        case 1
            u = 2*x/L(1);
            if x < min_x
                p = u;
            else
                p = 1;
            end
            dp = 2/L(1);
        case 2
            u = 2*x/L(1);
            v = 2*x/L(2);
            if x < min_x
                p = u + v - u*v;
            else
                p = 1;
            end
            dp = (-8*x + 2*L(1) + 2*L(2))/(L(1)*L(2));
        case 3
            u = 2*x/L(1);
            v = 2*x/L(2);
            w = 2*x/L(3);
            if x < min_x
                p = u + v + w - u*v - u*w - v*w + u*v*w;
            else
                p = 1;
            end
            dp = (24*x^2 - 8*(L(1)+L(2)+L(3))*x +...
                2*(L(1)*L(2) + L(1)*L(3) + L(2)*L(3)))...
                /(L(1)*L(2)*L(2));
        otherwise
            error('Too many entries for fracture spacing L');
    end


end

