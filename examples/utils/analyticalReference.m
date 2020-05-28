function res = analyticalReference(lambda, mu)
% Compute reference analytical solution

    % Analytical solution
    d1 = @(x, y) x .* (1 - x) .* sin(2 * pi * y); 
    d2 = @(x, y) sin(2 * pi * x) .* sin(2 * pi * y); 
    dvec = @(coord) [d1(coord( :, 1), coord( :, 2)), d2(coord( :, 1), coord( :, 2))]; 
    % divD = diff(d1, x) + diff(d2, y); 
    d1dx = @(x, y) (1 - 2*x).*sin(2*pi*y); 
    d1dxdx = @(x, y) - 2.*sin(2*pi*y); 
    d1dxdy = @(x, y) (1 - 2*x).*2*pi.*cos(2*pi*y); d1dydx = d1dxdy; 
    d2dy = @(x, y) sin(2*pi*x) .*2*pi.*cos(2*pi*y); 
    d2dydy = @(x, y) - sin(2*pi*x) .*(2*pi).^2.*sin(2*pi*y); 
    d2dydx = @(x, y) cos(2*pi*x) .*(2*pi).^2.*cos(2*pi*y); d2dxdy = d2dydx; 
    d1dy = @(x, y) x.*(1 - x).*2*pi.*cos(2*pi*y); 
    d1dydy = @(x, y) - x.*(1 - x).*(2*pi)^2.*sin(2*pi*y); 
    d2dx = @(x, y) cos(2*pi*x).*sin(2*pi*y)*2*pi; 
    d2dxdx = @(x, y) - sin(2*pi*x).*sin(2*pi*y)*(2*pi)^2; 

    divD = @(x, y) d1dx(x, y) + d2dy(x, y); 
    divDdx = @(x, y) d1dxdx(x, y) + d2dydx(x, y); 
    divDdy = @(x, y) d1dxdy(x, y) + d2dydy(x, y); 
    s11 = @(x, y) 2 * mu * d1dx(x, y) + lambda * divD(x, y); 
    s11dx = @(x, y) 2*mu *d1dxdx(x, y) + lambda * divDdx(x, y); 
    s12 = @(x, y)  mu * (d1dy(x, y) + d2dx(x, y)); 
    s12dy = @(x, y) mu * (d1dydy(x, y) + d2dxdy(x, y)); 
    s12dx = @(x, y)  mu * (d1dydx(x, y) + d2dxdx(x, y)); 
    s21 = @(x, y)  mu * (d1dy(x, y) + d2dx(x, y)); 
    s21dx = @(x, y) mu*(d1dydx(x, y) + d2dxdx(x, y)); 
    s21dy = @(x, y)  mu * (d1dydy(x, y) + d2dxdy(x, y)); 
    s22 = @(x, y) 2 * mu * d2dy(x, y) + lambda * divD(x, y); 
    s22dy = @(x, y) 2 * mu * d2dydy(x, y) + lambda * divDdy(x, y); 

    % Apply divergence
    mrhs1 = @(x, y) s11dx(x, y) + s21dy(x, y); 
    mrhs2 = @(x, y) s12dx(x, y) + s22dy(x, y); 

    res = struct('mrhs2' , mrhs2 , ... 
                 'mrhs1' , mrhs1 , ... 
                 'dvec'  , dvec  , ...
                 'd1'    , d1    , ... 
                 'd2'    , d2    , ...
                 's11'   , s11   , ...
                 's21'   , s21   , ...  
                 's12'   , s12   , ... 
                 's22'   , s22);
    
end