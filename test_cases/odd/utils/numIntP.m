function Pb = numIntP(zref, Pref, h, rho_fun, Tfun)
% Given a reference pressure (Pref) and height (zref),
% compute the hydrostatic pressure at 'zref'+'h', given
% the fluid density function 'rho_fun', which is a function of 
% pressure and temperature, and the temperature as a function of depth ('Tfun')

    Pb = ode23(@(z, p) norm(gravity) * rho_fun(p, Tfun(z)), ...
               [zref(:), zref(:)+h(:)], ...
               Pref(:));
    % Pb = ode113(@(z, p) norm(gravity) * rho_fun(p, Tfun(z)), ...
    %            [zref(:), zref(:)+h(:)], ...
    %            Pref(:));

end

    