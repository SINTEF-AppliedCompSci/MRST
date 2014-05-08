function h = columnHeight(Pcap, Tcap, Tgrad, CO2props, rhoWater, theta, mass, col_area, poro, quick)
% 
% Compute the column height necessary to contain 'mass' amount of CO2.  It
% is assumed that the reference surface \zeta_R is that of the CO2/water interface.
% SYNOPSIS:
%   function h = columnHeight(Pcap, Tcap, Tgrad, CO2Props)
%
% PARAMETERS:
%   Pcap       - hydrostatic pressure at CAPROCK level 
%   Tcap       - temperature at CAPROCK level
%   Tgrad      - global temperature gradient (in BASE UNITS, deg/m)
%   CO2props   - CO2 properties object
%   rhoWater   - density of brine (considered constant)
%   theta      - inclination of coordinate system
%   mass       - desired column mass
%   col_area   - column base area
%   poro       - column rock porosity (considered constant in height)    
%   quick      - if 'true', use Taylor expansion rather than solving the ODE    
%
% RETURNS:
%   h - height of the column
%
% EXAMPLE:
%
% SEE ALSO:
%

    % Determine mass per area, taking into account porosity
    mass1D = mass / col_area / poro;
    
    g_cos_t = norm(gravity) * cos(theta);
    G_cos_t = Tgrad * cos(theta);
    
    % defining functions to give reference pressure and temperature for a
    % given column height
    pressFun = @(x) Pcap + x * g_cos_t * rhoWater;
    tempFun  = @(x) Tcap + x * G_cos_t;
    
    % start value for height
    h0  = mass1D / CO2props.rho(Pcap, Tcap);

    % computing height
    if quick % approximate using Taylor expansion
        h = fzero(@(x) compute_column_mass_Taylor(pressFun(x), tempFun(x), x, CO2props, G_cos_t, g_cos_t) - mass1D, h0);
    else % exact (and slow) computation using ODE
        h = fzero(@(x) compute_column_mass_ODE(pressFun(x), tempFun, x, CO2props.rho, g_cos_t) - mass1D, h0);
    end
    
    
            
end

%% Helper functions

function mass = compute_column_mass_Taylor(Pref, Tref, h, CO2props, Gct, gct)
    Ieta = etaIntegrals(CO2props, Pref, Tref, Gct, gct);
    mass = CO2props.rho(Pref, Tref) * Ieta(-h) * h;
end

function mass = compute_column_mass_ODE(Pref, tfun, h, rhofun, g_cos_t)
    Ptop = numIntTopPress(Pref, tfun, h, rhofun, g_cos_t);
    mass = (Pref - Ptop) / g_cos_t;
end

function Ptop = numIntTopPress(Pref, tfun, h, rhofun, g_cos_t)

    res = ode23(@(z, p) g_cos_t * rhofun(p, tfun(z)), [h 0], Pref);
    Ptop = res.y(end);
end




% Pb = ode23(@(z, p) norm(gravity) * rho_fun(p, Tfun(z)), ...
%            [zref(:), zref(:)+h(:)], ...
%            Pref(:));
