function [ rho_b ] = density_brine(p,t,c)
%Compute density of a brine following Spivey et al., 2004
%
% SYNOPSIS:
%   rho_b = density_brine(p,t,c)
% 
% PARAMETERS: 
%   p   - pressure. One value per cell. 
%   t   - temperature. One value per cell. 
%   c   - NaCl mass fraction. One value per cell.
% 
% RETURNS:
%   rho_b - brine density. One value per cell 
% 
% SEE ALSO:
%   'density_pure_water', 'viscosity_pure_water', 'viscosity_brine'

%{
Copyright 2009-2022 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    t               = t-273.15;        % conversion from Kelvin to �Celsius
    MW_NaCl         = 0.0584428;       % [kg/mol] Molar weigth of NaCl
    X_NaCl          = c;               % mass fraction of NaCl  
    X_Water         = 1 - X_NaCl;      % mass fraction of water    
    cNaCl           = X_NaCl./X_Water; % concentration of NaCl per kg of water [kg/kg]   
    mol_NaCl        = cNaCl./MW_NaCl;  % molality [mol/kg H2O]  
    p0              = 70e6;            % [Pa] Reference pressure

    tMin            = 0;
    tMax            = 275;
    pMin            = 0;
    pMax            = 200e6;
    mol_NaClMin     = 0;
    mol_NaClMax     = 5.7;

    if max(value(t)) > 275 || min(value(t)) < 0
        warning('T out of range')
    end 

    if max(value(p)) > 200e6
        warning('p out of range')
    end

    if max(value(mol_NaCl)) > 5.7
        warning('mol_NaCl out of range')
    end

    % Caping for the ADI variable scheme
    t        = min(max(t, tMin), tMax);
    p        = min(max(p, pMin), pMax);
    mol_NaCl = min(max(mol_NaCl, mol_NaClMin), mol_NaClMax);

    % Get coefficients for pure water to compute density brine (rho_b0) and
    % coefficients for the brine (Eb, Fb)
    [rho_w0, Ew, Fw] = coefficients_pure_water( t );

    % Compute the density of the brine at the reference pressure
    rho_b0 = density_brine_refP(rho_w0,t,mol_NaCl);

    % Compute the coefficients Eb and Fb for the brine
    [Eb,Fb] = coefficients_brine(Ew,Fw,mol_NaCl,t);

    % Compute isothermal compressibility
     Ib  = isothermal_compressibility(p,p0,Eb,Fb);
     Ib0 = isothermal_compressibility(p0,p0,Eb,Fb);

    % Compute density
     rho_b = rho_b0 .* exp(Ib - Ib0); % Eq.(15) in Spivey et al., 2004
     rho_b = rho_b.*1000;

end

%--------------------------------------------------------------------------
%                             References
%--------------------------------------------------------------------------
% Spivey, J.P., Mccain, W.D., North, R. 2004. Estimating density, formation
% volume factor, compressibility, methane solubility, and viscosity for
% oilfield brines at temperatures from 0 to 275ºC, pressures to 200 Mpa,
% and salinities to 5.7 mol.kg-1, Journal of Canadian Petroleum Technology,
% 43 (7), 52-60.