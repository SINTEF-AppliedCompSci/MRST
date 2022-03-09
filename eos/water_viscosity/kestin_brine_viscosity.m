function mu_mPas = kestin_brine_viscosity(p, T, c)
%Compute the brine viscosity defined by Kestin et al. (1981) ref. Tables of
%the dynamic and kinematic viscosity of aqueous NaCl solution temperature
%range 20-150ºC and the pressure range 0.1-35MPa
% 
% Pressure is in MPa
% Temperature in ºC 
% Salt concentration in mol/kg H2O

    % Compute the beta factors
    % compute reduced concentration c_star = c/cs 
    cs = 6.044 + (0.28e-2).*T + (0.36e-4).*(T.^2); 
    c_star = c./cs; 

    % pressure coefficient for water
    beta_w = -1.297 + (0.574e-1).*T - (0.697e-3).*(T.^2) + (0.447e-5).*(T.^3) - (0.105e-7).*(T.^4);

    % excess pressure coefficient at the saturation 
    beta_sE = 0.545 + (0.28e-2).*T - beta_w; 

    % reduced coefficient for pressure
    beta_star = 2.5.*c_star - 2.0.*(c_star.^2) + 0.5.*(c_star.^3);

    % beta factor for viscosity equation 
    beta = beta_sE.*beta_star + beta_w;

    % Compute the different viscosities

    % Compute A and B factors 

    A = (0.3324e-1).*c + (0.3624e-2).*(c.^2) - (0.1879e-3).*(c.^3);
    B = -(0.396e-1).*c + (0.102e-1).*(c.^2) - (0.702e-3).*(c.^3);


    mu_w0_20 = 1002;                                                           % viscosity of mu_w0 at 20ºC in µPa s. 

    % compute µ_w0 (mu_w0) 
    log_mu_w0_ov_mu_w0_20 = (20-T)./(96+T).*[1.2378 - (1.303e-3).*(20-T) + (3.06e-6).*((20-T).^2) ...
                                            + (2.55e-8).*((20-T).^3)];

    mu_w0 = (10.^log_mu_w0_ov_mu_w0_20).*mu_w0_20;

    % compute µ_0 (mu_0)
    log_mu_0_ov_mu_W0 = A + B.*log_mu_w0_ov_mu_w0_20;                                    
    mu_0 = (10.^log_mu_0_ov_mu_W0).*mu_w0;                                  

    % finally compute the viscosity of brine 
    mu = mu_0.*(1+beta.*p.*1e-3);

    mu_mPas = mu;

end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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