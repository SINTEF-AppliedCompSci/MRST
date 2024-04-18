function [rho_brine, c_brine] = pvtBrineRoweChou1970(T, P, S)
% SYNOPSIS:
% [rho_brine, c_brine] = pvtBrineRoweChou1970(T, P, S)
%
% DESCRIPTION:
% Calculate brine density and/or compressibility using Rowe and Chou's
% (1970) correlation.
%
% RANGE:
% P <= 35 MPa
% 293.15 <= T <= 423.15 K
%
% INPUT:
% T:    Double with temperature value in Kelvin
% P:    Double with pressure value in bar
% S:    Salt mass fraction
%
% OUTPUT: 
% rho_brine:    Double with density value in kg/m3
% c_brine:      Double with compressibility value in 1/kPa
%

% Check if T, P conditions are within range
if P > 350
    warning(['Pressure out of range in ', mfilename])
end
if T < 293.15 || T > 423.15
    warning(['Temperature out of range in ', mfilename])
end

% Convert pressure to kg/cm^2
bar_to_kgfcm2 = 1.01972;
q             = P*bar_to_kgfcm2;                        % [kgf/cm^2]
p_kpa         = P*100;                                  % [kPa]

% Correlation factors
a1 = 5.916365 - 0.01035794*T + 0.9270048*10^(-5)*T^2 - 1127.522/T ...
     + 100674.1/T^2;
a2 = 0.520491*10^(-2) - 0.10482101*10^(-4)*T + 0.8328532*10^(-8)*T^2 ...
     - 1.1702939/T + 102.2783/T^2;
a3 = 0.118547*10^(-7) - 0.6599143*10^(-10)*T;
a4 = -2.5166 + 0.0111766*T - 0.170522*10^(-4)*T^2;
a5 = 2.84851 - 0.0154305*T + 0.223982*10^(-4)*T^2;
a6 = -0.0014814 + 0.82969*10^(-5)*T - 0.12469*10^(-7)*T^2;
a7 = 0.0027141 - 0.15391*10^(-4)*T + 0.22655*10^(-7)*T^2;
a8 = 0.62158*10^(-6) - 0.40075*10^(-8)*T + 0.65972*10^(-11)*T^2;

% Compute density and compressibility
rho_brine = 1/((a1 - q*a2 - q^2*a3 + a4*S + a5*S^2 - q*a6*S - q*a7*S^2 ...
               -0.5*q^2*a8*S)*0.001);                  % [kg/m^3]

p_ref     = 101.325;                                   % [kPa]
r         = (p_ref*0.01)*bar_to_kgfcm2;                % [kgf/cm^2]
rho_ref   = 1/((a1 - r*a2 - r^2*a3 + a4*S + a5*S^2 - r*a6*S - r*a7*S^2 ...
               -0.5*r^2*a8*S)*0.001);                  % [kgf/m^3]
c_brine   = (rho_brine - rho_ref)/(rho_brine*(p_kpa - p_ref)); % [1/kPa]

end