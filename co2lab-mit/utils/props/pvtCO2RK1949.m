function [V_m, rhox, rho] = pvtCO2RK1949(T, P, varargin)
% SYNOPSIS
% [V_m, rhox, rho] = pvtCO2RK1949(T, P)
% [V_m, rhox, rho] = pvtCO2RK1949(T, P, 'a_m', a_m, 'b_m', b_m)
%
% DESCRIPTION:
% Calculate CO2 molar volume and density using Redlich and Kwong (1949)
% EoS (= RK EoS).
%
% RANGE:
% Tested by Spycher et al. (2003) with constant intermolecular attraction
% parameters to yield accurate results in the T range ~10 to ~100C and P 
% range up to 600 bar, for (1) CO2 compressibility factor, (2) CO2 fugacity
% coefficient and (3) mutual solubilities of H2O and CO2 in the gas and 
% aqueous phase (respectively).
%
% INPUT:
% T:    Double with temperature value in Kelvin
% P:    Double with pressure value in bar
%
% Optional arguments:
% a_m:  Intermolecular attraction constant (of the mixture) in bar*cm^6*K^0.5/mol^2
% b_m:  Intermolecular repulsion constant (of the mixture) in cm^3/mol
%
% OUTPUT: 
% V_m:          Double with molar volume in [cm^3/mol]
% rhox:         Double with density in [mol/m^3]
% rho:          Double with density in [kg/m^3]
%

% Check if T, P conditions are within range
if P > 600
    warning(['Pressure out of tested range in ', mfilename])
end
if T < 283.15 || T > 373.15
    warning(['Temperature out of tested range in ', mfilename])
end

% Ideal gas constant
R           = 83.1447;                                  % bar*cm^3/mol*K

if nargin < 4
    % These are for CO2 after Spycher et al. (2003). Infinite dilution of
    % H2O in the gaseous phase is assumed (with y_h2o =0 and y_co2 = 1 in
    % the mixing rules).
   a_m       = 7.54*10^7 - 4.13*10^4*T;                 % bar*cm^6*K^0.5/mol^2
   b_m       = 27.8;                                    % cm^3/mol 
end

% Molar volume 
%   Make sure not to leave gaps between signs and each coefficient
%   when using the "roots" fcn. Also, more than one value below the 
%   critical point is returned. The gas phase volume is always given by the 
%   maximum root (see Appx B.2, Spycher et al., 2003).
c3 = 1;
c2 = (R*T/P);
c1 = ((R*T*b_m)/P - a_m/(P*T^0.5) + b_m^2);
c0 = a_m*b_m/(P*T^0.5);
r = roots([c3 -c2 -c1 -c0]);
V_m = r(abs(imag(r))<eps);                              % [cm^3/mol]
if numel(V_m) > 1   
    V_gas   = max(V_m);                                 % V gas phase  
    V_liq   = min(V_m);                                 % V of liquid phase
    w1      = P*(V_gas - V_liq);
    w21     = R*T*log((V_gas-b_m)/(V_liq-b_m));
    w22     = a_m/(b_m*T^0.5) * log(((V_gas+b_m)*V_liq)/((V_liq+b_m)*V_gas));
    w2 = w21 + w22;
    if w2-w1 > 0+eps,       V_m = V_gas;
    elseif w2-w1 < 0-eps,   V_m = V_liq;
    else,                   V_m = V_gas;
    end                                                   
end

% Molar density
rhox = (1/V_m)*10^6;                                    % [mol/m^3]

% Mass density
Mw_co2 = 44.0095/10^3;                                  % [kg/mol]
rho    = rhox*Mw_co2;                                   % [kg/m^3]

end