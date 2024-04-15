function rho_b = pvtBrineBW1992(T, P, w_nacl)
% SYNOPSIS:
% rho_b = pvtBrineBW1992(T, P, w_nacl)
%
% DESCRIPTION:
% Calculate the brine (H2O + NaCl) density based on Batzle & Wang (1992).
% These authors used the data of Rowe and Chou (1970), Zarembo & Fedorov
% (1975) and Potter & Brown (1977) to expand the P, T validity range. 
%
% RANGE:
% P valid from 5 to 100 MPa, T from 20 to 350 C (Adams & Bachu, 2002)
%
% INPUT:
% T:          Double with temperature value in Kelvin
% P:          Double with pressure value in bar
% w_nacl:     Salt (NaCl) mass fraction
%
% OUTPUT: 
% rho_b:      Double with brine density in kg/m3
%

% Check range
if P < 50 || P > 1000
    warning(['Pressure out of measured range in ', mfilename])
end
if T < 273.15 + 20 || T > 273.15 + 350
    warning(['Temperature out of tested range in ', mfilename])
end
if w_nacl > 0.3 % approx, limit is 320,000 mg/l
    warning(['Salinity of tested range in ', mfilename])
end

P = 0.1*P;
T = T-273.15;

% Water density
rho_w = 1 + 10^(-6)*(-80*T -3.3*T^2 +0.00175*T^3 +489*P ...
                     -2*T*P +0.016*T^2*P -1.3*10^(-5)*T^3*P  ...
                     -0.333*P^2 -0.002*T*P^2);                            % [g/cm^3]

% Brine density
fact1 = 300*P - 2400*P*w_nacl + T*(80 +3*T -3300*w_nacl -13*P +47*P*w_nacl);
rho_b = (rho_w + w_nacl*(0.668 + 0.44*w_nacl + 10^(-6)*fact1))*1000;      % [kg/m^3]


end
