function mu_b_co2 = viscBrineCO2IC2012(T, P, m_nacl, w_co2)
% SYNOPSIS
% mu_b_co2 = viscBrineCO2IC2012(T, P, m_nacl, w_co2)
%
% DESCRIPTION:
% Calculate the dynamic viscosity of a solution of H2O + NaCl (brine) with
% dissolved CO2.
%
% RANGE:
% For pure water + CO2, the model is based on experimental data by Kumagai
% et al. (1998), which is for p up to 400 bar and T up to 50 C, and Bando 
% et al. (2004), which is valid in 30-60C and 10-20MPa.
% The model of Mao & Duan (2009) for brine viscosity reaches 623K, 1000 bar
% and high ionic strength. However, the model used to determine the viscosity 
% when co2 dissolves in the brine (Islam & Carlson, 2012) is based on
% experimental data by Bando et al. (2004) and Fleury and Deschamps (2008),
% who provided experimental data up to P = 200 bar, T extrapolated to
% 100 C, and maximum salinity of 2.7M.
%
% INPUT:
% T:          Double with temperature value in Kelvin
% P:          Double with pressure value in bar
% m_nacl:     Salt molality (NaCl) in mol/kg solvent
% w_co2:      Mass fraction of CO2 in the aqueous solution (i.e. brine)
%
% OUTPUT: 
% mu_b_co2:   Double with dynamic viscosity in Pa*s
%

% Check range
if m_nacl == 0
    if P > 400
        warning(['Pressure out of measured range in ', mfilename])
    end
    if T < 273.15 || T > 333.15
        warning(['Temperature out of tested range in ', mfilename])
    end
else
    if P > 200
        warning(['Pressure out of measured range in ', mfilename])
    end
    if T < 273.15 + 35 || T > 373.15
        warning(['Temperature out of tested range in ', mfilename])
    end
    if m_nacl > 3.1 % m approx, limit of experimental data is 2.738M
        warning(['Salinity of tested range in ', mfilename])
    end
end

% Units
P_mpa = P/10;

% Pure water density (see Islam & Carlson, 2012)
a = 1.34136579*10^2;
b = [-4.07743800*10^3, 1.63192756*10^4, 1.37091355*10^3];
c = [-5.56126409*10^-3, -1.07149234*10^-2, -5.46294495*10^-4];
d = [4.45861703*10^-1, -4.51029739*10^-4];
cp = [1, 2];

rho_h2o = a + sum(b.*10.^(c.*T)) + sum(d.*P_mpa.^cp);                     % [kg/m3]
rho_h2o = rho_h2o/10^3;                                                   % [g/cm^3]

% Pure water viscosity (Mao & Duan, 2009)
d = [0.28853170*10^7, -0.11072577*10^5, -0.90834095*10, ...
     0.30925651*10^-1, -0.27407100*10^-4, -0.19283851*10^7, ...
     0.56216046*10^4, 0.13827250*10^2, -0.47609523*10^-1, ...
     0.35545041*10^-4];
coefs1 = (1:5)-3;
coefs2 = (6:10)-8;

mu_h2o = exp(sum(d(1:5).*T.^coefs1) + sum(d(6:10).*rho_h2o.*T.^coefs2));

if m_nacl == 0
    % Viscosity of H2O + CO2 (Islam & Carlson, 2012)
    % a = [7.632609119*1e2 -9.46077673946*1e3]; % Reported by IC2012
    a = [1.632609119*1e3 -9.46077673946*1e2];   % Actual fit to plots in paper + data by Bando et al. (2004)
    b = [-1.047187396332*1e4 3.68325597*1e1];
    c1 = [1 2];
    c2 = [0 1];

    mu_r = 1 + sum(a.*w_co2.^c1)/sum(b.*T.^c2);
    mu_b_co2 = mu_r*mu_h2o;
else
    % Brine viscosity (H2O + NaCl) ( Mao & Duan, 2009)
    A = -0.21319213 + 0.13651589*10^(-2)*T -0.12191756*10^(-5)*T^2;
    B = 0.69161945*10^-1 - 0.27292263*10^(-3)*T + 0.20852448*10^(-6)*T^2;
    C = -0.25988855*10^-2 + 0.77989227*10^(-5)*T;
    
    mu_b = exp(A*m_nacl + B*m_nacl^2 + C*m_nacl^3) * mu_h2o;
    
    % Brine viscosity with dissolved CO2 (H2O + NaCl + CO2) (Islam & Carlson,
    % 2012)
    mu_b_co2 = mu_b*(1 + 4.65*w_co2^1.0134);
end

% Data from Bando et al., J Chem Eng Data (2004) for test
% xcs = P_mpa/(36.1*P_mpa +3.87*T-1097.1+w_salt*(196*P_mpa + 26.9*T-8810))
% dat = [303.15 100 0.89e-3; 303.15 200 0.88e-3; 313.15 100 0.7e-3; 
%        313.15 200 0.71e-3; 323.15 100 0.57e-3; 323.15 200 0.57e-3;
%        333.15 100 0.47e-3; 333.15 200 0.48e-3];

end