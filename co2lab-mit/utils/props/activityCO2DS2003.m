function gamma_co2 = activityCO2DS2003(T, P, m_io)
% SYNOPSIS:
% gamma_co2 = activityCO2DS2003(T, P, m_io)
%
% DESCRIPTION:
% Calculate a CO2 pseudo activity coefficient based on a virial expansion
% of excess Gibbs energy.
%
% RANGE:
%
% INPUT:
% T:    Double with temperature value in Kelvin
% P:    Double with pressure value in bar
% m_io: Double array of size(1x6) where each entry corresponds to the
%       molality of a particular ion in the initial brine solution. The
%       order is as follows:
%       [ Na(+),   K(+),  Ca(2+), Mg(2+), Cl(-), SO4(2-)]
%
% OUTPUT: 
% V_m:          Double with molar volume in [cm^3/mol]
% rhox:         Double with density in [mol/m^3]
% rho:          Double with density in [kg/m^3]
%

% Check if T, P conditions are within range
if P > 2000
    warning(['Pressure out of tested range in ', mfilename])
end
if T < 273 || T > 533
    warning(['Temperature out of tested range in ', mfilename])
end
if any(m_io > 4.3)
    warning(['Ionic strength out of tested range in ', mfilename])
end

% Interaction parameter constants
c1_co2_C = -0.411370585; c2_co2_C = 6.07632013*10^(-4);
c3_co2_C = 97.5347708;   c4_co2_C = -0.0237622469;
c5_co2_C = 0.0170656236; c6_co2_C = 1.41335834*10^(-5);

c1_co2_C_A = 3.36389723*10^(-4);  c2_co2_C_A = -1.98298980*10^(-5);
c3_co2_C_A = 0;                   c4_co2_C_A = 2.122220830*10^(-3);
c5_co2_C_A = -5.24873303*10^(-3); c6_co2_C_A = 0;

% Compute interaction parameters
lam_co2_C    = c1_co2_C + c2_co2_C*T + c3_co2_C/T + ...
    c4_co2_C*P/T + c5_co2_C*P/(630-T) + ...
    c6_co2_C*T*log(P);

zet_co2_C_A = c1_co2_C_A + c2_co2_C_A*T + c3_co2_C_A/T + ...
    c4_co2_C_A*P/T + c5_co2_C_A*P/(630-T) + ...
    c6_co2_C_A*T*log(P);

% Compute activity coefficient
gamma_co2  = exp(2*lam_co2_C*(sum(m_io(1:2)) + 2*sum(m_io(3:4))) + ...
                 zet_co2_C_A*m_io(5)*(sum(m_io(1:4))) - 0.07*m_io(6));

end