function [Pb, Pi] = Pt2PbPi(Pt, CO2Properties, waterProperties, h, H, T, theta)
% Under the vertical equilibrium assumption, compute pressure at interface
% and pressure at bottom, given pressure at top. 
%
% SYNOPSIS:
%   function [Pb, Pi] = Pt2PbPi(Pt, CO2Properties, waterProperties, h, H)
%
% PARAMETERS:
%   Pt              - Pressure field at top
%   CO2Properties   - provides the functions 'rho', 'beta' and 'bder' for CO2
%   waterProperties - provides the functions 'rho', 'beta' and 'bder' for water
%   h               - height of CO2 plume
%   H               - total height of aquifer
%   T               - temperature (@@ scalar for now)
%   theta           - aquifer dip angle
%
% RETURNS:
%   Pb - pressure field at bottom
%   Pi - pressure field at interface
%
if ~isscalar(T)
    error(['T needs to be scalar in Pt2PbPi, as CO2props does not yet handle ' ...
           'temperature vectors']);
end



% computing pressure at interface
Pi = Pt + pressureDifference(CO2Properties, h, Pt, T, theta, true);

% computing pressure at bottom
Pb = Pi + pressureDifference(waterProperties, H-h, Pi, T, theta, true);

end



