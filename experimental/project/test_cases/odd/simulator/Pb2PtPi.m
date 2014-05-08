function [Pt, Pi] = Pb2PtPi(Pb, CO2Properties, waterProperties, h, H, T, theta)
% Under the vertical equilibrium assumption, compute pressure at interface
% and top, given pressure at bottom. 
%
% SYNOPSIS:
%   function [Pb, Pi] = Pb2PtPi(Pt, CO2Properties, waterProperties, h, H, T, theta)
%
% PARAMETERS:
%   Pt              - Pressure field at top
%   CO2Properties   - provides the functions 'rho', '
%   waterProperties - provides the functions 'rho', '
%   h               - height of CO2 plume            
%   H               - total height of aquifer        
%   T               - temperature (@@ scalar for now)
%   theta           - aquifer dip angle
%
% RETURNS:
%   Pb - pressure field at bottom   
%   Pi - pressure field at interface

if ~isscalar(T)
    error(['T needs to be scalar in Pt2PbPi, as CO2props does not yet handle ' ...
           'temperature vectors']);
end


gval = norm(gravity) * cos(theta); % scalar value of gravity in play

% computing pressure at interface
Pi = Pb + pressureDifference(waterProperties, -(H-h), Pb, T, theta, false);

% computing pressure at top
Pt = Pi + pressureDifference(CO2Properties, -h, Pi, T, theta, false);

end
