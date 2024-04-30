function viscMixture = viscGasMixtD1993(x, M, mu)
% SYNOPSIS:
% viscMixture = viscGasMixtD1993(x, M, mu)
%
% DESCRIPTION:
% Calculate the viscosity of a gas mixture following Davidson (1993).

% RANGE:
% In principle, valid for any range within which the individual components'
% viscosities are valid.
%
% INPUT:
% x:          Mole fraction of each component
% M:          Molar mass of each component
% mu:         Viscosity of each component in user chosen units
% 
% Each input should be a double array of size (1, n) where n is the total
% number of components
%
% OUTPUT: 
% viscMixture:    Double with viscosity of the mixture in same units as mu
%

% Momentum fractions
num = x.*sqrt(M);
den = sum(num);
y   = num/den;

% Mass fraction
M_mixt  = sum(x.*M);
w = x.*(M./M_mixt);

% Efficiency of momentum transfer
E = 2.*sqrt(w)'.*sqrt(w)./(w'+w); 
E(isnan(E)) = 0;                   % fix for when some component is absent.

% Fluidity
A = 0.375;  % Best fit (N = 164 mixtures) with RMS = 1.28% (Davidson, 1993)
f = sum(sum(y'.*y./(sqrt(mu)'.*sqrt(mu)).*E.^A));

% Viscosity
viscMixture = 1/f;

% (Slow way of computing E and F)
% N = numel(x);
% f = 0;
% for i = 1:N
%     for j = 1:N
%         E = 2*sqrt(w(i))*sqrt(w(j))/(w(i) + w(j));
%         f = f + (y(i)*y(j)/(sqrt(mu(i))*sqrt(mu(j)))*E^A);
%     end
% end
end