function [E, K, G, nu] = HS_bound(E_m, E_f, nu_m, nu_f, vol_m, vol_f, type)
%
% SYNOPSIS:
%   function [E, K, G, nu] =  HS_bound(E_m, E_f, nu_m, nu_f, vol_m, vol_f, type)
%
% DESCRIPTION:
%   Function to calculate the upscaled properties of a two phase, isotropic
%   material (in this case a fractured rock mass), according to the Hashin-
%   Shtrikman bounds (Hashin and Shtrikman 1962).
%
% PARAMETERS:
%   E_m     - matrix Young's modulus
%   E_f     - fracture Young's modulus
%   nu_m    - matrix Poisson's ratio
%   nu_f    - fracture Poisson's ratio
%   vol_m   - matrix volume fraction
%   vol_f   - fracture volume fraction
%   type    - upper or lower bound
%
% RETURNS:
%   E    - upscaled Young's modulus (2D)
%   K    - upscaled bulk modulus
%   G    - upscaled shear modulus
%   nu   - upscaled Poisson's modulus (2D)
%
% SEE ALSO: example_void_fractures
%

% Shear moduli
G_f = E_f/(2*(1+nu_f));
G_m = E_m/(2*(1+nu_m));

% Bulk moduli, with plane-stress
K_f = E_f/(3*(1-2*nu_f)); 
K_m = E_m/(3*(1-2*nu_m)); 

% Calculate 
if strcmp(type, 'lower')
    K = K_f + vol_m/((K_m-K_f)^-1 + vol_f*(K_f + (4/3)*G_f)^-1);
    G = G_f + vol_m/((G_m-G_f)^-1 + (2*vol_f*(K_f + 2*G_f))*(5*G_f*(K_f + (4/3)*G_f))^-1);
elseif strcmp(type, 'upper')
    K = K_m + vol_f/((K_f-K_m)^-1 + vol_m*(K_m + (4/3)*G_m)^-1);
    G = G_m + vol_f/((G_f-G_m)^-1 + (2*vol_m*(K_m + 2*G_m))*(5*G_m*(K_m + (4/3)*G_m))^-1);
end

K_2d = 9*K*G/(3*K+4*G);
nu = (K_2d - G)/(K_2d + G);
E = K_2d*2*(1-nu);
end
