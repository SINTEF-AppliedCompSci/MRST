function mu = viscCO2F1998(T, rho)
% SYNOPSIS:
% mu = viscCO2F1998(T, rho)
%
% DESCRIPTION:
% Calculate CO2 viscosity from Vesovic et al., J Phys Chem Ref Data (1990) 
% and Fenghour et al., J Phys Chem Ref Data (1998), as described in
% Hassanzadeh et al., IJGGC (2008). 
%
% RANGE:
%
% INPUT:
% T:    Double with temperature value in Kelvin
% rho:  Double with density value in kg/m^3
%
% OUTPUT: 
% mu:   Double with dynamic viscosity in Pa*s
%

% Check if T, P conditions are within range

% From Hassanzadeh et al., IJGGC (2008)
% Constants
e = [0.235156, -0.491266, 5.211155*10^-2, ...
     5.347906*1e-2, -1.537102*1e-2];
f = [5.5934*1e-3, 6.1757*1e-5, 0.0, 2.6430*1e-11];
g = [0.4071119*1e-2, 0.7198037*1e-4, 0.2411697*1e-16, ...
     0.2971072*1e-22, -0.1627888*1e-22];
exps = [0 1 2 3 4];

% Functions of scaled T
T_x    = T*(1/251.196);
psi_m  = exp(sum(e.*log(T_x).^exps));

% Compute viscosity
a1 = 1.00697*T^0.5/psi_m;
a2 = g(1)*rho;
a3 = g(2)*rho^2;
a4 = g(3)*rho^6/T_x^3;
a5 = g(4)*rho^8;
a6 = g(5)*rho^8/T_x;
if T < 313.15
    a7 = sum(f.*rho);
else
    a7 = 0;
end
mu = (a1 + a2 + a3 + a4 + a5 + a6 + a7)/10^6;

end