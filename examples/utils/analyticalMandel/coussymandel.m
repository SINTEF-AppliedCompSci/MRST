N = 20;
dx = N/1000;
x = (dx : dx : N)';

b  = 1;    % biot parameter
iM = 0;    % biot compressibility: iM = 1/M
K  = 1;    % Bulk modulus
nu = 0.4;  % Poisson coefficient

iMKu = K*iM + b^2; % iMKu = Ku/M; page 86 (97) 4.66
B = b/(iMKu);           % page 86 (97) 4.68
a = b*B*(1 - 2*nu)/3;
nuu = (a + nu)/(1 - a); % page 87 (98) 4.71;

eta = (1 - nu)./(nuu - nu)/2 % for comparison with verruijt

fx = tan(x) - (1 - nu)./(nuu - nu)*x;

plot(x, fx)
