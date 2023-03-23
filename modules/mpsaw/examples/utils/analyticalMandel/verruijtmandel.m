%% Verruijt

N = 20;
dx = N/1000;
x = (0 : dx : N)';

b = 1;    % Biot parameter
nu = 0.4; % Poisson coefficient

eta = (1 - nu)/(1 - 2*nu)

fx = tan(x) - 2*eta*x;

plot(x, fx);
axis([0, x(end), -10, 10])
