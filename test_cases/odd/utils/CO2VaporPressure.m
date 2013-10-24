function p = CO2VaporPressure(T)
%% Compute the CO2 vapor pressure for a given temperature
% Based on formula 3.13 (p. 1524) of the paper:
% "A new Equation of State for Carbon Dioxide [...]" by Span and Wagner

%% defining necessary constants

[P_c, T_c] = CO2CriticalPoint();
% P_c = 7.3773 * mega * Pascal; % critical pressure
% T_c = 304.1282;               % critical temperature (in Kelvin)

a = [-7.0602087, 1.9391218, -1.6463597, -3.2995634];
t = [1.0, 1.5, 2.0, 4.0];

%% Computing p according to the formula

% computing dimensionless temperature
dlT = T(:)/T_c;

% computing power terms
pterms = bsxfun(@power, 1-dlT, t);

% summing the terms
summed = pterms * a';

% dividing by dimensionless temperature to obtain the right-hand
% side of equation (3.13) in paper
rhs = summed ./ dlT;

% finally computing the resulting pressure
p = P_c * exp(rhs);

%% Ensuring valid range
% Formula only valid up to the critical point.  If temperatures above this
% have been given, replace computed result by NaN.
p(T>T_c) = NaN;

