mrstModule add static-modeling

% In this example, we generate some realizations of 1D and 2D Gaussian random field

%% Generate 1D gaussian random field

% In order to simulate a stationary Gaussian process, we need to specify the
% covariance function, which descibes the covariance between points as a
% function of their mutual distance.  A typical example is the use of an
% exponential correlation function, which we specify as a lambda function below.
num_samples = 1000;
corr_fun = @(x) exp(-abs(x)/0.3); % we use exponential correlation function

% We use the correlation function and number of samples as arguments to the
% 'GaussianProcess1D' function to produce a realization of the 1D random field.
field_1D_A = GaussianProcess1D(num_samples, corr_fun);

% We repeat the process using a different covariance function, namely the
% Gaussian correlation function.
corr_fun = @(x) exp(-(x/0.05).^2);
field_1D_B = GaussianProcess1D(num_samples, corr_fun);

% Plot the results
clf;
subplot(1,2,1); plot(field_1D_A); title('Exponential correlation');
subplot(1,2,2); plot(field_1D_B); title('Gaussian correlation');

%% Generate 2D gaussian random field

% We now generate some 2D realizations

num_samples = 400;

% generating field using a Gaussian covariance function.  The main
% difference from the 1D case is that the correlation function now takes a
% 2-component vector. Only the norm of the vector will be used by the
% correlation function, though.
corr_fun = @(xy) exp(-sum(xy.^2, 2)/0.0001);
fieldA = GaussianProcessND([num_samples, num_samples], corr_fun);

% Re-generating field, changing correlation length but keeping the shape of
% the correlation function (Gaussian)
corr_fun = @(xy) exp(-sum(xy.^2, 2)/0.01);
fieldB = GaussianProcessND([num_samples, num_samples], corr_fun);

% Re-generating field using exponential correlation function
corr_fun = @(xy) exp(-sqrt(sum(xy.^2, 2))/0.1);
fieldC = GaussianProcessND([num_samples, num_samples], corr_fun);

% plotting the results and comparing

figure;
subplot(1,3,1); surf(fieldA, 'edgealpha', 0); view(0, 90); 
title("Gaussian correllation, sigma=1e-2");
                                                  
subplot(1,3,2); surf(fieldB, 'edgealpha', 0); view(0, 90);
title("Gaussian correllation, sigma=1e-1");

subplot(1,3,3); surf(fieldC, 'edgealpha', 0); view(0, 90);
title("Exponential correllation");

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
