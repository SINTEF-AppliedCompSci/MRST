mrstModule add static-modeling

% This script demonstrates the use of the function 'estimateSemivariogram1D',
% which uses a set of input samples, presumably generated from a
% one-dimensional, stationary Gaussian process with unknown covariance, and
% estimates the associated semivariogram, i.e. the function describing the
% degree of spatial dependence between values (variance as function of the
% spatial distance between sampled points).  

%% create artificial sample data
% We first generate a set of sample data to test the
% 'estimateSemivariogram1D' function on

num_measurements = 3000; % number of samples to generate
D = 6000 * meter; % total length of the 1D domain where we realize the
                  % Gaussian process.
sigma = 10 * meter; % parameter describing the correllation distance (for a
                    % Gaussian correllation function)
xvals = linspace(0, D, num_measurements);  % Abscissa of the sample points

%% defining correlation function and simulating sample data

% We define a correlation function (correlation between two samples as a
% function of distance between the samples)
corr_fun = @(x) exp(-(x/(sigma/D)).^2); % normalize sigma by distance 

% We generate the set of samples using the 'GaussianProcess1D' with the
% specified correlation function.
measurements = GaussianProcess1D(num_measurements, corr_fun);

%% Estimating semi-variogram

% We now use the 'estimateSemivariogram1D' function to estimate the
% empirical semivariogram based on the samples (which describes the
% covariance of the process - the relationship between semivariogram \gamma
% and covariance function C is:  C = sill - \gamma.
% 'sill' here refers to the height of the semivariogram where it levels off,
% or in other words the variance between points 'infinitely' far apart.

curvepts = 500;  % we evaluate the empirical semivariogram at this many points
max_reldist = 15 * sigma / D; % how far to evaluate the semivariogram
rel_radius = 3; % width of the averaging window used when estimating the semivariogram
[curve, sill, range] = estimateSemivariogram1D({xvals}, {measurements}, curvepts, ...
                                           max_reldist, rel_radius);

%% Plotting results

figure;
subplot(2,1,1); 
plot(xvals, measurements); 
xlabel('x-position'); ylabel('measurement value'); title("measurements");
subplot(2,1,2);
plot(linspace(0, D * max_reldist, curvepts+1), 2 * curve); 
xlabel('distance'); ylabel('mean squared difference'); title('estimated variogram');

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
