mrstModule add static-modeling

% Use maximum likelihood estimation to determine optimal parameters of a 
% probability distribution.  Whereas any probability distribution with any
% type and number of parameters could in principle be used, we here
% demonstrate mle estimation on a very simple example involving a univariate
% normal probability distribution.

% We generate a set of samples from a normal distribution with a specified
% mean (mu) and variance (sigma^2).  The task is to use maximum
% likelihood to estimate these parameters given only the samples.  

%% Define the probability distribution function and its gradient wrt. parameters

% We specify the normal distribution function in 'x', with mean (mu) and
% standard deviation (sigma) specified by the vector p = [mu, sigma]
normal_dist = @(x, p) 1/sqrt(2 * pi * p(2)) * exp( - (x - p(1)).^2/(2*p(2)));

% We then specify the partial derivative of the normal distribution function
% with respect to the mean 'mu' (as specified by p(1))
normal_dist_d_mu = @(x, p) (x-p(1)) .* exp(- (x-p(1)).^2 ./ (2 * p(2))) ./ ...
                    (sqrt(2*pi) .* p(2).^(3./2));

% We also specify the partial derivative of the normal distribution function 
% with respect to the variance 'sigma'^2, (as specified by the square of p(2))
normal_dist_d_sigma2 = @(x, p) ...
    exp( -(p(1) - x).^2./(2 * p(2))) .* ( (p(1)-x).^2 - p(2)) ./ ...
    (2 * sqrt(2*pi) .* p(2).^(5/2));

% Combining the two partial derivative functions above, we can specify
% the 'gradient' with respect to parameters mu and sigma2
normal_dist_grad = @(x, p) [normal_dist_d_mu(x, p), normal_dist_d_sigma2(x, p)];


%% Generate samples from the distribution, for a set of 'true' parameters
mu = 1.5; % this is the 'true' mean (which we will aim to estimate from the
          % samples)
sigma2 = 0.45; % this is the true variance (which we will aim to estimate
               % from the samples)

N = 500; % number of samples

samples = sqrt(sigma2) * randn(N, 1) + mu;

%% Estimating lambda using MLE

% We pretend we do not know the actual distribution sampled.  We assume it is
% a normal distribution (as we specified above with 'normal_dist' and its
% gradient 'normal_dist_grad'), and we aim to estimate the parameters for
% which the available samples would be most likely. 


init_guess = [0, 0.1]; % This is our initial guess of the unknown parameters
                       % [mu_init, sigma2_init];
pmin = [-2, 1e-4]; % We will not search for parameter values below this range 
pmax = [ 2, 4];  % We will not search for parameter values above this range

% Run 'mle_prob_param' to estimate the unknown parameters from the provided samples
[params, hist] = mle_prob_param(samples, normal_dist, normal_dist_grad, ...
                                init_guess, pmin, pmax);
   
%% Validating results

fprintf('Estimated mean:%f, sample mean: %f, actual mean: %f.\n', ...
        params(1), mean(samples), mu);
fprintf('Estimated variance:%f, sample variance: %f, actual variance: %f.\n', ...
        params(2), mean((samples - mean(samples)).^2), sigma2);

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
